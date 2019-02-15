#!/usr/bin/env nextflow

//For pretty-printing nested maps etc
import static groovy.json.JsonOutput.*


//RETURNS DNA ALIGNER NAMES/LABELS IF BOTH INDEXING AND ALIGNMENT TEMPLATES PRESENT
Channel.fromFilePairs("${workflow.projectDir}/templates/{index,dna}/*_{index,align}.sh", maxDepth: 1, checkIfExists: true)
  .filter{ params.alignersDNA == 'all' || it[0].matches(params.alignersDNA) }
  .map { [it[0], "DNA"] }
  .set {alignersDNA}

//RETURNS RNA ALIGNER NAMES/LABELS IF BOTH INDEXING AND ALIGNMENT TEMPLATES PRESENT
Channel.fromFilePairs("${workflow.projectDir}/templates/{index,rna}/*_{index,align}.sh", maxDepth: 1, checkIfExists: true)
  .filter{ params.alignersRNA == 'all' || it[0].matches(params.alignersRNA) }
  .map { [it[0], "RNA"] }
  .set { alignersRNA }

//Pre-computed BEERS datasets (RNA)
datasetsSimulatedRNA = Channel.from(['human_t1r1','human_t1r2','human_t1r3','human_t2r1','human_t2r2','human_t2r3','human_t3r1','human_t3r2','human_t3r3'])
  .filter{ !params.debug || it == params.debugDataset }
  .filter{ (it[-1] as Integer) <= params.replicates}

//Download reference: hg19
url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz'


/*
  Generic method for extracting a string tag or a file basename from a metadata map
 */
 def getTagFromMeta(meta, delim = '_') {
  return meta.species+delim+meta.version //+(trialLines == null ? "" : delim+trialLines+delim+"trialLines")
}

def helpMessage() {
  log.info"""
  Usage:

  nextflow run csiro-crop-informatics/biokanga-manuscript -profile singularity

  Default params:
  """.stripIndent()
  println(prettyPrint(toJson(params)))
  // println(prettyPrint(toJson(config)))
  // println(prettyPrint(toJson(config.process)))
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}


process downloadReferenceRNA {
  storeDir {executor == 'awsbatch' ? "${params.outdir}/downloaded" : "downloaded"}
  scratch false
  // scratch true

  input:
    val(url)

  output:
    file('chromFa.tar.gz') into downloadedRefsRNA

  when:
    'simulatedRNA'.matches(params.mode) || 'realRNA'.matches(params.mode)

  script:
  """
  wget ${url}
  """
}

process downloadDatasetsRNA {
  // label 'download'
  tag("${dataset}")
  storeDir {executor == 'awsbatch' ? "${params.outdir}/downloaded" : "downloaded"} // storeDir "${workflow.workDir}/downloaded" put the datasets there and prevent generating cost to dataset creators through repeated downloads on re-runs
  scratch false

  input:
    val(dataset) from datasetsSimulatedRNA

  output:
    file("${dataset}.tar.bz2") into downloadedDatasets
    // set val("${dataset}"), file("${dataset}.tar.bz2") into downloadedDatasets  //not possible when using storeDir

  when:
    'simulatedRNA'.matches(params.mode)

  script:
    """
    wget http://bp1.s3.amazonaws.com/${dataset}.tar.bz2
    """
}

process extractDatasetsRNA {
  label 'slow'
  tag("${dataset}")
  echo true

  input:
    // set val(dataset), file("${dataset}.tar.bz2") from downloadedDatasets
    file(datasetfile) from downloadedDatasets

  output:
    set val(dataset), file("${dataset}") into extractedDatasets

  script:
    dataset = datasetfile.name.replaceAll("\\..*","")
    """
    mkdir -p ${dataset} \
    && pbzip2 --decompress --stdout -p${task.cpus} ${datasetfile} | tar -x --directory ${dataset}
    """
}

process convertReferenceRNA {
  label 'slow'

  input:
    file(downloadedRef) from downloadedRefsRNA

  output:
    set val(meta), file(ref) into refsRNA
    // file('ucsc.hg19.fa') into reference

  script:
  meta = [:]
  meta.mode = 'RNA'
  if(params.debug) {
    ref="${params.debugChromosome}.fa"
    """
    tar xzvf ${downloadedRef} ${ref}
    """
  } else {
    ref='ucsc.hg19.fa'
    """
    tar xzvf ${downloadedRef}
    cat \$(ls | grep -E 'chr([0-9]{1,2}|X|Y)\\.fa' | sort -V)  > ${ref}
    """
  }
}

/*
 1. Input pointers to FASTA converted to files, NF would fetch remote as well and create tmp files,
    but avoiding that as may not scale with large genomes, prefer to do in process.
 2. Conversion would not have been necessary and script could point directly to meta.fasta
    but local files might not be on paths automatically mounted in the container.
*/
referencesDNA = Channel.from(params.references).map { (it.fasta).matches("^(https?|ftp)://.*\$") ? [it, file(workDir+'/REMOTE')] : [it, file(it.fasta)] }

process fetchReferenceForDNAAlignment {
  tag{meta.subMap(['species','version'])}
  //as above, storeDir not mounted accessible storeDir { (fasta.name).matches("REMOTE") ? (executor == 'awsbatch' ? "${params.outdir}/downloaded" : "downloaded") : null }

  input:
    set val(meta), file(fasta) from referencesDNA

  output:
    set val(meta), file("${basename}.fasta") into references4rnfSimReads, referencesForAlignersDNA

  when:
    'simulatedDNA'.matches(params.mode)

  script:
    //Abbreviate Genus_species name to G_species
    meta.species = (meta.species =~ /^./)[0]+(meta.species =~ /_.*$/)[0]
    meta.mode = 'DNA'
    basename=getTagFromMeta(meta)
    if((fasta.name).matches("REMOTE")) { //REMOTE FILE
      decompress = (meta.fasta).matches("^.*\\.gz\$") ?  "| gunzip --stdout " :  " "
      """
      curl ${meta.fasta} ${decompress} > ${basename}.fasta
      """
    } else if((fasta.name).matches("^.*\\.gz\$")){ //LOCAL GZIPPED
      """
      gunzip --stdout  ${fasta}  > ${basename}.fasta
      """
    } else { //LOCAL FLAT
      """
      cp -s  ${fasta} ${basename}.fasta
      """
    }
}


//DNA and RNA aligners in one channel as single indexing process defined
// alignersDNA.println { "$it DNA" }
// alignersRNA.println { "$it RNA" }
alignersDNA.join(alignersRNA , remainder: true)//.println { it }
  .map { [tool: it[0], dna: it[1]!=null, rna: it[2]!=null] }
  .set { aligners }

// aligners.combine(referencesForAlignersDNA).println { it }

// // // referencesForAlignersDNA.println { it }
// // // aligners.println { it }
process indexGenerator {
  label 'index'
  //label "${tool}" // it is currently not possible to set dynamic process labels in NF, see https://github.com/nextflow-io/nextflow/issues/894
  container { this.config.process.get("withLabel:${alignermeta.tool}" as String).get("container") }
  tag("${alignermeta.tool} << ${ref}")

  input:
    set val(alignermeta), val(refmeta), file(ref) from aligners.combine(refsRNA.mix(referencesForAlignersDNA))

  output:
    set val(meta), file("*") into indices4simulatedRNA, indices4realRNA, indices4simulatedDNA

  when: //check if dataset intended for {D,R}NA alignment reference and tool available for that purpose
    (refmeta.mode == 'DNA' && alignermeta.dna) || (refmeta.mode == 'RNA' && alignermeta.rna)
  // exec: //dev
  // meta =  alignermeta+refmeta//[target: "${ref}"]
  // println(meta)
  script:
    meta = [tool: "${alignermeta.tool}", target: "${ref}", mode: refmeta.mode]
    template "index/${alignermeta.tool}_index.sh" //points to e.g. biokanga_index.sh under templates/
}

process prepareDatasetsRNA {
  tag("${dataset}")

  input:
    set val(dataset), file(dataDir) from extractedDatasets

  output:
    set val(meta), file(r1), file(r2), file(cig) into preparedDatasets, prepareDatasetsForAdapters

  script:
  meta = [dataset: dataset, adapters: false]
  if(params.debug) { //FOR QUICKER RUNS, ONLY TAKE READS AND GROUDG-THRUTH FOR A SINGLE CHROMOSOME
    """
    awk '\$2~/${params.debugChromosome}\$/' ${dataDir}/*.cig \
      | sort -k1,1V --parallel ${task.cpus} \
      | tee >(awk -vOFS="\\t" 'NR%2==1{n++};{gsub(/[0-9]+/,n,\$1);print}' > cig) \
      | cut -f1 > debug.ids \
    && paste \
       <(paste - - < ${dataDir}/*.forward.fa) \
       <(paste - - < ${dataDir}/*.reverse.fa) \
       | awk -vOFS='\\t' 'NR==FNR{a[">"\$1]}; NR!=FNR && \$1 in a {n++; gsub(/[0-9]+/,n,\$1); gsub(/[0-9]+/,n,\$3); print}' \
         debug.ids - \
       | tee >(cut -f1,2 | tr '\\t' '\\n' > r1) | cut -f3,4 | tr '\\t' '\\n' > r2
    """
  } else {
    """
    ln -sf "\$(readlink -f ${dataDir}/*.forward.fa)" r1
    ln -sf "\$(readlink -f ${dataDir}/*.reverse.fa)" r2
    ln -sf "\$(readlink -f ${dataDir}/*.cig)" cig
    """
  }
}


process addAdaptersRNA {
  tag("${meta.dataset}")

  input:
    set val(inmeta), file(r1), file(r2), file(cig) from prepareDatasetsForAdapters

  output:
    set val(meta), file(a1), file(a2), file(cig)  into datasetsWithAdapters

  when:
    !params.debug || params.adapters //omitting this process to speed things up a bit for debug runs

  script:
    meta = inmeta.clone()
    meta.adapters = true
    """
    add_adapter2fasta_V3.pl ${r1} ${r2} a1 a2
    """
}

process alignSimulatedReadsRNA {
  label 'align'
  // label("${idxmeta.tool}") // it is currently not possible to set dynamic process labels in NF, see https://github.com/nextflow-io/nextflow/issues/894
  container { this.config.process.get("withLabel:${idxmeta.tool}" as String).get("container") }
  tag("${idxmeta} << ${readsmeta}")
  //GRAB CPU MODEL
  //afterScript 'hostname > .command.cpu; fgrep -m1 "model name" /proc/cpuinfo | sed "s/.*: //"  >> .command.cpu'


  input:
    set val(idxmeta), file("*"), val(readsmeta), file(r1), file(r2), file(cig) from indices4simulatedRNA.combine(datasetsWithAdapters.mix(preparedDatasets))

  output:
    set val(meta), file("*sam"), file(cig), file('.command.trace') into alignedDatasets

  script:
    meta = idxmeta.clone() + readsmeta.clone()
    template "rna/${idxmeta.tool}_align.sh"  //points to e.g. biokanga_align.sh in templates/
}


process nameSortSAM {
  label 'sort'
  label 'samtools'
  tag("${meta}")
  input:
    set val(meta), file(sam), file(cig) from alignedDatasets.map { meta, sam, cig, trace ->
        // meta.'aligntrace' = trace.splitCsv( header: true, limit: 1, sep: ' ')
        // meta.'aligntrace'.'duration' = trace.text.tokenize('\n').last()
        //meta.'aligntime' = trace.text.tokenize('\n').last()
        trace.splitEachLine("=", { record ->
        if(record.size() > 1 && record[0]=='realtime') { //to grab all, remove second condition and { meta."${record[0]}" = record[1] }
          meta.'aligntime'  = record[1]
        }
        new Tuple(meta, sam, cig)
      }

  output:
    set val(meta), file(sortedsam), file(cig) into sortedSAMs

  script:
    """
    samtools sort --threads ${task.cpus} -n --output-fmt BAM  ${sam} > sortedsam
    """
}


//Repeat downstream processes by either  leaving SAM as is or removing secondary & supplementary alignments
uniqSAM = Channel.from([false, true])

process fixSAM {
  label 'benchmark'
  tag("${meta}")

  input:
    set val(inmeta), file(sortedsam), file(cig), val(uniqed) from sortedSAMs.combine(uniqSAM)

  output:
    set val(meta), file(fixedsam), file(cig) into fixedSAMs

  when:
    uniqed == false || (params.uniqed == true && uniqed == true) //FILTERING SECONDARY&SUPPLEMENTARY IS OPTIONAL - GENERATES ADDITIONAL PLOTS

  script:
  meta = inmeta.clone() + [uniqed: uniqed]
  INSAM = uniqed ? "<(samtools view -F 2304 ${sortedsam})" : "<(samtools view ${sortedsam})"
  if(params.debug) {
    //1. should probably get exect value not just for --debug run. Otherwise --nummer fixed to 10 mil (?!)
    """
    fix_sam.rb --nummer \$(paste - - < ${cig} | wc -l) ${INSAM} | gzip -1c > fixedsam
    """
  } else {
    """
    fix_sam.rb ${INSAM} | gzip -1c > fixedsam
    """
  }
}

actions = Channel.from(['unique', 'multi'])
process compareToTruth {
  label 'benchmark'
  // label 'stats'
  tag("${outmeta}")

  input:
    set val(meta), file(fixedsam), file(cig), val(action) from fixedSAMs.combine(actions)
    // each action from actions

  output:
    set val(outmeta), file(stat) into stats

  script:
    outmeta = meta.clone() + [type : action]
    if(action == 'multi') {
      """
      gzip -dkc ${fixedsam} > sam
      compare2truth_multi_mappers.rb ${cig} sam > stat
      rm sam
      """
    } else {
      """
      gzip -dkc ${fixedsam} > sam
      compare2truth.rb ${cig} sam > stat
      rm sam
      """
    }
}

process tidyStats {
  label 'rscript'
  tag("${inmeta}")

  input:
    set val(inmeta), file(instats) from stats

  output:
    file 'tidy.csv' into tidyStats

  exec:
    meta = inmeta.clone()
    meta.replicate = meta.dataset[-1] //replicate num is last char
    meta.dataset = meta.dataset[0..-3] //strip of last 2 chars, eg. r1
    keyValue = meta.toMapString().replaceAll("[\\[\\],]","").replaceAll(':true',':TRUE').replaceAll(':false',':FALSE')


  shell:
    '''
    < !{instats} stats_parser.R !{meta.type} > tidy.csv
    for KV in !{keyValue}; do
      sed -i -e "1s/$/,${KV%:*}/" -e "2,\$ s/$/,${KV#*:}/" tidy.csv
    done
    '''
  //sed adds key to the header line and the value to each remaining line
}


process ggplot {
  tag 'figures'
  errorStrategy 'finish'
  label 'rscript'
  label 'figures'

  input:
    file csv from tidyStats.collectFile(name: 'all.csv', keepHeader: true)

  output:
    set file(csv), file('*.pdf') into plots

  shell:
    '''
    < !{csv} stats_figures.R
    '''
}



process downloadSRA {
  storeDir {executor == 'awsbatch' ? "${params.outdir}/downloaded" : "downloaded"} // storeDir "${workflow.workDir}/downloaded" put the datasets there and prevent generating cost to dataset creators through repeated downloads on re-runs
  scratch false

  output:
    file('*.sra') into downloadedSRA

  when:
    'realRNA'.matches(params.mode)

  script:
    """
    wget ftp://ftp.ddbj.nig.ac.jp/ddbj_database/dra/sralite/ByExp/litesra/SRX/SRX215/SRX2155547/SRR4228250/SRR4228250.sra
    """
}

process fromSRAtoFASTA {
  label 'sra'
  label 'slow'
  tag("${SRA}")

  input:
    file SRA from downloadedSRA

  output:
    // set file('*_1.fasta.gz'), file('*_2.fasta.gz') into sraFASTA
    set val(readsmeta), file('*_1.fasta'), file('*_2.fasta') into sraFASTA
    // val(readsmeta) into sraFASTA

  script:
  readsmeta = [sra: SRA.name[0..-5]]
  MAX_READS = params.debug ? '--maxSpotId 10000' : ''
  """
  fastq-dump --fasta 0 --split-files --origfmt --readids ${MAX_READS} ${SRA}
  """
}

process alignRealReadsRNA {
  label 'align'
  container { this.config.process.get("withLabel:${idxmeta.tool}" as String).get("container") }
  tag("${idxmeta} << ${readsmeta}")

  input:
    // set file(r1), file(r2) from sraFASTA
    // val(readsmeta) from sraFASTA
    // set val(readsmeta), file(r1), file(r2) from sraFASTA
    set val(idxmeta), file("*"), val(readsmeta), file(r1), file(r2) from indices4realRNA.combine(sraFASTA)

  output:
    set val(meta), file('*sam'), file('.command.trace') into alignedRealRNA

  script:
    meta = idxmeta.clone() + readsmeta.clone()
    template "rna/${idxmeta.tool}_align.sh"  //points to e.g. biokanga_align.sh in templates/
}

process samStatsRealRNA {
  echo true
  executor 'local'
  label 'samtools'
  tag("${inmeta}")

  input:
    set val(inmeta), file(sam) from alignedRealRNA.map { inmeta, sam, trace ->
        inmeta.'aligntime' = trace.text.tokenize('\n').last()
        new Tuple(inmeta, sam)
      }

  output:
    file 'csv' into statsRealRNA

  exec:
  // alntime = meta.tool+" "+(meta.aligntime.toLong()*10**-3/60)+" minutes"
  // """
  // echo ${alntime}
  // samtools flagstat sam \
  //   | sed -n '5p;9p' \
  //   | sed 's/^/${meta.tool}/g'
  // """
  meta = inmeta.clone()
  keyValue = meta.toMapString().replaceAll("[\\[\\],]","").replaceAll(':true',':TRUE').replaceAll(':false',':FALSE')
  shell:
    '''
    echo "aligned,paired" > csv
    samtools view -hF 2304 !{sam} | samtools flagstat - \
      | sed -n '5p;9p'  | cut -f1 -d' ' | paste - - | tr '\t' ',' >> csv
    for KV in !{keyValue}; do
      sed -i -e "1s/$/,${KV%:*}/" -e "2,\$ s/$/,${KV#*:}/" csv
    done
    '''
}

process ggplotRealRNA {
  tag 'figures'
  errorStrategy 'finish'
  label 'rscript'
  label 'figures'

  input:
    file csv from statsRealRNA.collectFile(name: 'real_RNA.csv', keepHeader: true)

  output:
    set file(csv), file('*.pdf') into plotsRealRNA

  shell:
  '''
  #!/usr/bin/env r

  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)

  stats <- read_csv("!{csv}")
  head(stats)

  ggplot(stats) +
    aes(aligned, aligntime*10^-3/60,colour=tool) +
    geom_point() +
    geom_label_repel(aes(label=tool)) +
    labs(title = "Accuracy vs alignment run times ",
        subtitle = "using 10 logical cores",
        x = "Aligned",
        y = "Run time (minutes)") +
    guides(label=FALSE, color=FALSE) #+
    #facet_wrap(adapters~dataset)
  ggsave(file="realRNA_aligned-runtime.pdf", width=16, height=9);

  ggplot(stats) +
    aes(paired, aligntime*10^-3/60,colour=tool) +
    geom_point() +
    geom_label_repel(aes(label=tool)) +
    labs(title = "Accuracy vs alignment run times ",
        subtitle = "using 10 logical cores",
        x = "Aligned as pairs",
        y = "Run time (minutes)") +
    guides(label=FALSE, color=FALSE)
  ggsave(file="realRNA_aligned-paired-runtime.pdf", width=16, height=9)
  '''
}



// ----- =======                   ======= -----
//                 DNA alignment
// ----- =======                   ======= -----





process indexReferences4rnfSimReadsDNA {
  tag{meta}
  label 'samtools'

  input:
    set val(meta), file(ref) from references4rnfSimReads

  output:
    set val(meta), file(ref), file('*.fai') into referencesWithIndex4rnfSimReads

  when:
    'simulatedDNA'.matches(params.mode) //only needed referencesLocal is a separate channel,

  script:
  """
  samtools faidx ${ref}
  """
}

process rnfSimReadsDNA {
  tag{simmeta}
  label 'rnftools'

  input:
    set val(meta), file(ref), file(fai) from referencesWithIndex4rnfSimReads
    each nsimreads from params.simreads.nreads.toString().tokenize(",")*.toInteger()
    each length from params.simreads.length.toString().tokenize(",")*.toInteger()
    each simulator from params.simreads.simulator
    each mode from params.simreads.mode //PE, SE
    each distance from params.simreads.distance //PE only
    each distanceDev from params.simreads.distanceDev //PE only

  output:
    set val(simmeta), file("*.fq.gz") into readsForAlignersDNA

  when:
    !(mode == "PE" && simulator == "CuReSim")

  script:
    tag=meta.species+"_"+meta.version+"_"+simulator
    simmeta = meta.subMap(['species','version'])+["simulator": simulator, "nreads":nsimreads, "mode": mode, "length": length ]
    len1 = length
    if(mode == "PE") {
      //FOR rnftools
      len2 = length
      tuple = 2
      dist="distance="+distance+","
      distDev= "distance_deviation="+distanceDev+","
      //FOR meta
      simmeta.dist = distance
      simmeta.distanceDev = distanceDev
    } else {
      len2 = 0
      tuple = 1
      dist=""
      distDev=""
    }
    """
    echo "import rnftools
    rnftools.mishmash.sample(\\"${tag}_reads\\",reads_in_tuple=${tuple})
    rnftools.mishmash.${simulator}(
            fasta=\\"${ref}\\",
            number_of_read_tuples=${nsimreads},
            ${dist}
            ${distDev}
            read_length_1=${len1},
            read_length_2=${len2}
    )
    include: rnftools.include()
    rule: input: rnftools.input()
    " > Snakefile
    snakemake \
    && for f in *.fq; do \
      paste - - - - < \${f} \
      | awk 'BEGIN{FS=OFS="\\t"};{gsub("[^ACGTUacgtu]","N",\$2); print}' \
      | tr '\\t' '\\n' \
      | gzip --stdout  --fast \
      > \${f}.gz \
      && rm \${f};
    done \
    && find . -type d -mindepth 2 | xargs rm -r
    """
}

process alignSimulatedReadsDNA {
  echo true
  label 'align'
  container { this.config.process.get("withLabel:${idxmeta.tool}" as String).get("container") } // label("${idxmeta.tool}") // it is currently not possible to set dynamic process labels in NF, see https://github.com/nextflow-io/nextflow/issues/894
  tag("${idxmeta} << ${readsmeta}")


input:
    set val(simmeta), file("?.fq.gz"), val(idxmeta), file('*') from readsForAlignersDNA.combine(indices4simulatedDNA) //cartesian product i.e. all input sets of reads vs all dbs

  // output:
  //   set val(alignmeta), file('out.bam') into bwaBAMs

  // when: //only align reads to the corresponding genome
  //   simmeta.species == dbmeta.species && simmeta.version == dbmeta.version

  """
  ls -la
  """

  // script:
  //   dbBasename=getTagFromMeta(dbmeta)
  //   alignmeta = dbmeta + simmeta
  //   alignmeta.aligner = "bwa"
  //   if(simmeta.mode == 'SE') {
  //     """
  //     bwa mem -t ${task.cpus} ${dbBasename} 1.fq.gz | samtools view -bS > out.bam
  //     """
  //   } else {
  //     """
  //     bwa mem -t ${task.cpus} ${dbBasename} 1.fq.gz 2.fq.gz | samtools view -bS > out.bam
  //     """
  //   }


  // script:
  //   meta = idxmeta.clone() + readsmeta.clone()
  //   template "dna/${idxmeta.tool}_align.sh"  //points to e.g. biokanga_align.sh in templates/
}


//WRAP-UP
writing = Channel.fromPath("${baseDir}/writing/*")
process render {
  tag 'manuscript'
  label 'rrender'
  label 'paper'
  stageInMode 'copy'
  scratch = true //hack, otherwise -profile singularity (with automounts) fails with FATAL:   container creation failed: unabled to {task.workDir} to mount list: destination ${task.workDir} is already in the mount point list

  input:
    file('*') from plots.flatten().toList()
    file('*') from plotsRealRNA.flatten().toList()
    file('*') from writing.collect()

  output:
    file '*'

  script:
  """
  #!/usr/bin/env Rscript

  library(rmarkdown)
  library(rticles)
  library(bookdown)

  rmarkdown::render(Sys.glob("*.Rmd"))
  """

}



// workflow.onComplete {
//     // any workflow property can be used here
//     println "Pipeline complete"
//     println "Command line: $workflow.commandLine"
//     println(workflow)
// }

// workflow.onError {
//     println "Oops .. something when wrong"
// }


