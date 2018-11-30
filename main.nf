#!/usr/bin/env nextflow

//RETURNS ALIGNER NAMES/LABELS IF BOTH INDEXING AND ALIGNMENT TEMPLATES PRESENT
aligners = Channel.fromFilePairs("${workflow.projectDir}/templates/*_{index,align}.sh", maxDepth: 1, checkIfExists: true)
  .map { it[0] }
  .filter{ params.aligners == 'all' || it.matches(params.aligners) }
  // .filter{ !it.matches("subread\$") } //temp
  // .filter{ !params.debug || it.matches("(biokanga|bowtie2|bwa|dart|hisat2|star)\$") }

//Pre-computed BEERS datasets
datasets = Channel.from(['human_t1r1','human_t1r2','human_t1r3','human_t2r1','human_t2r2','human_t2r3','human_t3r1','human_t3r2','human_t3r3'])
  .filter{ !params.debug || it == params.debugDataset }
  .filter{ (it[-1] as Integer) <= params.replicates}

//Download reference: hg19
url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz'

//For pretty-printing nested maps etc
import static groovy.json.JsonOutput.*

def helpMessage() {
  log.info"""
  Usage:

  nextflow run csiro-crop-informatics/biokanga-manuscript -profile singularity

  Default params:
  """.stripIndent()
  // println(prettyPrint(toJson(params)))
  println(prettyPrint(toJson(config)))
  // println(prettyPrint(toJson(config.process)))
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

process downloadReference {
  storeDir {executor == 'awsbatch' ? "${params.outdir}/downloaded" : "downloaded"}
  scratch false

  input:
    val(url)

  output:
    file('chromFa.tar.gz') into refs

  script:
  """
  wget ${url}
  """
}

process downloadDatasets {
  // label 'download'
  tag("${dataset}")
  storeDir {executor == 'awsbatch' ? "${params.outdir}/downloaded" : "downloaded"} // storeDir "${workflow.workDir}/downloaded" put the datasets there and prevent generating cost to dataset creators through repeated downloads on re-runs
  scratch false

  input:
    val(dataset) from datasets

  output:
    file("${dataset}.tar.bz2") into downloadedDatasets
    // set val("${dataset}"), file("${dataset}.tar.bz2") into downloadedDatasets  //not possible when using storeDir

  script:
    """
    wget http://bp1.s3.amazonaws.com/${dataset}.tar.bz2
    """
}

process extractDatasets {
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

process convertReference {
  input:
    file(downloadedRef) from refs

  output:
    file ref
    // file('ucsc.hg19.fa') into reference

  script:
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

process indexGenerator {
  label 'index'
  //label "${tool}" // it is currently not possible to set dynamic process labels in NF, see https://github.com/nextflow-io/nextflow/issues/894
  container { this.config.process.get("withLabel:${tool}" as String).get("container") }
  tag("${tool} << ${ref}")

  input:
    file ref
    val(tool) from aligners

  output:
    set val(meta), file("*") into indices

  script:
    meta = [tool: "${tool}", target: "${ref}"]
    template "${tool}_index.sh" //points to e.g. biokanga_index.sh in templates/
}


process prepareDatasets {
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


process addAdapters {
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

// optimised settings
// hisat2 : default-1-20-0.5-25-5-12-1-0-3-0
//  hisat2 -N 1 -L 20 -i S,1,0.5 -D 25 -R 5 --pen-noncansplice 12  --mp 1,0  --sp 3,0
//  hisat2 --end-to-end -N <NUM_MISMATCH> -L <SEED_LENGTH> -i S,1,<SEED_INTERVAL> -D <SEED_EXTENSION> -R <RE_SEED>
     //--pen-noncansplice <PENALITY_NONCANONICAL> --mp <MAX_MISMATCH_PENALITY>,<MIN_MISMATCH_PENALITY>
     //--sp <MAX_SOFTCLIPPING_PENALITY>,<MIN_SOFTCLIPPING_PENALITY>--time --reorder --known-splicesite-infile <output index path>/<genome name>.splicesites.txt --novel-splicesite-outfile splicesites.novel.txt --novel-splicesite-infile splicesites.novel.txt -f -x <index name> -1 <read file 1> -2 <read file 2> -S <output sam file>
//  ENDTOEND_MODE - NUM_MISMATCH - SEED_LENGTH - SEED_INTERVAL - SEED_EXTENSION - RE_SEED - PENALITY_NONCANONICAL - MAX_MISMATCH_PENALITY - MIN_MISMATHC_PENALITY - MAX_SOFTCLIPPING_PENALITY - MIN_SOFTCLIPPING_PENALITY
// gsnap : 15-1-10-221-41
//  gsnap --max-mismatches <MAX_MISMATCHES> --indel-penalty <INDEL_PENALITY> --gmap-min-match-length <GMAP_MIN_MATCH_LENGTH> --pairexpect <PAIR_EXPECT> --pairdev <PAIR_DEV> --merge-distant-samechr --ordered --novelsplicing 1 --use-splicing <index name>.splicesites --nthreads 16 --batch 5 --expand-offsets 1 <read file 1> <read file 2> > <output sam file>
//  MAX_MISMATCHES - INDEL_PENALITY - GMAP_MIN_MATCH_LENGTH - PAIR_EXPECT - PAIR_DEV
// star: 1000000-1000000-100-33-0.3-12-15-Local-0-0.3-50-3-BySJout
//  STAR --twopassMode Basic --outSAMunmapped Within --limitOutSJcollapsed <NUM_COLLAPSED_JUNCTIONS> --limitSjdbInsertNsj <NUM_INSERTED_JUNCTIONS> --outFilterMultimapNmax <NUM_MULTIMAPPER> --outFilterMismatchNmax <NUM_FILTER_MISMATCHES> --outFilterMismatchNoverLmax <RATIO_FILTER_MISMATCHES> --seedSearchStartLmax <SEED_LENGTH> --alignSJoverhangMin <OVERHANG> --alignEndsType <END_ALIGNMENT_TYPE> --outFilterMatchNminOverLread <NUM_FILTER_MATCHES> --outFilterScoreMinOverLread <NUM_FILTER_SCORE> --winAnchorMultimapNmax <NUM_ANCHOR> --alignSJDBoverhangMin <OVERHANG_ANNOTATED> --outFilterType <OUT_FILTER>
//  STAR --twopassMode Basic --outSAMunmapped Within --limitOutSJcollapsed 1000000 --limitSjdbInsertNsj 1000000 --outFilterMultimapNmax <NUM_MULTIMAPPER> --outFilterMismatchNmax <NUM_FILTER_MISMATCHES> --outFilterMismatchNoverLmax <RATIO_FILTER_MISMATCHES> --seedSearchStartLmax <SEED_LENGTH> --alignSJoverhangMin <OVERHANG> --alignEndsType <END_ALIGNMENT_TYPE> --outFilterMatchNminOverLread <NUM_FILTER_MATCHES> --outFilterScoreMinOverLread <NUM_FILTER_SCORE> --winAnchorMultimapNmax <NUM_ANCHOR> --alignSJDBoverhangMin <OVERHANG_ANNOTATED> --outFilterType <OUT_FILTER>
//  NUM_COLLAPSED_JUNCTIONS - NUM_INSERTED_JUNCTIONS - NUM_MULTIMAPPER - NUM_FILTER_MISMATCHES - RATIO_FILTER_MISMATCHES - SEED_LENGTH - OVERHANG - END_ALIGNMENT_TYPE - NUM_FILTER_MATCHES - NUM_FILTER_SCORE - NUM_ANCHOR - OVERHANG_ANNOTATED - OUT_FILTER


process align {
  label 'align'
  // label("${idxmeta.tool}") // it is currently not possible to set dynamic process labels in NF, see https://github.com/nextflow-io/nextflow/issues/894
  container { this.config.process.get("withLabel:${idxmeta.tool}" as String).get("container") }
  tag("${idxmeta} << ${readsmeta}")
  // cache { idxmeta.tool == 'star' ? 'deep' : 'lenient'}
  // time { idxmeta.tool == 'gsnap' ? '4.h' : '2.h'}
  //GRAB CPU MODEL
  //afterScript 'hostname > .command.cpu; fgrep -m1 "model name" /proc/cpuinfo | sed "s/.*: //"  >> .command.cpu'

  input:
    set val(idxmeta), file("*"), val(readsmeta), file(r1), file(r2), file(cig) from indices.combine(datasetsWithAdapters.mix(preparedDatasets))

  output:
    set val(meta), file("*sam"), file(cig), file('.command.trace') into alignedDatasets

  script:
    meta = idxmeta.clone() + readsmeta.clone()
    template "${idxmeta.tool}_align.sh"  //points to e.g. biokanga_align.sh in templates/
}


process nameSortSAM {
  label 'sort'
  label 'samtools'
  tag("${meta}")
  input:
    set val(meta), file(sam), file(cig) from alignedDatasets.map { meta, sam, cig, trace ->
        // meta.'aligntrace' = trace.splitCsv( header: true, limit: 1, sep: ' ')
        // meta.'aligntrace'.'duration' = trace.text.tokenize('\n').last()
        meta.'aligntime' = trace.text.tokenize('\n').last()
        new Tuple(meta, sam, cig)
      }

  output:
    set val(meta), file(sortedsam), file(cig) into sortedSAMs

  script:
    """
    samtools sort --threads ${task.cpus} -n --output-fmt BAM  ${sam} > sortedsam
    """
    // """
    // samtools view -F2304  ${sam} | grep -v '^@' | sort -t'.' -k2,2n --parallel ${task.cpus} > sortedsam
    // """
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

writing = Channel.fromPath("${baseDir}/writing/*")
process render {
  tag 'manuscript'
  label 'rrender'
  label 'paper'
  stageInMode 'copy'

  input:
    file('*') from plots.flatten().toList()
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
