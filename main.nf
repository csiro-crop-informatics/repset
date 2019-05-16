#!/usr/bin/env nextflow

//For pretty-printing nested maps etc
import static groovy.json.JsonOutput.*

//RETURNS DNA ALIGNER NAMES/LABELS IF BOTH INDEXING AND ALIGNMENT TEMPLATES PRESENT
Channel.fromFilePairs("${workflow.projectDir}/templates/{index,dna}/*_{index,align}.sh", maxDepth: 1, checkIfExists: true)
  .filter{ params.alignersDNA == 'all' || it[0].matches(params.alignersDNA) }
  .map {
    params.defaults.alignersParams.DNA.putIfAbsent(it[0], [default: ''])  //make sure empty default param set available for every templated aligner
    params.defaults.alignersParams.DNA.(it[0]).putIfAbsent('default', '') //make sure empty default param set available for every templated aligner
    [it[0], "DNA"]
  }
  .set {alignersDNA}

//RETURNS RNA ALIGNER NAMES/LABELS IF BOTH INDEXING AND ALIGNMENT TEMPLATES PRESENT
Channel.fromFilePairs("${workflow.projectDir}/templates/{index,rna}/*_{index,align}.sh", maxDepth: 1, checkIfExists: true)
  .filter{ params.alignersRNA == 'all' || it[0].matches(params.alignersRNA) }
  .map {
    params.defaults.alignersParams.RNA.putIfAbsent(it[0], [default: ''])  //make sure empty default param set available for every templated aligner
    params.defaults.alignersParams.RNA.(it[0]).putIfAbsent('default', '') //make sure empty default param set available for every templated aligner
    [it[0], "RNA"]
  }
  .set { alignersRNA }

//DNA and RNA aligners in one channel as single indexing process defined
alignersDNA.join(alignersRNA , remainder: true)
  .map { [tool: it[0], dna: it[1]!=null, rna: it[2]!=null] }
  .set { aligners }

/*
 * Add to or overwrite map content recursively
 * Used to enable the use of NF -params-file opt such that params can be added and not just overwritten
 */
Map.metaClass.addNested = { Map rhs ->
    def lhs = delegate
    rhs.each { k, v -> lhs[k] = lhs[k] in Map ? lhs[k].addNested(v) : v }
    lhs
}

//Combine default and user parmas maps, then transform into a list and read into a channel to be consumed by alignment process(es)
alignersParamsList = []
params.defaults.alignersParams.addNested(params.alignersParams).each { seqtype, rnaOrDnaParams ->
  rnaOrDnaParams.each { tool, paramsets ->
    paramsets.each { paramslabel, ALIGN_PARAMS ->
      alignersParamsList << [tool: tool, paramslabel: paramslabel, seqtype: seqtype, ALIGN_PARAMS:ALIGN_PARAMS]
    }
  }
}
Channel.from(alignersParamsList).into {alignersParams4realDNA; alignersParams4SimulatedDNA}


/*
  Generic method for extracting a string tag or a file basename from a metadata map
 */
 def getTagFromMeta(meta, delim = '_') {
  return meta.species+delim+meta.version //+(trialLines == null ? "" : delim+trialLines+delim+"trialLines")
}

/*
  Simplistic method for checking if String is URL
 */
String.metaClass.isURL() {
   delegate.matches("^(https?|ftp)://.*\$")
}

def helpMessage() {
  log.info"""
  Usage:

  nextflow run csiro-crop-informatics/biokanga-manuscript -profile singularity
  nextflow run csiro-crop-informatics/biokanga-manuscript -profile docker

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



/*
 1. Input pointers to FASTA converted to files, NF would fetch remote as well and create tmp files,
    but avoiding that as may not scale with large genomes, prefer to do in process.
 2. Conversion would not have been necessary and script could point directly to meta.fasta
    but local files might not be on paths automatically mounted in the container.
*/
// referencesDNA = Channel.from(params.referencesDNA).map {
//   (it.fasta).matches("^(https?|ftp)://.*\$") ? [it, file(workDir+'/REMOTE')] : [it, file(it.fasta)]
// }
Channel.from(params.referencesDNA).map {
  new Tuple(
    it,
    it.fasta.isURL() ? file(workDir+'/REMOTE') : file(it.fasta),
    (it.containsKey('gff') ? (it.gff.isURL() ? file(workDir+'/REMOTE') : file(it.gff)) : file(workDir+'/NULL1')),
    (it.containsKey('bed') ? (it.bed.isURL() ? file(workDir+'/REMOTE') : file(it.bed)) : file(workDir+'/NULL2'))
  )
}.set { referencesDNA }

process fetchReference {
  echo true
  tag{meta.subMap(['species','version'])}
  //as above, storeDir not mounted accessible storeDir { (fasta.name).matches("REMOTE") ? (executor == 'awsbatch' ? "${params.outdir}/downloaded" : "downloaded") : null }

  input:
    set val(meta), file(fasta), file(gff), file(bed) from referencesDNA

  output:
    set val(meta), file("${basename}.fasta"), file("${basename}.bed"), file("${basename}.gff") into references4rnfSimReads, referencesForAlignersDNA, referencesForTranscriptomeExtraction

  when:
    'simulatedDNA'.matches(params.mode)

  script:
  // exec:
    //Abbreviate Genus_species name to G_species
    meta.species = (meta.species =~ /^./)[0]+(meta.species =~ /_.*$/)[0]
    meta.seqtype = 'DNA' // NEED TO HANDLE DIFFERENTLY
    basename=getTagFromMeta(meta)
    CMD = 'echo'
    [fasta: fasta, gff: gff, bed:bed].each { k, v ->
      meta."${k}" = "${basename}.${k}"
      if((v.name).matches("^NULL[0-9]+\$")) { //OPTIONAl FILE NOT SPECIFIED
        meta."${k}" = null
        CMD += " && touch ${basename}.${k}" //DUMMY FILE
      } else if((v.name).matches("REMOTE")) { //REMOTE FILE
        decompress = (meta.(k)).matches("^.*\\.gz\$") ?  "| gunzip --stdout " :  " "
        CMD += " && curl ${meta.(${k})} ${decompress} > ${basename}.${k}"
      } else if((v.name).matches("^.*\\.gz\$")){ //LOCAL GZIPPED
        CMD += " && gunzip --stdout  ${v}  > ${basename}.${k} "
      } else { //LOCAL FLAT
        CMD += " && cp -s  ${v} ${basename}.${k}"
      }
  }
  // println "$CMD"
  // println meta
  """
  ${CMD}
  """
}



// System.exit 0


process extarctTranscriptome {
  echo true
  scratch false
  container null
  module 'bedtools/2.28.0:samtools/1.9.0'
  input:
    set val(meta), file(ref), file(bed), file(gff) from referencesForTranscriptomeExtraction

  // output:
  //   set val(meta), file('trans.fa'), file('trans.fa.fai') into referencesWithIndex4rnfSimReads2

  script:
  println(prettyPrint(toJson(meta)))
  """
  ls -la
  """
  // """
  // bedtools getfasta \
  //   -fi ${ref} \
  //   -bed <(awk '\$8 ~ /mRNA/ ' ${bed}) \
  //   -name \
  //   | head -20 \
  //   > trans.fa \
  // && samtools faidx trans.fa
  // """
}

// // // // referencesForAlignersDNA.println { it }
// // // // aligners.println { it }
// process indexGenerator {
//   label 'index'
//   //label "${tool}" // it is currently not possible to set dynamic process labels in NF, see https://github.com/nextflow-io/nextflow/issues/894
//   container { this.config.process.get("withLabel:${alignermeta.tool}" as String).get("container") }
//   tag("${alignermeta.tool} << ${refmeta}")

//   input:
//     set val(alignermeta), val(refmeta), file(ref) from aligners.combine(referencesForAlignersDNA)

//   output:
//     set val(meta), file("*") into indices

//   when: //check if dataset intended for {D,R}NA alignment reference and tool available for that purpose
//     (refmeta.seqtype == 'DNA' && alignermeta.dna) || (refmeta.seqtype == 'RNA' && alignermeta.rna)
//   // exec: //dev
//   // meta =  alignermeta+refmeta//[target: "${ref}"]
//   // println(meta)
//   script:
//     meta = [tool: "${alignermeta.tool}", target: "${ref}"]+refmeta.subMap(['species','version','seqtype'])
//     template "index/${alignermeta.tool}_index.sh" //points to e.g. biokanga_index.sh under templates/
// }

// // process indexReferences4rnfSimReads {
// //   tag{meta}
// //   label 'samtools'

// //   input:
// //     set val(meta), file(ref) from references4rnfSimReads

// //   output:
// //     set val(meta), file(ref), file('*.fai') into referencesWithIndex4rnfSimReads

// //   when:
// //     'simulatedDNA'.matches(params.mode) //only needed referencesLocal is a separate channel,

// //   script:
// //   """
// //   samtools faidx ${ref}
// //   """
// // }

// // process optionalExtarctTranscriptome {
// //   scratch false
// //   container null
// //   module 'bedtools/2.28.0:samtools/1.9.0'
// //   input:
// //     set val(meta), file(ref), file(fai), file(bed) from referencesWithIndex4rnfSimReads.combine(Channel.fromPath(params.referencesDNA[0].bed))

// //   output:
// //     set val(meta), file('trans.fa'), file('trans.fa.fai') into referencesWithIndex4rnfSimReads2

// //   script:
// //   """
// //   bedtools getfasta \
// //     -fi ${ref} \
// //     -bed <(awk '\$8 ~ /mRNA/ ' ${bed}) \
// //     -name \
// //     | head -20 \
// //     > trans.fa \
// //   && samtools faidx trans.fa
// //   """
// // }

// process rnfSimReadsDNA {
//   tag{simmeta}
//   label 'rnftools'

//   input:
//     // set val(meta), file(ref), file(fai) from referencesWithIndex4rnfSimReads
//     set val(meta), file(ref) from references4rnfSimReads
//     each nsimreads from params.simreadsDNA.nreads.toString().tokenize(",")*.toInteger()
//     each length from params.simreadsDNA.length.toString().tokenize(",")*.toInteger()
//     each simulator from params.simreadsDNA.simulator
//     each mode from params.simreadsDNA.mode //PE, SE
//     each distance from params.simreadsDNA.distance //PE only
//     each distanceDev from params.simreadsDNA.distanceDev //PE only

//   output:
//     set val(simmeta), file("*.fq.gz") into readsForAlignersDNA

//   when:
//     !(mode == "PE" && simulator == "CuReSim")

//   script:
//     tag=meta.species+"_"+meta.version+"_"+simulator
//     simmeta = meta.subMap(['species','version'])+["simulator": simulator, "nreads":nsimreads, "mode": mode, "length": length]
//     len1 = length
//     if(mode == "PE") {
//       //FOR rnftools
//       len2 = length
//       tuple = 2
//       dist="distance="+distance+","
//       distDev= "distance_deviation="+distanceDev+","
//       //FOR meta
//       simmeta.dist = distance
//       simmeta.distanceDev = distanceDev
//     } else {
//       len2 = 0
//       tuple = 1
//       dist=""
//       distDev=""
//     }
//     """
//     echo "import rnftools
//     rnftools.mishmash.sample(\\"${tag}_reads\\",reads_in_tuple=${tuple})
//     rnftools.mishmash.${simulator}(
//             fasta=\\"${ref}\\",
//             number_of_read_tuples=${nsimreads},
//             ${dist}
//             ${distDev}
//             read_length_1=${len1},
//             read_length_2=${len2}
//     )
//     include: rnftools.include()
//     rule: input: rnftools.input()
//     " > Snakefile
//     snakemake -p \
//     && for f in *.fq; do \
//       paste - - - - < \${f} \
//       | awk 'BEGIN{FS=OFS="\\t"};{gsub("[^ACGTUacgtu]","N",\$2); print}' \
//       | tr '\\t' '\\n' \
//       | gzip --stdout  --fast \
//       > \${f}.gz \
//       && rm \${f};
//     done \
//     && find . -type d -mindepth 2 | xargs rm -r
//     """
// }

// process alignSimulatedReads {
//   label 'align'
//   container { this.config.process.get("withLabel:${idxmeta.tool}" as String).get("container") } // label("${idxmeta.tool}") // it is currently not possible to set dynamic process labels in NF, see https://github.com/nextflow-io/nextflow/issues/894
//   tag("${idxmeta.subMap(['tool','species'])} << ${simmeta.subMap(['simulator','nreads'])} @ ${paramsmeta.subMap(['paramslabel'])}")

//   input:
//     set val(simmeta), file("?.fq.gz"), val(idxmeta), file('*'), val(paramsmeta) from readsForAlignersDNA.combine(indices).combine(alignersParams4SimulatedDNA) //cartesian product i.e. all input sets of reads vs all dbs

//   output:
//     set val(alignmeta), file('out.?am') into alignedSimulatedDNA

//   when: //only align DNA reads to the corresponding genome, using the corresponding params set
//     idxmeta.seqtype == 'DNA' && idxmeta.species == simmeta.species && idxmeta.version == simmeta.version && paramsmeta.tool == idxmeta.tool && paramsmeta.seqtype == 'DNA'

//   script:
//     alignmeta = idxmeta.clone() + simmeta.clone() + paramsmeta.clone()
//     // if(simmeta.mode == 'PE') {
//     ALIGN_PARAMS = paramsmeta.ALIGN_PARAMS
//       template "dna/${idxmeta.tool}_align.sh"  //points to e.g. biokanga_align.sh in templates/
//     // } else {
//     //   template "dna/${idxmeta.tool}_align.sh"  //points to e.g. biokanga_align.sh in templates/
//     // }
// }

// process rnfEvaluateSimulated {
//   label 'rnftools'
//   tag{alignmeta.subMap(['tool','simulator','target','paramslabel'])}


//   input:
//     set val(alignmeta), file(samOrBam) from alignedSimulatedDNA

//   output:
//      set val(alignmeta), file(summary) into summariesSimulatedDNA
//      set val(alignmeta), file(detail) into detailsSimulatedDNA

//   script:
//   // println prettyPrint(toJson(alignmeta))
//   """
//   paste \
//     <( rnftools sam2es -i ${samOrBam} -o - | awk '\$1 !~ /^#/' \
//       | tee ES \
//       | tee >( awk -vOFS="\\t" '{category[\$7]++}; END{for(k in category) {print k,category[k]}}' > summary ) \
//     ) \
//     <( samtools view ${samOrBam} ) \
//   | awk -vOFS="\\t" '{if(\$1 == \$9 && \$5 == \$12){print \$11,\$12,\$7} else {print "BAM - ES mismatch, terminating",\$0 > "/dev/stderr"; exit 1}}' > detail
//   """

// // rnftools sam2es OUTPUT header
// // # RN:   read name
// // # Q:    is mapped with quality
// // # Chr:  chr id
// // # D:    direction
// // # L:    leftmost nucleotide
// // # R:    rightmost nucleotide
// // # Cat:  category of alignment assigned by LAVEnder
// // #         M_i    i-th segment is correctly mapped
// // #         m      segment should be unmapped but it is mapped
// // #         w      segment is mapped to a wrong location
// // #         U      segment is unmapped and should be unmapped
// // #         u      segment is unmapped and should be mapped
// // # Segs: number of segments
// // #
// // # RN    Q       Chr     D       L       R       Cat     Segs
// }

// process collateDetailsSimulatedDNA {
//   label 'stats'
//   executor 'local' //explicit to avoid a warning being prined. Either way must be local exec as no script block for this process just nextflow/groovy exec

//   input:
//     val collected from detailsSimulatedDNA.collect()

//   output:
//     file 'details.tsv' into collatedDetailsSimulatedDNA

//   exec:
//   def outfileTSV = task.workDir.resolve('details.tsv')
//   i = 0;
//   sep = "\t"
//   header = "Species\tChromosome\tPosition\tClass\tSimulator\tAligner\tMode\n"
//   // outfileTSV << header
//   outfileTSV.withWriter { target ->
//     target << header
//     collected.each {
//       if(i++ %2 == 0) {
//         meta = it
//       } else {
//         common = meta.simulator+sep+meta.tool+sep+meta.mode+"\n"
//         it.withReader { source ->
//           String line
//           while( line=source.readLine() ) {
//             StringBuilder sb = new StringBuilder()
//             sb.append(meta.species).append(sep).append(line).append(sep).append(common)
//             target << sb
//             // target << meta.species+sep+line+sep+common
//           }
//         }
//       }
//       // it.eachLine { line ->
//       //   outfileTSV << meta.species+sep+line+sep+common
//       // }
//     }
//   }
// }

// process collateSummariesSimulatedDNA {
//   label 'stats'
//   executor 'local' //explicit to avoid a warning being prined. Either way must be local exec as no script block for this process just nextflow/groovy exec

//   input:
//     val collected from summariesSimulatedDNA.collect()

//   output:
//     // set file('summaries.csv'), file('summaries.json') into collatedSummariesSimulatedDNA
//     set file('summaries.json'), file('categories.json') into collatedSummariesSimulatedDNA

//   exec:
//   def outfileJSON = task.workDir.resolve('summaries.json')
//   def categoriesJSON = task.workDir.resolve('categories.json')
//   // def outfileCSV = task.workDir.resolve('summaries.csv')
//   categories = ["M_1":"First segment is correctly mapped", "M_2":"Second segment is correctly mapped",
//   "m":"segment should be unmapped but it is mapped", "w":"segment is mapped to a wrong location",
//   "U":"segment is unmapped and should be unmapped", "u":"segment is unmapped and should be mapped"]
//   categoriesJSON << prettyPrint(toJson(categories))
//   entry = null
//   entries = []
//   // entries << [categories: categories]
//   i=0;
//   TreeSet headersMeta = []
//   TreeSet headersResults = []
//   collected.each {
//     if(i++ %2 == 0) {
//       if(entry != null) {
//         entries << entry
//         entry.meta.each {k,v ->
//           headersMeta << k
//         }
//       }
//       entry = [:]
//       entry.meta = it.clone()
//     } else {
//       entry.results = [:]
//       it.eachLine { line ->
//         (k, v) = line.split()
//         entry.results << [(k) : v ]
//         //entry.results << [(categories[(k)]) : v ]
//         headersResults << (k)
//         //headersResults << (categories[(k)])
//       }
//     }
//   }
//   entries << entry
//   outfileJSON << prettyPrint(toJson(entries))

//   // //GENERATE CSV OUTPUT
//   // SEP=","
//   // outfileCSV << headersMeta.join(SEP)+SEP+headersResults.join(SEP)+"\n"
//   // entries.each { entry ->
//   //   line = ""
//   //   headersMeta.each { k ->
//   //     val = "${entry.meta[k]}".isNumber() ? entry.meta[k] :  "\"${entry.meta[k]}\""
//   //     line += line == "" ? val : (SEP+val)
//   //   }
//   //   headersResults.each { k ->
//   //     value = entry.results[k]
//   //     line += SEP
//   //     // println(k + ' -> ' + value)
//   //     line += value == null ? 0 : (value.isNumber() ? value : "\"${value}\"") //NOT QUITE RIGHT, ok for 'w' not for 'u'
//   //   }
//   //   outfileCSV << line+"\n"
//   // }

// }

// process plotDetailSimulatedDNA {
//   label 'rscript'
//   label 'figures'

//   input:
//     file '*' from collatedDetailsSimulatedDNA

//   output:
//     file '*' into collatedDetailsPlotsSimulatedDNA

//   script:        //============================ TODO : move under bin/
//   binWidth='1E5'
//   """
//   touch plotPlaceholderD
//   #< !{csv} plot_details_simulatedDNA.R
//   """

// }

// process plotSummarySimulatedDNA {
//   label 'rscript'
//   label 'figures'

//   input:
//     // set file(csv), file(json) from collatedSummariesSimulatedDNA
//     set file(json), file(categories) from collatedSummariesSimulatedDNA

//   output:
//     file '*' into collatedSummariesPlotsSimulatedDNA

//   shell:
//   '''
//   < !{json} plot_simulatedDNA.R
//   '''
// }

// //WRAP-UP
// writing = Channel.fromPath("$baseDir/report/*").mix(Channel.fromPath("$baseDir/manuscript/*")) //manuscript dir exists only on manuscript branch

// process render {
//   tag {"Render ${Rmd}"}
//   label 'rrender'
//   label 'report'
//   stageInMode 'copy'
//   //scratch = true //hack, otherwise -profile singularity (with automounts) fails with FATAL:   container creation failed: unabled to {task.workDir} to mount list: destination ${task.workDir} is already in the mount point list

//   input:
//     // file('*') from plots.flatten().toList()
//     // file('*') from plotsRealRNA.flatten().toList()
//     file(Rmd) from writing
//     file('*') from collatedDetailsPlotsSimulatedDNA.collect()
//     file('*') from collatedSummariesPlotsSimulatedDNA.collect()

//   output:
//     file '*'

//   script:
//   """
//   #!/usr/bin/env Rscript

//   library(rmarkdown)
//   library(rticles)
//   library(bookdown)

//   rmarkdown::render("${Rmd}")
//   """
// }