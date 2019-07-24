#!/usr/bin/env nextflow

//For pretty-printing nested maps etc
import static groovy.json.JsonGenerator.*
// import static groovy.json.JsonSlurper as JsonSlurper

//Otherwise JSON generation triggers stackoverflow when encountering Path objects
jsonGenerator = new groovy.json.JsonGenerator.Options()
                .addConverter(java.nio.file.Path) { java.nio.file.Path p, String key -> p.toUriString() }
                .build()


//Input validation specified elswhere
def validators = new GroovyShell().parse(new File("${baseDir}/groovy/Validators.groovy"))

//Read, parse, validate and sanitize alignment/mapping tools config
def allRequired = ['tool','version','container','index'] //Fields required for each tool in config
def allModes = 'dna2dna|rna2rna|rna2dna' //At leas one mode has to be defined as supported by each tool
def allVersions = validators.validateMappersDefinitions(params.mappersDefinitions, allRequired, allModes)

//Check if specified template files exist
validators.validateTemplatesAndScripts(params.mappersDefinitions, (['index']+(allModes.split('\\|') as List)), "${baseDir}/templates")

//Read, sanitize and validate alignment/mapping param sets
validators.validateMapperParamsDefinitions(params.mapperParamsDefinitions, allVersions)

//Validated now, so gobble up mappers and their params definitions
mappersChannel = Channel.from(params.mappersDefinitions)
  .filter{ params.mappers == 'all' || it.tool.matches(params.mappers) } //TODO Could allow :version

mappersParamsChannel = Channel.from(params.mapperParamsDefinitions)

//one or more mapping mode
mapModesChannel = Channel.from(params.mapmode.split('\\||,'))





// //RETURNS DNA2DNA ALIGNER NAMES/LABELS IF BOTH INDEXING AND ALIGNMENT TEMPLATES PRESENT
// Channel.fromFilePairs("${workflow.projectDir}/templates/{index,dna2dna}/*_{index,align}.sh", maxDepth: 1, checkIfExists: true)
//   .map {
//     params.defaults.alignersParams.DNA2DNA.putIfAbsent(it[0], [default: ''])  //make sure empty default param set available for every templated aligner
//     params.defaults.alignersParams.DNA2DNA.(it[0]).putIfAbsent('default', '') //make sure empty default param set available for every templated aligner
//     [it[0], "DNA2DNA"]
//   }
//   .set {alignersDNA2DNA}

// //RETURNS RNA2DNA ALIGNER NAMES/LABELS IF BOTH INDEXING AND ALIGNMENT TEMPLATES PRESENT
// Channel.fromFilePairs("${workflow.projectDir}/templates/{index,rna2dna}/*_{index,align}.sh", maxDepth: 1, checkIfExists: true)
//   .map {
//     params.defaults.alignersParams.RNA2DNA.putIfAbsent(it[0], [default: ''])  //make sure empty default param set available for every templated aligner
//     params.defaults.alignersParams.RNA2DNA.(it[0]).putIfAbsent('default', '') //make sure empty default param set available for every templated aligner
//     [it[0], "RNA2DNA"]
//   }
//   .set { alignersRNA2DNA }

// //RETURNS RNA2RNA ALIGNER NAMES/LABELS IF BOTH INDEXING AND ALIGNMENT TEMPLATES PRESENT
// Channel.fromFilePairs("${workflow.projectDir}/templates/{index,rna2rna}/*_{index,align}.sh", maxDepth: 1, checkIfExists: true)
//   .map {
//     params.defaults.alignersParams.RNA2RNA.putIfAbsent(it[0], [default: ''])  //make sure empty default param set available for every templated aligner
//     params.defaults.alignersParams.RNA2RNA.(it[0]).putIfAbsent('default', '') //make sure empty default param set available for every templated aligner
//     [it[0], "RNA2RNA"]
//   }
//   .set { alignersRNA2RNA }

// //DNA and RNA aligners in one channel as single indexing process defined
// alignersDNA2DNA.mix(alignersRNA2DNA).mix(alignersRNA2RNA)
//   .groupTuple(size:3, remainder: true, sort:true)
//   .map { tool, alnModes ->
//     [tool: tool, modes: [dna2dna: alnModes.contains('DNA2DNA'), rna2dna: alnModes.contains('RNA2DNA'), rna2rna: alnModes.contains('RNA2RNA')]]
//   }
//  .filter{ params.aligners == 'all' || it.tool.matches(params.aligners) }
//  .set { aligners }

/*
 * Add to or overwrite map content recursively
 * Used to enable the use of NF -params-file opt such that params can be added and not just overwritten
 */
Map.metaClass.addNested = { Map rhs ->
    def lhs = delegate
    rhs.each { k, v -> lhs[k] = lhs[k] in Map ? lhs[k].addNested(v) : v }
    lhs
}

// //Combine default and user parmas maps, then transform into a list and read into a channel to be consumed by alignment process(es)
// alignersParamsList = []
// params.defaults.alignersParams.addNested(params.alignersParams).each { alignMode, rnaOrDnaParams ->
//   rnaOrDnaParams.each { tool, paramsets ->
//     paramsets.each { paramslabel, ALIGN_PARAMS ->
//       alignersParamsList << [tool: tool, paramslabel: paramslabel, alignMode: alignMode, ALIGN_PARAMS:ALIGN_PARAMS]
//     }
//   }
// }
// Channel.from(alignersParamsList).set {alignersParams}
// // Channel.from(alignersParamsList).into {alignersParams4realDNA; alignersParams4SimulatedDNA}


// // println(groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(alignersParamsList)))



/*
  Generic method for extracting a string tag or a file basename from a metadata map
 */
def getTagFromMeta(meta, delim = '_') {
  return meta.species+delim+meta.version //+(trialLines == null ? "" : delim+trialLines+delim+"trialLines")
}

/*
  Given a file with '=' delimited key value pairs on each line
  (this could e.g. be .command.trace)
  parse and store in map provided,
 */
def parseFileMap(filemap, map, subset = false) {
  filemap.splitEachLine("=", { record ->
      if(record.size() > 1 && (subset == false || record[0] in subset)) {
        v = record[1]
        map."${record[0]}" = v.isInteger() ? v.toInteger() : v.isDouble() ? v.toDouble() : v
      }
    })
}

// /*
//   Delayed interpolation
// */
// def String bindTemplate(template, binding) {
//   def engine = new groovy.text.SimpleTemplateEngine()
//   engine.createTemplate(template).make(binding).toString()
// }

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
fastaChn = Channel.create()
gffChn = Channel.create()

Channel.from(params.references).separate (fastaChn, gffChn) { it ->
  onlyFasta = it.clone()
  onlyFasta.remove('gff')
  onlyGff = it.clone()
  onlyGff.remove('fasta')
  fasta = [onlyFasta, it.fasta.isURL() ? file(workDir+'/REMOTE1') : file(it.fasta)]
  it.containsKey('gff') ? [fasta, [onlyGff, it.gff.isURL() ? file(workDir+'/REMOTE2') : file(it.gff)]] : [fasta, [onlyGff, file(workDir+'/NULL')]]
}

process stageInputFile {
  tag{meta.subMap(['species','version'])+[fileType: fileType]}

  //as above, storeDir not mounted accessible
  // storeDir { (infile.name).matches("REMOTE") ? (executor == 'awsbatch' ? "${params.outdir}/downloaded" : "downloaded") : null }
  // storeDir { (infile.name).matches("REMOTE.*") ? "${params.outdir}/downloaded" : null }
  storeDir { (infile.name).matches("REMOTE.*") ? "${workDir}/downloaded" : null } //- perhaps more robust as workdir is mounted in singularity unlike outdir?

  input:
    set val(meta), file(infile) from fastaChn.mix(gffChn)

  output:
    set val(outmeta), file(outfile)  into stagedFiles
    // set val(outmeta), file(outfile), val(meta.seqtype)  into stagedFiles

  when:
    !infile.name.matches('NULL')

  script:
    basename=getTagFromMeta(meta)
    // println(prettyPrint(toJson(meta)))
    fileType = meta.containsKey('fasta') ? 'fasta' : 'gff'
    outfile =  "${basename}.${fileType}"
    outmeta = meta.subMap(['species', 'version','seqtype'])
    if(infile.name.matches("REMOTE.*")) { //REMOTE FILE
      remoteFileName = meta."${fileType}"
      decompress = remoteFileName.matches("^.*\\.gz\$") ?  "| gunzip --stdout " :  " "
      """curl ${remoteFileName} ${decompress} > ${outfile}"""
    } else if((infile.name).matches("^.*\\.gz\$")){ //LOCAL GZIPPED
      """gunzip --stdout  ${infile}  > ${outfile} """
    } else { //LOCAL FLAT
      """cp -s  ${infile} ${outfile}"""
    }
}

// referencesOnly = Channel.create()
// referencesForTranscriptomeExtraction = Channel.create()
stagedFiles
  // .view()
  .groupTuple() //match back fasta with it's gffif available
  // .view { meta, files, seqtype -> "meta: ${meta}\nfiles: ${files}\nseqtype: ${seqtype}"}
  // .map { meta, files ->
  //   files.sort { a,b -> a.name.substring(a.name.lastIndexOf('.')+1) <=> b.name.substring(b.name.lastIndexOf('.')+1) } //ensure gff goes after fasta based on extension
  //   [meta, files]
  // }
  // .view { it -> println(groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it)))}
  .set { stagedReferences }


process faidxGenomeFASTA {
  tag("${refmeta}")
  label 'samtools'

  input:
    set val(refmeta), file(infiles) from stagedReferences

  output:
    set val(refmeta), file("${ref}.fai") into genomeIndicesForReadCoordinateConversion
    set val(refmeta), file(ref), file("${ref}.fai") into genomesForIndexing, genomesForRnfSimReads
    set val(refmeta), file(ref), file("${ref}.fai"), file("${gff}") optional true into referencesForTranscriptomeExtraction //refsForIndexing

  script:
  ref = infiles[infiles.findIndexOf { fname -> fname =~ /\.fasta$/ }]
  gffIdx = infiles.findIndexOf { fname -> fname =~ /\.gff$/ }
  gff = gffIdx >= 0 ? infiles[gffIdx] : 'NO_GFF_HERE'
  """
  samtools faidx ${ref}
  """
}

// referencesOnly.view {println(groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it)))}

// referencesOnly
//   // .map { [it[0], it[1][0]]} //meta, fasta_file
//   // .view()
//   .map { meta, files ->
//     // meta.seqtype = meta.seqtype[0] //un-list
//     // files.sort { a,b -> a.name.substring(a.name.lastIndexOf('.')+1) <=> b.name.substring(b.name.lastIndexOf('.')+1) } //ensure gff goes after fasta based on extension
//     [meta, files[0]] //meta, fasta_file
//   }
//   // .view()
//   .into {referencesForAligners; references4rnfSimReads}


process extarctTranscripts {
  echo true
  label 'gffread'
  label 'slow'
  tag{meta.subMap(['species','version'])}
  scratch false

  //SLOW? add fasta fai

  input:
    // set val(meta), file(ref), file(features) from referencesForTranscriptomeExtraction
    set val(meta), file(ref), file(fai), file(features) from referencesForTranscriptomeExtraction
                    // .filter { it[1].size() == 2 } //2 files needed aka skip if fasta-only
                    // .map { meta, files ->
                    //   // files.sort { a,b -> a.name.substring(a.name.lastIndexOf('.')+1) <=> b.name.substring(b.name.lastIndexOf('.')+1) } //ensure gff goes after fasta
                    //   [meta, files[0], files[1]] }
  output:
    set val(outmeta), file(outfile) into extractedTranscriptomes //transcripts4indexing, transcripts4rnfSimReads

  when:
    'rna2rna'.matches(params.mapmode) || 'rna2dna'.matches(params.mapmode)

  shell:
    // println(prettyPrint(toJson(meta)))
    basename=getTagFromMeta(meta)
    outmeta = meta.subMap(['species','version']) //meta.clone()
    outmeta.seqtype = 'RNA'
    outfile = "${basename}.transcripts.fa"
    // println(prettyPrint(toJson(outmeta)))
    // FEATURE_FIELD = meta.featfmt == 'bed' ? 8 : 3 //BED OR GFF3
    // '''
    // gffread --merge -W -w !{outfile} -g !{ref} !{features}
    // '''
    //set -eo pipefail
    '''
    gffread -W -w- -g !{ref} !{features} \
      | awk '/^>/ { if(NR>1) print "";  printf("%s\\t",$0); next; } { printf("%s",$0);} END {printf("\\n");}' \
      | tee tmp.fa \
      | awk 'NR==FNR{all[$1]+=1}; NR!=FNR{if(all[$1]==1){print}}' - tmp.fa  \
      | tr '\\t' '\\n' \
      > !{outfile} && rm tmp.fa
    '''
    // #-w- AND | awk '/^>/ { if(NR>1) print "";  printf("%s\\t",$0); next; } { printf("%s",$0);} END {printf("\\n");}' \
    // #| sort -k1,1V | tr '\\t' '\\n' > !{outfile}
    // '''
    //     #awk '/^>/ { if(NR>1) print "";  printf("%s\\n",$0); next; } { printf("%s",$0);} END {printf("\\n");}' tmp \
    // #| paste - - | sort -k1,1V | tr '\\t' '\\n' > !{outfile}
}

process faidxTranscriptomeFASTA {
  tag("${refmeta}")
  label 'samtools'

  input:
    set val(refmeta), file(fa) from extractedTranscriptomes

  output:
    set val(refmeta), file(fa), file("${fa}.fai") into transcriptomesForIndexing, transcriptomesForRnfSimReads

  script:
  """
  samtools faidx ${fa}
  """
}

/*
Resolve variables emebeded in single-quoted strings
*/
def String resolveScriptVariables(String template, Map binding) {
  def engine = new groovy.text.SimpleTemplateEngine()
  engine.createTemplate(template).make(binding).toString()
}


process indexGenerator {
  label 'index'
  maxForks 1
  container { "${mapper.container}" }
  // tag("${alignermeta.tool} << ${refmeta}")
  tag { [refmeta.subMap(['species','version','seqtype']), mapper.subMap(['tool','version'])] }

  input:
    // set val(alignermeta), val(refmeta), file(ref), file(fai) from aligners.combine(genomesForIndexing.mix(transcriptomesForIndexing))
    set val(mapper), val(refmeta), file(ref), file(fai) from mappersChannel.combine(genomesForIndexing.mix(transcriptomesForIndexing))

  output:
    set val(mapper), val(refmeta), file(ref), file(fai), file("*") into indices


  when: //check if reference intended for {D,R}NA alignment reference and tool has a template declared for that purpose which is also included in mapmode
    [['RNA','rna2rna'],['DNA','rna2dna'],['DNA','dna2dna']].any { refmeta.seqtype == it[0] && mapper.containsKey(it[1]) && it[1].matches(params.mapmode) }

  exec:
    //meta = [toolmodes: alignermeta.modes, tool: "${alignermeta.tool}", target: "${ref}"]+refmeta.subMap(['species','version','seqtype'])
    // meta = [mapper: mapper, target: refmeta+[file: ref]]
    // println(groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(mapper+refmeta+[ref: ref])))
    def binding = [ref: "${ref.name}", task: task.clone()]
    script:
      if('index' in mapper.templates) { //Indexing template expected
        template mapper.index == true ? "index/${mapper.tool}.sh" : "index/${mapper.index}" //either default or explicit template file name
      } else { //indexing script defined in config
        resolveScriptVariables(mapper.index, binding)
      }
}

process rnfSimReads {
  // echo true
  tag{simmeta}
  label 'rnftools'
  label 'slow'

  input:
    // set val(meta), file(ref), file(fai) from referencesWithIndex4rnfSimReads
    set val(meta), file(ref), file(fai) from genomesForRnfSimReads.mix(transcriptomesForRnfSimReads)
    // set val(meta), file(ref) from transcripts4rnfSimReads
    // each nsimreads from params.simreadsDNA.nreads.toString().tokenize(",")*.toInteger()
    each coverage from params.simreadsDNA.coverage
    each length from params.simreadsDNA.length.toString().tokenize(",")*.toInteger()
    each simulator from params.simreadsDNA.simulator
    each mode from params.simreadsDNA.mode //PE, SE
    each distance from params.simreadsDNA.distance //PE only
    each distanceDev from params.simreadsDNA.distanceDev //PE only

  output:
    // set val(simmeta), file("*.fq.gz") into readsForAlignment
    // set val(simmeta), file(ref), file("*.fq.gz") into readsForCoordinateConversion
    set val(simmeta), file(ref), file(simStats), file("*.fq.gz") into simulatedReads


  when:
    !(mode == "PE" && simulator == "CuReSim") && \
    (meta.seqtype == 'RNA' || (meta.seqtype == 'DNA' && 'dna2dna'.matches(params.alnmode) ))
    // ((meta.seqtype == 'mRNA' && 'rna2rna'.matches(params.alnmode)) || (meta.seqtype == 'DNA' && 'rna2rna'.matches(params.alnmode))


  // exec:
  //   println(prettyPrint(toJson(meta)))

  script:
    basename=meta.species+"_"+meta.version+"_"+simulator
    simmeta = meta.subMap(['species','version','seqtype'])+["simulator": simulator, "coverage":coverage, "mode": mode, "length": length, 'coordinates': meta.seqtype]
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
    rnftools.mishmash.sample(\\"${basename}_reads\\",reads_in_tuple=${tuple})
    rnftools.mishmash.${simulator}(
            fasta=\\"${ref}\\",
            coverage=${coverage},
            ${dist}
            ${distDev}
            read_length_1=${len1},
            read_length_2=${len2}
    )
    include: rnftools.include()
    rule: input: rnftools.input()
    " > Snakefile
    snakemake -p \
    && paste --delimiters '=' <(echo -n nreads) <(sed -n '1~4p' *.fq | wc -l) > simStats \
    && time sed -i '2~4 s/[^ACGTUacgtu]/N/g' *.fq \
    && time gzip --fast *.fq \
    && find . -type d -mindepth 2 | xargs rm -r
    """
}

//extract simulation stats from file (currently number of reads only), reshape and split to different channels
readsForCoordinateConversion = Channel.create()
simulatedReads.map { simmeta, ref, simStats, simReads ->
    // simStats.splitEachLine("=", { record ->
    //   if(record.size() > 1) {
    //     v = record[1]
    //     simmeta."${record[0]}" = v.isInteger() ? v.toInteger() : v.isDouble() ? v.toDouble() : v
    //   }
    // })
    parseFileMap(simStats, simmeta)
    new Tuple(simmeta, ref, simReads)
  }
  .tap ( readsForCoordinateConversion )
  .map { simmeta, ref, simReads  ->
    new Tuple(simmeta, simReads)
  }
  // .view { it -> println(groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))) }
  .set{ readsForAlignment }




// // process simStats{
// //   input:
// //     set val(simmeta), file(reads) from readsForSimStats

// //   output:
// //     set val(simmeta), stdout(count) into simCounts

// //   """
// //   zcat ${reads[0]} | sed -n '1~4p' | wc -l
// //   """
// // }

process convertReadCoordinates {
  label 'groovy'
  echo true
  tag{simmeta.subMap(['species','version'])}


  input:
    set val(simmeta), file(ref), file(reads), val(refmeta), file(fai) from readsForCoordinateConversion.combine(genomeIndicesForReadCoordinateConversion)

  output:
    set val(outmeta), file('*.fq.gz') into convertedCoordinatesReads

  when:
    simmeta.seqtype == 'RNA' && 'rna2dna'.matches(params.alnmode) \
    && simmeta.species == refmeta.species && simmeta.version == refmeta.version

  // exec:
  //   println(prettyPrint(toJson(simmeta)))
  //   println(prettyPrint(toJson(refmeta)))

  script:
  out1 = reads[0].name.replace('.1.fq.gz','.R1.fq.gz')
  out2 = reads[1].name.replace('.2.fq.gz','.R2.fq.gz')
  outmeta = [:]
  outmeta.putAll(simmeta)
  outmeta.remove('coordiantes')
  outmeta.coordinates = 'DNA'
  """
  tct_rnf.groovy \
    --genome-index ${fai} \
    --transcriptome ${ref} \
    --in-forward ${reads[0]} --in-reverse ${reads[1]} \
    --out-forward ${out1} --out-reverse ${out2}
  """
}

// // // convertedCoordinatesReads.view()

// // // convertedCoordinatesReads.mix(readsForAlignment).combine(indices).combine(alignersParams).view { it -> println(groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it)))}

process mapSimulatedReads {
  label 'align'
  container { "${mapper.container}" }
  // tag {"${mapper.tool}@${mapper.version} ${refmeta.subMap(['species','version','seqtype'])} << ${simmeta.subMap(['simulator','seqtype'])} @ ${mode} params: ${paramsmeta.label}"}

  input:
    set val(simmeta), file(reads), val(mapper), val(refmeta), file(ref), file(fai), file('*'), val(paramsmeta), val(mode) from convertedCoordinatesReads.mix(readsForAlignment).combine(indices).combine(mappersParamsChannel).combine(mapModesChannel)

  // output:
  //   set val(mapmeta), file(ref), file(fai), file('*.?am'), file('.command.trace') into alignedSimulated

  when: //only align simulated reads to the corresponding genome, using the corresponding params set, in the correct mode: DNA2DNA, RNA2DNA, RNA2RNA
    mapper.tool == paramsmeta.tool && mapper.version.matches(paramsmeta.version.join('|')) \
    && mapper.containsKey(mode) && mode.matches(paramsmeta.mode) \
    && refmeta.species == simmeta.species && refmeta.version == simmeta.version \
    && mode.startsWith(simmeta.seqtype.toLowerCase()) \
    && mode.endsWith(simmeta.coordinates.toLowerCase()) \
    && mode.endsWith(refmeta.seqtype.toLowerCase())

    //&& paramsmeta.mode.endsWith(simmeta.coordinates.toLowerCase()) &&  paramsmeta.mode.endsWith(refmeta.seqtype.toLowerCase()))
  //   // (paramsmeta.alignMode.startsWith(simmeta.seqtype) && paramsmeta.alignMode.endsWith(simmeta.coordinates) &&  paramsmeta.alignMode.endsWith(idxmeta.seqtype))

  // exec:
  //   println "Reads size: "+reads[0].size()
  //   println "Ref size: "+ref.size()
  // println mapper.version.matches(paramsmeta.version.join('|'))
    // mapmeta = [tool: [name: mapper.tool, version: mapper.version, container: task.container],
    //              target: refmeta,
    //              query: simmeta,
    //              params: paramsmeta.params,
    //              mode: mode]
    // println(groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(mapmeta)))
    // println(groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson([simmeta, mapper, refmeta, paramsmeta, [mode:mode]])))


  // println(prettyPrint(toJson(simmeta))+'\n'+prettyPrint(toJson(idxmeta))+'\n'+prettyPrint(toJson(paramsmeta)))
  // script:
  script:
    ALIGN_PARAMS = paramsmeta.params
    def binding = [ref: ref, reads: reads, task: task.clone(), ALIGN_PARAMS: paramsmeta.params]
    if(mode in mapper.templates) { //Indexing template expected
      template mapper."${mode}" == true ? "${mode}/${mapper.tool}.sh" : "${mode}/${mapper.${mode}}" //either default or explicit template file name
    } else { //indexing script defined in config
      resolveScriptVariables(mapper."${mode}", binding)
    }
}

// // // alignedSimulated.view { it -> println(groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it)))}



// process evaluateAlignmentsRNF {
//   label 'samtools'
//   // label 'ES'
//   tag{alignmeta.tool.subMap(['name'])+alignmeta.target.subMap(['species','version'])+alignmeta.query.subMap(['seqtype','nreads'])+alignmeta.params.subMap(['paramslabel'])}
//   // tag{alignmeta.params.subMap(['paramslabel'])}
//   // tag{alignmeta.subMap(['tool','simulator','target.species','alignMode','paramslabel'])}

//   input:
//     set val(alignmeta), file(ref), file(fai), file(samOrBam) from alignedSimulated.map { meta, ref, fai, samOrBam, trace ->
//         meta.trace = [:]
//         parseFileMap(trace, meta.trace) //could be parseFileMap(trace, meta.trace, 'realtime') or parseFileMap(trace, meta, ['realtime','..']) or parseFileMap(trace, meta) to capture all fields
//         new Tuple(meta, ref, fai, samOrBam)
//       }

//   output:
//      set val(alignmeta), file ('*.json') into evaluatedAlignmentsRNF
//     //  set val(alignmeta), file('ES.gz'),  into esChannel  //add to script: --es-output ES.gz

//   // exec:
//   script:
//   // println prettyPrint(toJson(alignmeta))
//   // println alignmeta.inspect()

//   """
//   set -eo pipefail
//   samtools view ${samOrBam} \
//   | eval_rnf.groovy \
//       --allowed-delta ${params.allowedDelta} \
//       --faidx ${fai} \
//       --output summary.json \
//   """
// }


// /**
// 1. Embed evaluation results JSON in META JSON.
// 2. Collect all datapoints in one JSON for output.
// **/
// def slurper = new groovy.json.JsonSlurper()
// evaluatedAlignmentsRNF.map { META, JSON ->
//       META.evaluation = slurper.parseText(JSON.text)
//       META
//   }
//   .collect()
//   .map {
//     file("${params.outdir}").mkdirs()
//     outfile = file("${params.outdir}/allstats.json")
//     outfile.text = groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it.sort( {k1,k2 -> k1.tool.name <=>  k2.tool.name} ) ))
//   }


// // process plotSummarySimulated {
// //   label 'rscript'
// //   label 'figures'

// //   input:
// //     // set file(csv), file(json) from collatedSummariesSimulatedDNA
// //     set file(json), file(categories) from collatedSummariesSimulated

// //   output:
// //     file '*' into collatedSummariesPlotsSimulated

// //   shell:
// //   '''
// //   < !{json} plot_simulatedDNA.R
// //   '''
// // }



// // process plotSummarySimulated {
// //   label 'rscript'
// //   label 'figures'

// //   input:
// //     // set file(csv), file(json) from collatedSummariesSimulatedDNA
// //     set file(json), file(categories) from collatedSummariesSimulatedDNA

// //   output:
// //     file '*' into collatedSummariesPlotsSimulated

// //   shell:
// //   '''
// //   < !{json} plot_simulatedDNA.R
// //   '''
// // }

// // //WRAP-UP
// // writing = Channel.fromPath("$baseDir/report/*").mix(Channel.fromPath("$baseDir/manuscript/*")) //manuscript dir exists only on manuscript branch

// // process render {
// //   tag {"Render ${Rmd}"}
// //   label 'rrender'
// //   label 'report'
// //   stageInMode 'copy'
// //   //scratch = true //hack, otherwise -profile singularity (with automounts) fails with FATAL:   container creation failed: unabled to {task.workDir} to mount list: destination ${task.workDir} is already in the mount point list

// //   input:
// //     // file('*') from plots.flatten().toList()
// //     // file('*') from plotsRealRNA.flatten().toList()
// //     file(Rmd) from writing
// //     file('*') from collatedDetailsPlotsSimulatedDNA.collect()
// //     file('*') from collatedSummariesPlotsSimulatedDNA.collect()

// //   output:
// //     file '*'

// //   script:
// //   """
// //   #!/usr/bin/env Rscript

// //   library(rmarkdown)
// //   library(rticles)
// //   library(bookdown)

// //   rmarkdown::render("${Rmd}")
// //   """
// // }
// // }

// // //WRAP-UP
// // writing = Channel.fromPath("$baseDir/report/*").mix(Channel.fromPath("$baseDir/manuscript/*")) //manuscript dir exists only on manuscript branch

// // process render {
// //   tag {"Render ${Rmd}"}
// //   label 'rrender'
// //   label 'report'
// //   stageInMode 'copy'
// //   //scratch = true //hack, otherwise -profile singularity (with automounts) fails with FATAL:   container creation failed: unabled to {task.workDir} to mount list: destination ${task.workDir} is already in the mount point list

// //   input:
// //     // file('*') from plots.flatten().toList()
// //     // file('*') from plotsRealRNA.flatten().toList()
// //     file(Rmd) from writing
// //     file('*') from collatedDetailsPlotsSimulatedDNA.collect()
// //     file('*') from collatedSummariesPlotsSimulatedDNA.collect()

// //   output:
// //     file '*'

// //   script:
// //   """
// //   #!/usr/bin/env Rscript

// //   library(rmarkdown)
// //   library(rticles)
// //   library(bookdown)

// //   rmarkdown::render("${Rmd}")
// //   """
// // }
