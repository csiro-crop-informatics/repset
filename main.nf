#!/usr/bin/env nextflow

//For pretty-printing nested maps etc
import groovy.json.JsonGenerator
import groovy.json.JsonSlurper
import groovy.json.JsonOutput
//as JsonGenerator
// import static groovy.json.JsonGenerator.*

//various helper methods in lib/Helpers.groovy
def hlp = new Helpers()

//Otherwise JSON generation triggers stackoverflow when encountering Path objects
// jsonGenerator = new JsonGenerator.Options()
//                 .addConverter(java.nio.file.Path) { java.nio.file.Path p, String key -> p.toUriString() }
//                 .build()

//Preventing stack overflow on Path objects and other  when map -> JSON
JsonGenerator jsonGenerator = new JsonGenerator.Options()
                .addConverter(java.nio.file.Path) { java.nio.file.Path p, String key -> p.toUriString() }
                .addConverter(Duration) { Duration d, String key -> d.durationInMillis }
                .addConverter(java.time.OffsetDateTime) { java.time.OffsetDateTime dt, String key -> dt.toString() }
                .addConverter(nextflow.NextflowMeta) { nextflow.NextflowMeta m, String key -> m.toJsonMap() }  //incompatible with Nextflow <= 19.04.0
                .excludeFieldsByType(java.lang.Class) // .excludeFieldsByName('class')
                // .excludeNulls()
                .build()

//Input validation specified elswhere
def validators = new Validators() //from lib/ instead of new GroovyShell().parse(new File("${baseDir}/groovy/Validators.groovy"))

//Read, parse, validate and sanitize alignment/mapping tools config
def allRequired = ['tool','version','container','index'] //Fields required for each tool in config
def allOptional = ['versionCall']
def allModes = 'dna2dna|rna2rna|rna2dna' //At leas one mode has to be defined as supported by each tool
def allVersions = validators.validateMappersDefinitions(log, params.mappersDefinitions, allRequired, allOptional, allModes)


//Check if specified template files exist
validators.validateTemplatesAndScripts(params.mappersDefinitions, (['index']+(allModes.split('\\|') as List)), "${baseDir}/templates")

//Read, sanitize and validate alignment/mapping param sets
validators.validateMapperParamsDefinitions(params.mapperParamsDefinitions, allVersions, allModes)

//Parse, sanitize and validate input dataset definitions
def requiredInputFields = ['species','version','fasta','seqtype']
validators.validateInputDefinitions(params.references, requiredInputFields, ['gff'])

//Some of the real_reads spec may point to sra ids, other to local files
Channel.from(params.rreads)
  .branch {
    sra   : it.containsKey('sra') && !(it.containsKey('r1') || it.containsKey('r2'))
    path :  !(it.containsKey('sra')) && (it.containsKey('r1') || it.containsKey('r2')) //only considering paired
    sink  : true //anything else goes into a black hole
  }.set { realReadsDefinitionsChannel }

realReadsDefinitionsChannel.sink.subscribe {
   log.warn "Malformed real reads entry will be ignored: ${it}"
}

//Prepare real read def for downloading: duplicate entry if multiuple SRR specified and for R1 and R2
realReadsDefinitionsChannel.sra
  .flatMap { it ->
    if(it.sra instanceof List) { //multiple SRR entries
      def elems = []
      it.sra.each { currentSRA ->
        elems << [sra: currentSRA]+it.subMap(it.keySet()-'sra')
      }
      elems
    } else {
      [it] //flatMap() will unwrap this
    }
  }
  .combine(Channel.from([1,2])) //download each mate sepearately
  .set { srrDownloadChannel }

// .path.view() //LOCAL READS NOT USED YET
// .sink.view {

process asperaDownload {
  tag { "${SRR} R${MATE}" }
  label 'aspera'
  storeDir { executor == 'awsbatch' ? null : "downloaded" }

  input:
    tuple val(META), val(MATE) from srrDownloadChannel

  output:
    tuple val(META), file("${SRR}_${MATE}.fastq.gz") into SraDownloadsChannel

  script:
  SRR = META.sra
  """
  ascp -T --policy=fair -P33001 -i /home/aspera/.aspera/cli/etc/asperaweb_id_dsa.openssh \
    era-fasp@fasp.sra.ebi.ac.uk:vol1/fastq/${SRR[0..5]}/${SRR}/${SRR}_${MATE}.fastq.gz ./
  """
}
SraDownloadsChannel
.groupTuple().view()


if(params.justvalidate) {
  log.info "Finished validating input config, exiting. Run without --justvalidate to proceed further."
  System.exit 0
}

//Validated now, so gobble up mappers...
// mappersChannel = Channel.from(params.mappersDefinitions)
Channel.from(params.mappersDefinitions)
  .filter{ params.mappers == 'all' || it.tool.matches(params.mappers) } //TODO Could allow :version
  // .tap { mappersMapChannel }
  // .map { it.subMap(allRequired)} //Exclude mapping specific fields from indexing process to avoid re-indexing e.g. on changes made to a mapping template
  .set { mappersChannel }
  // .into { mappersIdxChannel; mappersMapChannel; mappersVersionChannel }

// mappersVersionChannel.view{ it -> JsonOutput.prettyPrint(jsonGenerator.toJson(it))}

mappersChannel.filter {
  if(it.containsKey('versionCall')) {
    true
  } else {
    log.warn """
    versionCall not specified for ${it.tool} ${it.version}
    it will not be included in this run
    """
  }
}
.set { mappersVithVersionCallChannel }

process parseMapperVersion {
  container { "${mapmeta.container}" }
  tag { mapmeta.subMap(['tool','version']) }

  input:  val(mapmeta) from mappersVithVersionCallChannel

  output: tuple val(mapmeta), stdout into mappersCapturedVersionChannel

  script: "${mapmeta.versionCall}"
  // script: "set -o pipefail; ${mapmeta.versionCall}"
}

mappersCapturedVersionChannel
.map { meta, ver ->
  if(meta.version != ver.trim()) {
    log.warn """
    Decalred version ${meta.version} for ${meta.tool}
    does not match version ${ver.trim()}
    obtained from versionCall: ${meta.versionCall}
    Updating version in metadata to ${ver.trim()}
    """
    //Please correct your mapper configuration file(s).
    // throw new RuntimeException('msg') or //
    // session.abort(new Exception())
    meta.version = ver.trim()
    if(!meta.container.contains(meta.version)) {
      log.error """
      Updated tool version string ${meta.version}
      not found in container image spec ${meta.container}.
      Please correct your mapper configuration file(s).

      Aborting...
      """
      session.abort(new Exception())  // throw new RuntimeException('msg')
    }
  }
  meta
}
// .view{ it -> JsonOutput.prettyPrint(jsonGenerator.toJson(it))}
.into { mappersIdxChannel; mappersMapChannel }

//...and their params definitions
mappersParamsChannel = Channel.from(params.mapperParamsDefinitions)

//one or more mapping mode
mapModesChannel = Channel.from(params.mapmode.split('\\||,'))

/*
 * Add to or overwrite map content recursively
 * Used to enable the use of NF -params-file opt such that params can be added and not just overwritten
 */
Map.metaClass.addNested = { Map rhs ->
    def lhs = delegate
    rhs.each { k, v -> lhs[k] = lhs[k] in Map ? lhs[k].addNested(v) : v }
    lhs
}

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

/*
  Simplistic method for checking if String is URL
*/
String.metaClass.isURL() {
   delegate.matches("^(https?|ftp)://.*\$")
}

def helpMessage() {
  log.info"""
  Usage:

  nextflow run csiro-crop-informatics/repset -profile singularity
  nextflow run csiro-crop-informatics/repset -profile docker

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
Channel.from(params.references)
.take( params.subset ) //only process n data sets (-1 means all)
.combine(Channel.from('fasta','gff')) //duplicate each reference record
.filter { meta, fileType -> meta.containsKey(fileType)} //exclude gff record if no gff declared
.tap { refsToStage } //download if URL
.filter { meta, fileType ->  !(meta."${fileType}").isURL() } //Exclude URLs
.map { meta, fileType ->  [meta, fileType, file(meta."${fileType}")] } //file declaration required for correct binding of source path
// .view { it -> groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))}
.set { refsToStageLocal }


process stageRemoteInputFile {
  tag{meta.subMap(['species','version'])+[fileType: fileType]}
  errorStrategy 'finish'
  /*
  storeDir can be problematic on s3 - leads to "Missing output file(s)" error
  workDir should be more robust as it is mounted in singularity unlike outdir?
  */
  storeDir { executor == 'awsbatch' ? null : "downloaded" }


  input:
    set val(meta), val(fileType) from refsToStage //fastaChn.mix(gffChn)

  output:
    set val(outmeta), file(outfile)  into stagedFilesRemote

  when:
    (meta."${fileType}").isURL()

  script:
    basename=getTagFromMeta(meta)
    outfile =  "${basename}.${fileType}"
    outmeta = meta.subMap(['species', 'version','seqtype'])
    fpath = meta."${fileType}"
    decompress = fpath.matches("^.*\\.gz\$") ?  "| gunzip --stdout " :  " "
    """curl ${fpath} ${decompress} > ${outfile}"""
}

process stageLocalInputFile {
  tag{meta.subMap(['species','version'])+[fileType: fileType]}
  errorStrategy 'finish'

  input:
    set val(meta), val(fileType), file(infile) from refsToStageLocal

  output:
    set val(outmeta), file(outfile)  into stagedFilesLocal

  script:
    basename=getTagFromMeta(meta)
    outfile = "${basename}.${fileType}"
    outmeta = meta.subMap(['species', 'version','seqtype'])
    if((infile.name).matches("^.*\\.gz\$")){ //GZIPPED
      """gunzip --stdout  ${infile}  > ${outfile} """
    } else { //FLAT
      """cp -s  ${infile} ${outfile}"""
    }
}

// referencesOnly = Channel.create()
// referencesForTranscriptomeExtraction = Channel.create()
stagedFilesRemote.mix(stagedFilesLocal)
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


process extractTranscripts {
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

mappersIdxChannel.combine(genomesForIndexing.mix(transcriptomesForIndexing))
.filter { mapper, refmeta, ref, fai->
  [['RNA','rna2rna'],['DNA','rna2dna'],['DNA','dna2dna']].any { refmeta.seqtype == it[0] && mapper.containsKey(it[1]) && it[1].matches(params.mapmode) }
}
// .view{ it -> groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))}
.map { mapper, refmeta, ref, fai ->
  [
    mapper.subMap(allRequired)+[idxTemplate: ('index' in mapper.templates) ], //second part shuld have been done at validation
    refmeta,
    ref,
    fai
  ]
} //Exclude mapping specific fields from indexing process to avoid re-indexing e.g. on changes made to a mapping template
// .view{ it -> groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))}
.set { forIndexing }

process indexGenerator {
  label 'index'
  container { "${mapper.container}" }
  // tag("${alignermeta.tool} << ${refmeta}")
  tag { [refmeta.subMap(['species','version','seqtype']), mapper.subMap(['tool','version'])] }

  input:
    // set val(alignermeta), val(refmeta), file(ref), file(fai) from aligners.combine(genomesForIndexing.mix(transcriptomesForIndexing))
    // set val(mapper), val(refmeta), file(ref), file(fai) from mappersIdxChannel.combine(genomesForIndexing.mix(transcriptomesForIndexing))
    set val(mapper), val(refmeta), file(ref), file(fai) from forIndexing

  output:
    set val(idxmeta), file(ref), file(fai), file("*") into indices

  // when: //check if reference intended for {D,R}NA alignment reference and tool has a template declared for that purpose which is also included in mapmode
  //   [['RNA','rna2rna'],['DNA','rna2dna'],['DNA','dna2dna']].any { refmeta.seqtype == it[0] && mapper.containsKey(it[1]) && it[1].matches(params.mapmode) }

  exec:
    //meta = [toolmodes: alignermeta.modes, tool: "${alignermeta.tool}", target: "${ref}"]+refmeta.subMap(['species','version','seqtype'])
    // meta = [mapper: mapper, target: refmeta+[file: ref]]
    // println(groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(mapper+refmeta+[ref: ref])))
    def binding = [ref: "${ref.name}", task: task.clone()]
    idxmeta = [mapper: mapper.subMap(['tool','version']), reference: refmeta]
  script:
    if(mapper.idxTemplate == true) { //Indexing template file declared
      template mapper.index == true ? "index/${mapper.tool}.sh" : "index/${mapper.index}" //either default or explicit template file name
    } else { //indexing script string embeded in config
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
    (meta.seqtype == 'RNA' || (meta.seqtype == 'DNA' && 'dna2dna'.matches(params.mapmode) ))
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
    set -eo pipefail
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
    && awk 'END{print "nreads="NR/4}' *.fq > simStats \
    && for f in *.fq; do
        paste - - - - < \$f \
        | awk -vFS="\\t" -vOFS="\\n" '{gsub(/[^ACGTUacgtu]/,"N",\$2);print}' \
        | gzip -c > \${f}.gz
    done && rm *.fq \
    && find . -type d -mindepth 2 | xargs rm -r
    """
}
  //  && paste --delimiters '=' <(echo -n nreads) <(sed -n '1~4p' *.fq | wc -l) > simStats \
  //  && time sed -i '2~4 s/[^ACGTUacgtu]/N/g' *.fq \
// && time gzip --fast *.fq \

//extract simulation stats from file (currently number of reads only), reshape and split to different channels
// readsForCoordinateConversion = Channel.create()
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
  .tap { readsForCoordinateConversion }
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
    simmeta.seqtype == 'RNA' && 'rna2dna'.matches(params.mapmode) \
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

// mappersMapChannel  .view { it -> groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))}

/*
 This is where we combine and match (and filter out inapropriate combinations)
 * simulated read sets
 * mapper indices
 * mapper specs
 * mapper params sets
 * mapping modes (dna2dna, rna2dna, rna2rna)
*/
convertedCoordinatesReads.mix(readsForAlignment)
.combine(indices)
.filter { simmeta, reads, idxmeta, ref, fai, idx ->
  //simulated reads vs reference species & version must match
  idxmeta.reference.species == simmeta.species && idxmeta.reference.version == simmeta.version
}
.combine(mappersMapChannel) //mappers definitions
.combine(mappersParamsChannel) //mappers params definitions
.filter { simmeta, reads, idxmeta, ref, fai, idx, mapper, paramsmeta -> //tool & version must match between mapper and a params set
  [mapper.tool, idxmeta.mapper.tool].every { it == paramsmeta.tool }  \
  && mapper.version == idxmeta.mapper.version \
  && mapper.version in paramsmeta.version
}
.combine(mapModesChannel) //one or more of the 3 possible mapping modes
.filter { simmeta, reads, idxmeta, ref, fai, idx, mapper, paramsmeta, mode ->  //map mode check is mapper able to work in that mode mand is params set aimed at this mode
  mapper.containsKey(mode) && mode.matches(paramsmeta.mode) \
  && mode.startsWith(simmeta.seqtype.toLowerCase()) \
  && [simmeta.coordinates, idxmeta.reference.seqtype].every { mode.endsWith( it.toLowerCase() ) }
}
// .view{ groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))}
.map { simmeta, reads, idxmeta, ref, fai, idx, mapper, paramsmeta, mode ->
  def template = (mode in mapper.templates) ? (mapper."${mode}" == true ? "${mode}/${mapper.tool}.sh" : "${mode}/${mapper.${mode}}") : false;
  [
    [mapper: mapper.subMap(['tool','version','container']), query: simmeta, target: idxmeta.reference, params: paramsmeta.subMap(['label','params'])],
    reads, //as is
    ref, fai, idx, //as is
    [template: template, script: (template ? false: mapper."${mode}")],
    paramsmeta.params
  ]
}
// .view{ groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))}
// .count()
.set{ combinedToMap }


// .view{ groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))}

process mapSimulatedReads {
  label 'align'
  container { "${meta.mapper.container}" }
  tag {"${meta.target.seqtype}@${meta.target.species}@${meta.target.version} << ${meta.query.nreads}@${meta.query.seqtype}; ${meta.mapper.tool}@${meta.mapper.version}@${meta.params.label}"}
  // beforeScript meta.mapper.containsKey('versionCall') ? "${meta.mapper.versionCall} > .mapper.version" : ''

  input:
    set val(meta), file(reads), file(ref), file(fai), file('*'), val(run), val(ALIGN_PARAMS) from combinedToMap

  output:
    set val(meta), file(ref), file(fai), file('*.?am'), file('.command.trace') into alignedSimulated

  script:
    def binding = [ref: ref, reads: reads, task: task.clone(), ALIGN_PARAMS: ALIGN_PARAMS]
    meta.resources = task.subMap(['cpus','memory','time'])
    if(run.template) { //if template file specified / declared
      template run.template //either default or explicit template file name
    } else { //indexing script defined in config
      resolveScriptVariables(run.script, binding)
    }
}

process evaluateAlignmentsRNF {
  label 'groovy_samtools'
  // label 'ES'
  // tag{alignmeta.tool.subMap(['name'])+alignmeta.target.subMap(['species','version'])+alignmeta.query.subMap(['seqtype','nreads'])+alignmeta.params.subMap(['paramslabel'])}
  // tag{alignmeta.params.subMap(['paramslabel'])}
  // tag{alignmeta.subMap(['tool','simulator','target.species','alignMode','paramslabel'])}
  tag {"${meta.target.seqtype}@${meta.target.species}@${meta.target.version} << ${meta.query.nreads}@${meta.query.seqtype}; ${meta.mapper.tool}@${meta.mapper.version}@${meta.params.label}"}

  input:
    set val(meta), file(ref), file(fai), file(samOrBam) from alignedSimulated.map { mapmeta, ref, fai, samOrBam, trace ->
        def traceMap = [:]
        parseFileMap(trace, traceMap) //could be parseFileMap(trace, meta.trace, 'realtime') or parseFileMap(trace, meta, ['realtime','..']) or parseFileMap(trace, meta) to capture all fields
        [mapmeta+[trace: traceMap], ref, fai, samOrBam]
      }

  output:
     set val(meta), file ('*.json') into evaluatedAlignmentsRNF
    //  set val(alignmeta), file('ES.gz'),  into esChannel  //add to script: --es-output ES.gz

  // exec:
  script:
  // println prettyPrint(toJson(alignmeta))
  // println alignmeta.inspect()
  // println groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(meta))
  """
  set -eo pipefail
  samtools view ${samOrBam} \
  | eval_rnf.groovy \
      --allowed-delta ${params.allowedDelta} \
      --faidx ${fai} \
      --output summary.json \
  """
}


/**
1. Embed evaluation results JSON in META JSON.
2. Collect all datapoints in one JSON for output.
**/
def slurper = new JsonSlurper()
evaluatedAlignmentsRNF.map { META, JSON ->
      [ META + [evaluation: slurper.parseText(JSON.text)] ]
  }
  .collect()
  .map { //Generate:  [ allstats.json, runmetapart.json]
    [
      JsonOutput.prettyPrint(jsonGenerator.toJson(
        it.sort( { k1,k2 -> k1.mapper.tool <=>  k2.mapper.tool }
      ))),
      JsonOutput.prettyPrint(jsonGenerator.toJson([
        workflow : workflow.getProperties(),
        params   : params
      ]))
    ]
  }
  .set { resultsJsonChannel }

// //WRAP-UP --TODO Manuscript rendering to be separated
// // writing = Channel.fromPath("$baseDir/report/*.Rmd").mix(Channel.fromPath("$baseDir/manuscript/*")) //manuscript dir exists only on manuscript branch

process renderReport {
  tag {"Render ${Rmd}"}
  label 'rrender'
  label 'report'
  stageInMode 'copy'
  cache false //Input includes run metadata so cache would not work anyway

  input:
    file(Rmd) from Channel.fromPath("$baseDir/report/report.Rmd")
    tuple file('allstats.json'), file('runmetapart.json') from resultsJsonChannel

  output:
    file '*'

  when:
    !(workflow.profile.contains('CI')) //until leaner container

  script:
  """
  #!/usr/bin/env Rscript

  library(rmarkdown)
  library(rticles)
  library(bookdown)
  library(tidyverse)
  library(jsonlite)
  library(kableExtra)

  rmarkdown::render("${Rmd}")
  """
}


workflow.onComplete {
    log.info "Pipeline complete, collecting metadata..."

    def runmeta = [:]
    runmeta['os']  = [:]
    runmeta['os']['Architecture']  = System.getProperty("os.arch")
    runmeta['os']['Name']  = System.getProperty("os.name")
    runmeta['os']['Version']  = System.getProperty("os.version")
    runmeta['java']  = [:]
    runmeta['java']['VM name']  = System.getProperty("java.vm.name")
    runmeta['java']['VM version']  = System.getProperty("java.vm.version")
    runmeta['java']['Vendor name']  = System.getProperty("java.vendor")
    runmeta['java']['Runtime Environment Version']  = System.getProperty("java.runtime.version")
    runmeta['java']['Version']  = System.getProperty("java.version")

    runmeta['workflow'] = workflow.getProperties()
    runmeta['params'] = params

    // println  !(runmeta['workflow']['nextflow'] instanceof java.util.LinkedHashMap)

    def runmetaJSON = new File("${params.infodir}/runmeta.json")
    runmetaJSON.text = JsonOutput.prettyPrint(jsonGenerator.toJson(runmeta))

    //  def command = "tree -h -D all"
    // def proc = command.execute()
    // proc.waitFor()

    // println "Process exit code: ${proc.exitValue()}"
    // println "Std Err: ${proc.err.text}"
    // println "Std Out: ${proc.in.text}"
    // evaluate(new File("$baseDir/conf/ApiCalls.groovy"))
    // apiCalls = new ApiCalls()


  // workflow.repository = 'rsuchecki/repset'
  // workflow.commitId = 'c9ac00dfa9b67b00e00a9bc71063b6bc76675d36'
  // workflow.revision = 'master'
  // IF --release requested by the user and execution from GH repo
  if(params.release && workflow.repository && workflow.commitId && workflow.revision && workflow.scriptName == 'main.nf') {
    //if(workflow.revision ==~ /^v?([0-9]+)\.([0-9]+)\.?([0-9]+)?$/ ) {
      GroovyShell shell = new GroovyShell()
      def apiCalls = shell.parse(new File("$baseDir/groovy/ApiCalls.groovy"))

      // def instant = Instant.now()
      // println instant
      // def utc =  LocalDateTime.ofInstant(instant, ZoneOffset.UTC)
      // def local = LocalDateTime.ofInstant(instant, ZoneId.systemDefault())
      def formatter = java.time.format.DateTimeFormatter.ofPattern("yyyy-MM-dd'T'HHmmssx");

      releaseArgs = [
        REPO : workflow.repository.replaceFirst("^(http[s]?://github\\.com/|git@github\\.com:)","").replaceFirst("\\.git\$",""),
        COMMIT : workflow.commitId,
        LOCAL_FILES : [
          // "${params.outdir}/report.html",
          // "${params.outdir}/repset-manuscript.pdf",
          "${params.outdir}/report.html",
          "${params.outdir}/allstats.csv",
          "${params.outdir}/allstats.json",
          "${params.infodir}/runmeta.json",
          "${params.infodir}/trace.tsv"
        ],
        RELEASE_TAG: "${workflow.revision}_${workflow.complete.format(formatter)}_${workflow.runName}_${workflow.sessionId}",
        RELEASE_NAME: "${workflow.revision} - results and metadata for run '${workflow.runName}'",
        RELEASE_BODY:
"""
- revision          `${workflow.revision}`
- commit ID          ${workflow.commitId}
- session ID        `${workflow.sessionId}`
- profile           `${workflow.profile}`
- started at        `${workflow.start}`
- completed at      `${workflow.complete}`

see assets for more details.
""".replace('\n','<br />')
      ]
      apiCalls.gitHubRelease(log, releaseArgs, params.draft)
    //} else {
    //  log.warn "Automated GH release generation only aimed at semantically tagged revisions (e.g. v1.5.4), current revision: ${workflow.revision}"
    //  log.warn "Note that this restriction can be lifted without adversely affecting functionality"
    //}
  }
}
