aligners = Channel.from(['bbmap', 'biokanga','dart','gsnap','hisat2','star','subread']) //hera? mapsplice2? Subread
// aligners = Channel.from(['biokanga','gsnap','star']) //hera? mapsplice2? Subread
//datasets = Channel.from(['human_t1r1','human_t1r2','human_t1r3','human_t2r1','human_t2r2','human_t2r3','human_t3r1','human_t3r2','human_t3r3'])
datasets = Channel.from(['human_t1r1','human_t2r1','human_t3r1']).filter{ !params.debug || it == 'human_t2r1' } //Only one ds if debug
url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz'

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

  input:
    val(url)

  output:
    file('*') into refs

  script:
  """
  wget ${url}
  """
}

process downloadDatasets {
  label 'download'
  tag("${dataset}")
  // storeDir "${workflow.workDir}/downloaded" //and put the downloaded datasets there and prevent generating cost to dataset creators through repeated downloads

  input:
    val(dataset) from datasets

  output:
    set val("${dataset}"), file("${dataset}.tar.bz2") into downloadedDatasets

  script:
    """
    wget http://bp1.s3.amazonaws.com/${dataset}.tar.bz2
    """
}



process extractDatasets {
  tag("${dataset}")

  input:
    set val(dataset), file("${dataset}.tar.bz2") from downloadedDatasets

  output:
    set val(dataset), file("${dataset}") into extractedDatasets

  script:
    """
    mkdir -p ${dataset}
    pbzip2 --decompress --stdout -p${task.cpus} ${dataset}.tar.bz2 | tar -x --directory ${dataset}
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
      | cut -f1 > debug.ids
    fgrep --no-filename --no-group-separator -A1 -wf debug.ids ${dataDir}/*.forward.fa \
      | paste - - | sort -k1,1V --parallel ${task.cpus} \
      | awk -vOFS="\\n" '{gsub(/[0-9]+/,++n,\$1);print \$1,\$2}' > r1
    fgrep --no-filename --no-group-separator -A1 -wf debug.ids ${dataDir}/*.reverse.fa \
      | paste - - | sort -k1,1V --parallel ${task.cpus} \
      | awk -vOFS="\\n" '{gsub(/[0-9]+/,++n,\$1);print \$1,\$2}' > r2
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
    !params.debug //omitting this process to speed things up a bit for debug runs

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
    samtools view -F2304  ${sam} | grep -v '^@' | sort -t'.' -k2,2n --parallel ${task.cpus} > sortedsam
    """
}

process fixSAM {
  label 'benchmark'
  label 'slow'
  tag("${meta}")

  input:
    set val(meta), file(sortedsam), file(cig) from sortedSAMs

  output:
    set val(meta), file(fixedsam), file(cig) into fixedSAMs

  script:
  if(params.debug) {
    //1. should probably get exect value not just for --debug run. Otherwise --nummer fixed to 10 mil (?!)
    //2. awk reading sam twice to exclude multi-aligners
    // """
    // fix_sam.rb \
    //   --nummer \$(paste - - < ${cig} | wc -l) \
    //   <(awk 'NR==FNR{occurances[\$1]+=1};NR!=FNR && occurances[\$1]==1{print}' ${sortedsam} ${sortedsam} | tee dedupedsam ) \
    //   > fixedsam
    // """
    """
    fix_sam.rb --nummer \$(paste - - < ${cig} | wc -l) ${sortedsam} > fixedsam
    """
  } else {
    """
    fix_sam.rb <(awk 'NR==FNR{occurances[\$1]+=1};NR!=FNR && occurances[\$1]==1{print}' ${sortedsam} ${sortedsam}) > fixedsam
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
      compare2truth_multi_mappers.rb ${cig} ${fixedsam} > stat
      """
    } else {
      """
      compare2truth.rb ${cig} ${fixedsam} > stat
      """
    }
}

process tidyStats {
  label 'rscript'
  // label 'stats'
  // echo true
  tag("${meta}")

  input:
    set val(meta), file(instats) from stats

  output:
    file 'tidy.csv' into tidyStats

  exec:
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

process aggregateStats {
  // label 'rscript'
  label 'stats'

  input:
    file '*.csv' from tidyStats.collect()

  output:
    file 'all.csv' into aggregatedStats

  script:
  """
  awk 'FNR>1 || NR==1' *.csv > all.csv
  """
  // """
  // #!/usr/bin/env Rscript

  // location <- "~/local/R_libs/"; dir.create(location, recursive = TRUE  )
  // if(!require(tidyverse)){
  //   install.packages("tidyverse", lib = location, repos='https://cran.csiro.au')
  //   library(tidyverse, lib.loc = location)
  // }
  // """
}
process ggplot {
  echo true
  errorStrategy 'finish'
  label 'rscript'
  label 'stats'

  input:
    file csv from aggregatedStats

  output:
    file '*.pdf'

  script:
  """
  #!/usr/bin/env Rscript

  location <- "~/local/R_libs/"; dir.create(location, recursive = TRUE  )
  if(!require(tidyverse)){
    install.packages("tidyverse", lib = location, repos='https://cran.csiro.au')
    library(tidyverse, lib.loc = location)
  }
  #if(!require(hrbrthemes)){
  #  install.packages("hrbrthemes", lib = location, repos='https://cran.csiro.au')
  #  library(hrbrthemes, lib.loc = location)
  #}
  stats <- read_csv('${csv}') %>%
    mutate(adapters = case_when(adapters ~ "With adapters",  !adapters  ~ "No adapters"))

  #RUN-TIMES
  ggplot(stats)+ aes(tool, aligntime/1000) +
    geom_point(aes(colour = dataset))
  ggsave(file="runtimes.pdf", width=16, height=9);

  #ALIGNMENT RATES
  #Get those that report percentages
stat_perc <- stats %>%
  filter(perc)

stats_prop <- stats %>%
  filter(var %in% c("total_read_accuracy",
                    "perc_reads_incorrect",
                    "perc_reads_unaligned",
                    "perc_reads_ambiguous"),
         type == "unique")

ggplot(stats_prop, aes(tool, value_dbl)) +
  geom_bar(aes(fill = var), stat = "identity") +
  facet_wrap(adapters~dataset) +
  #scale_fill_viridis(discrete = TRUE, direction = -1) +
  labs(title = "Read alignment statistics ",
       subtitle = "uniquely aligned reads",
       x = "Tool",
       y = "Percentage",
       fill = "Alignment classification") #+  theme_ipsum_rc()
  ggsave(file="align-rates.pdf", width=16, height=9);
  """
}
