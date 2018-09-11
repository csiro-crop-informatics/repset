datasets = ['human_t1r1','human_t1r2','human_t1r3']
url = 'http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz'

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

process convertReference {
  input:
    file(ref) from refs

  output:
    file('ucsc.hg19.fa') into kangaRef
    file('ucsc.hg19.fa') into hisat2Ref

  script:
  """
  tar xzvf ${ref}
  cat \$(ls | grep -E 'chr([0-9]{1,2}|X|Y)\\.fa' | sort -V)  > ucsc.hg19.fa
  """
}

process kangaIndex {
  label 'index'
  tag("${ref}")
  module = 'biokanga/4.3.9'

  input:
    file(ref) from kangaRef

  output:
    file("*.sfx") into kangaRefs

  script:
    """
    biokanga index --threads ${task.cpus} -i ${ref} -o ${ref}.sfx --ref human
    """
}

process hisat2Index {
  label 'index'
  tag("${ref}")
  module = 'hisat/2.0.5'

  input:
    file(ref) from hisat2Ref

  output:
    set val("${ref}"), file("${ref}.*.ht2") into hisat2Refs

  script:
    """
    hisat2-build ${ref} ${ref} -p ${task.cpus}
    """
}



process downloadDatasets {
  //storeDir "${workflow.workDir}/downloaded/${dataset}" and put the downloaded datasets there and prevent generating cost to dataset creators through repeated downloads
  tag("${dataset}")

  input:
    each dataset from datasets

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
    set val(dataset), file("${ds}")  into datasetsForKanga, datasetsForHisat2
    // file('*') into extractedDatasets

  script:
    ds = dataset.replaceFirst("human","dataset")
    """
    mkdir -p ${ds}
    pbzip2 --decompress --stdout -p${task.cpus} ${dataset}.tar.bz2 | tar -x --directory ${ds}
    """
}

process kangaAlign {
  label 'align'
  tag("${dataset}"+" VS "+"${ref}")
  module = 'biokanga/4.3.9'

  input:
    set file(ref), val(dataset), file(dataDir) from kangaRefs.combine(datasetsForKanga) //cartesian product i.e. all input sets of reads vs all dbs - easy way of repeating ref for each dataset
    //each pemode from [2,3]

  output:
    set val(meta), file(dataDir), file(alignDir) into kangaAlignedDatasets

  script:
    meta = [tool: 'biokanga', id: dataset.replaceFirst("human_","")]
    outfile='Aligned.out.sam'
    alignDir='aligndir'
    alignPath="${alignDir}/${meta.tool}/alignment/dataset_human_hg19_RefSeq_${meta.id}"
    CMD = "biokanga align --sfx ${ref} \
      --mode 0 \
      --format 5 \
      --maxns 2 \
      --pemode 2 \
      --pairmaxlen 50000 \
      --in ${dataDir}/*.forward.fa \
      --pair ${dataDir}/*.reverse.fa  \
	    --out ${alignPath}/${outfile} \
      --log ${alignPath}/${meta.tool}.log \
	    --threads ${task.cpus} "
    CMD += "--substitutions 5 \
      --minchimeric 50"
    """

    mkdir -p ${alignPath}
    ${CMD}
    """
}


process hisat2Align {
  label 'align'
  tag("${dataset}"+" VS "+"${ref}")
  module = 'hisat/2.0.5'

  input:
    set val(ref), file("${ref}.*.ht2"), val(dataset), file(dataDir) from hisat2Refs.combine(datasetsForHisat2) //cartesian product i.e. all input sets of reads vs all dbs - easy way of repeating ref for each dataset

  output:
    set val(meta), file(dataDir), file(alignDir) into hisat2AlignedDatasets

  script:
    meta = [tool: 'hisat2', id: dataset.replaceFirst("human_","")]
    //guess what, it appears these file names need to be hardcoded for the framework to work
    outfile='output.sam'
    alignDir='aligndir'
    alignPath="${alignDir}/${meta.tool}/alignment/dataset_human_hg19_RefSeq_${meta.id}"
    """
    mkdir -p ${alignPath}
    hisat2 -x ${ref} -1 ${dataDir}/*.forward.fa -2 ${dataDir}/*.reverse.fa \
      --time \
      --threads ${task.cpus} \
      --reorder \
      -f \
      -S ${alignPath}/${outfile} \
      2> ${alignPath}/${meta.tool}.log
    """
}

process benchmark {
  //echo true
  // publishDir 'results/'
  tag("${meta}")
  module = 'singularity/2.5.0'
  beforeScript = "SINGULARITY_CACHEDIR=${workflow.workDir}/singularity"

  input:
    set val(meta), file(dataDir), file(alignDir) from kangaAlignedDatasets.mix(hisat2AlignedDatasets) //kangaAlignedDatasets.first() //hisat2AlignedDatasets.first() //kangaAlignedDatasets.mix(hisat2AlignedDatasets)

  output:
    file("statistics") into benchmarkedStats
    // set val(meta), file("statistics/human_${meta.id}/${meta.tool}/*.txt") into benchmarkedStats

//work/tool_results/biokanga/alignment/dataset_human_hg19_RefSeq_t1r1/Aligned.out.sam
  script:
  """
  mkdir -p statistics
  singularity exec --writable \
    --bind ${dataDir}:/project/itmatlab/aligner_benchmark/dataset/human/${dataDir} \
    --bind ${alignDir}:/project/itmatlab/aligner_benchmark/tool_results/ \
    --bind \${PWD}/statistics:/project/itmatlab/aligner_benchmark/statistics/ \
    docker://rsuchecki/biokanga_benchmark:0.4 /bin/bash -c \
    "cd /project/itmatlab/aligner_benchmark && ruby master.rb -v ${meta.id} ${meta.id} /project/itmatlab/aligner_benchmark -a${meta.tool}"
  """
}

process collectStats {
  echo true
  module = 'singularity/2.5.0'
  beforeScript = "SINGULARITY_CACHEDIR=${workflow.workDir}/singularity"

  input:
    file("statistics*") from benchmarkedStats.collect()

  script:
  """
  mkdir collected
  for d in statistics*; do
    rsync -a --exclude '*.sam' \${d}/ statistics/
  done
  #tree -hl statistics
  singularity exec --writable \
    --bind \${PWD}/statistics:/project/itmatlab/aligner_benchmark/statistics/ \
    docker://rsuchecki/biokanga_benchmark:0.4 /bin/bash -c \
    "cd /project/itmatlab/aligner_benchmark \
      && find . -name comp_res.txt |sort | xargs ruby ./read_stats.rb >> statistics/default_summary.txt"
  """

}


// process benchmark2 {
//   echo true
//   tag("${meta}")
//   module = 'singularity/2.5.0'
//   beforeScript = "mkdir -p statistics biokanga/alignment/dataset_human_hg19_RefSeq_${meta.id} \
//                   cp --preserve=links ${sam} biokanga/alignment/dataset_human_hg19_RefSeq_${meta.id}/"
//   singularity {
//       enabled = true
//       // autoMounts = true
//       cacheDir = "${workflow.workDir}/singularity"
//       runOptions = "--writable \
//                     --bind ${dataDir}:/project/itmatlab/aligner_benchmark/dataset/human/${dataDir} \
//                     --bind \${PWD}:/project/itmatlab/aligner_benchmark/tool_results/ \
//                     --bind \${PWD}/statistics:/project/itmatlab/aligner_benchmark/statistics/"
//       container = 'rsuchecki/biokanga_benchmark:0.3.4'
//   }

//   input:
//     set val(meta), file(dataDir), file(sam) from kangaAlignedDatasets //.mix(hisat2AlignedDatasets)

//   // output:
//   //   set val(meta), file("dataset_${dataset}") into benchmarkedStats

// //work/tool_results/biokanga/alignment/dataset_human_hg19_RefSeq_t1r1/Aligned.out.sam
//   script:
//   """
//   cd /project/itmatlab/aligner_benchmark && ruby master.rb -v ${meta.id} ${meta.id} /project/itmatlab/aligner_benchmark -a${meta.tool}
//   """
// }

// process plot {

//   input:
//     set val(tool), val(dataset), file(statsDir) from benchmarkedStats

//   script:
//   """
//   ls -l
//   """
// }

// process benchmark {
//   storeDir "${workflow.workDir}/statistics/human_${dataset}/${tool}"
//   tag("${dataset} ${tool}")

//   input:
//     set val(tool), val(dataset) from kangaAlignedDatasets //.mix(hisat2AlignedDatasets)

//   output:
//     file('comp_res_multi_mappers.txt')
//     file('comp_res.txt')
//     file('*.sam')

//   script:
//   """
//   mkdir -p ${workflow.workDir}/statistics
//   SINGULARITY_CACHEDIR=${workflow.workDir}/singularity
//   singularity exec --writable \
//   --bind ${workflow.workDir}/dataset/human/:/project/itmatlab/aligner_benchmark/dataset/human/ \
//   --bind ${workflow.workDir}/tool_results/:/project/itmatlab/aligner_benchmark/tool_results/ \
//   --bind ${workflow.workDir}/statistics/:/project/itmatlab/aligner_benchmark/statistics/ \
//   ${params.container} /bin/bash -c \
//   "cd /project/itmatlab/aligner_benchmark && ruby master.rb -v ${dataset} ${dataset} /project/itmatlab/aligner_benchmark -a${tool}"
//   """

// }


/*
This script:
Pulls test datasets if not locally present
  TODO: define input.conf and/or input.json
    BEERS DATASETS
    hg19 download (genome, annotations)
      ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz
      ftp://ftp.ensembl.org/pub/grch37/update/fasta/homo_sapiens/dna/
      http://itmat.rum.s3.amazonaws.com/indexes/hg19_genome_one-line-seqs.fa.gz


      ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
      or ftp://ftp.ensembl.org/pub/grch37/update/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gtf.gz


Executes experimental alignemts
Generates figures and/or tables for the manuscript.
Due to... unusal way the experimental pipeline is implemented, what we need/have to do here is to ensure portability.
Multiple hardcoded paths need to be inserted/replaced.
  TODO: turn appropriate script files into templates to be used at pipeline runtime
Required software to be made available via container(s).
Custom scripts to be included in same repo - under bin/ or inline
 */

 /*
  Processes required

  With and without retained adapters?
  biokanga_align
    defaults
    optimised
  hista2_align
    defaults
    optimised

  */