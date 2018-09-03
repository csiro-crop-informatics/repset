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
  storeDir "${workflow.workDir}/dataset/human/genome"

  input:
    file(ref) from refs

  output:
    file('ucsc.hg19.fa') into kangaRef
    file('ucsc.hg19.fa') into hisat2Ref

  script:
  """
  tar xzvf ${ref}
  cat \$(ls | grep -E 'chr([0-9]{1,2}|X|Y)\\.fa' | sort -V) > ucsc.hg19.fa
  """
}

process kangaIndex {
  label 'index'
  tag("${ref}")
  storeDir "${workflow.workDir}/tool_results/biokanga/index"
  //stageInMode 'link'
  // docker.runOptions = "--volume ${workflow.workDir}/dataset/human/:/project/itmatlab/aligner_benchmark/dataset/human/"
  // container = 'rsuchecki/biokanga_benchmark:0.1.2'
  // config.singularity.runOptions = "--writable --bind ${workflow.workDir}/dataset/human/:/project/itmatlab/aligner_benchmark/dataset/human/"
  // println(runOptions)

  input:
    file(ref) from kangaRef //not used from workdir

  output:
    file('*.sfx') into kangaRefs

   script:
  """
    mkdir -p ${workflow.workDir}/tool_results
    SINGULARITY_CACHEDIR=${workflow.workDir}/singularity
    singularity exec --writable --bind \
    ${workflow.workDir}/dataset/human/:/project/itmatlab/aligner_benchmark/dataset/human/ \
    --bind ${workflow.workDir}/tool_results/:/project/itmatlab/aligner_benchmark/tool_results/ \
    ${params.container} /bin/bash -c \
    "/bin/bash /project/itmatlab/aligner_benchmark/jobs/biokanga/biokanga-index.sh ${task.cpus}  \
    /project/itmatlab/aligner_benchmark/jobs/settings/dataset_human_hg19_t1r1.sh"
  """
}

// process hisat2Index {
//   label 'index'
//   tag("${ref}")
//   //storeDir "${workflow.workDir}/tool_results/hisat2/index"

//   input:
//     file(ref) from hisat2Ref //not used from workdir

//   output:
//     file('*') into hisat2Refs

//    script:
//   """
//     mkdir -p ${workflow.workDir}/tool_results
//     SINGULARITY_CACHEDIR=${workflow.workDir}/singularity
//     singularity exec --writable --bind \
//     ${workflow.workDir}/dataset/human/:/project/itmatlab/aligner_benchmark/dataset/human/ \
//     --bind ${workflow.workDir}/tool_results/:/project/itmatlab/aligner_benchmark/tool_results/ \
//     ${params.container} /bin/bash -c \
//     "/bin/bash /project/itmatlab/aligner_benchmark/jobs/hisat2/hisat2-index.sh ${task.cpus}  \
//     /project/itmatlab/aligner_benchmark/jobs/settings/dataset_human_hg19_t1r1.sh"
//   """
// }

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
  storeDir "${workflow.workDir}/dataset/human/${ds}"
  tag("${dataset}")

  input:
    set val(dataset), file("${dataset}.tar.bz2") from downloadedDatasets

  output:
    val(dataset) into datasetsForKanga
    file('*') //into extractedDatasets

  script:
    ds = dataset.replaceFirst("human","dataset")
    """
    tar xjvf ${dataset}.tar.bz2
    """
}

process kangaAlign {
  label 'align'
  storeDir "${workflow.workDir}/tool_results/biokanga/alignment/dataset_human_hg19_RefSeq_${ds}"
  tag("${dataset}"+" VS "+"${ref}")

  input:
    set val(dataset), file(ref) from datasetsForKanga.combine(kangaRefs) //cartesian product i.e. all input sets of reads vs all dbs

  output:
    set val(tool), val(ds) into kangaAlignedDatasets
    file('Aligned.out.sam')

  script:
    tool='biokanga'
    ds = dataset.replaceFirst("human_","")
  """
    mkdir -p ${workflow.workDir}/tool_results
    SINGULARITY_CACHEDIR=${workflow.workDir}/singularity
    singularity exec --writable \
    --bind ${workflow.workDir}/dataset/human/:/project/itmatlab/aligner_benchmark/dataset/human/ \
    --bind ${workflow.workDir}/tool_results/:/project/itmatlab/aligner_benchmark/tool_results/ \
    ${params.container} /bin/bash -c \
    "/bin/bash /project/itmatlab/aligner_benchmark/jobs/biokanga/biokanga-align-PE.sh ${task.cpus}   \
    /project/itmatlab/aligner_benchmark/jobs/settings/dataset_human_hg19_${ds}.sh -#100 "
  """
}


process benchmark {
  storeDir "${workflow.workDir}/statistics/human_${dataset}/${tool}"
  tag("${dataset} ${tool}")

  input:
    set val(tool), val(dataset) from kangaAlignedDatasets //.mix(hisat2AlignedDatasets)

  output:
    file('comp_res_multi_mappers.txt')
    file('comp_res.txt')
    file('*.sam')

  script:
  """
  mkdir -p ${workflow.workDir}/statistics
  SINGULARITY_CACHEDIR=${workflow.workDir}/singularity
  singularity exec --writable \
  --bind ${workflow.workDir}/dataset/human/:/project/itmatlab/aligner_benchmark/dataset/human/ \
  --bind ${workflow.workDir}/tool_results/:/project/itmatlab/aligner_benchmark/tool_results/ \
  --bind ${workflow.workDir}/statistics/:/project/itmatlab/aligner_benchmark/statistics/ \
  ${params.container} /bin/bash -c \
  "cd /project/itmatlab/aligner_benchmark && ruby master.rb -v ${dataset} ${dataset} /project/itmatlab/aligner_benchmark -a${tool}"
  """

}

// #!/bin/bash
// # indexing, generating alignments as PE only from the biokanga jobs directory and lastly processing for alignment statistics
// cd /project/itmatlab/aligner_benchmark/jobs/biokanga

// #generating index for human - only needs to be done once
// # ./biokanga-index.sh 48 /project/itmatlab/aligner_benchmark/jobs/settings/dataset_human_hg19_t1r1.sh

// #generating PE alignments for human t1r1
// ./biokanga-align-PE.sh 48 /project/itmatlab/aligner_benchmark/jobs/settings/dataset_human_hg19_t1r1.sh
// ./biokanga-align-PE.sh 48 /project/itmatlab/aligner_benchmark/jobs/settings/dataset_human_hg19_t2r1.sh
// ./biokanga-align-PE.sh 48 /project/itmatlab/aligner_benchmark/jobs/settings/dataset_human_hg19_t3r1.sh

// #now back to the benchmarking root directory and run the ruby scripts for generating stats on the alignments
// cd /project/itmatlab/aligner_benchmark
// ruby master.rb -v t1r1 t1r1 /project/itmatlab/aligner_benchmark -abiokanga
// ruby master.rb -v t2r1 t2r1 /project/itmatlab/aligner_benchmark -abiokanga
// ruby master.rb -v t3r1 t3r1 /project/itmatlab/aligner_benchmark -abiokanga



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