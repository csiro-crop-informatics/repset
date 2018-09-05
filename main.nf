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
  //storeDir "${workflow.workDir}/dataset/human/genome"

  input:
    file(ref) from refs

  output:
    file('ucsc.hg19.fa') into kangaRef
    file('ucsc.hg19.fa') into hisat2Ref

  script:
  """
  tar xzvf ${ref}
  cat \$(ls | grep -E 'chr([0-9]{1,2}|X|Y)\\.fa' | sort -V) | head -10000 > ucsc.hg19.fa
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
  tag("${dataset}")

  input:
    set val(dataset), file("${dataset}.tar.bz2") from downloadedDatasets

  output:
    set val(dataset), file("${dataset}")  into datasetsForKanga
    // file('*') into extractedDatasets

  script:
    //ds = dataset.replaceFirst("human","dataset")
    """
    mkdir -p ${dataset}
    pbzip2 --decompress --stdout -p${task.cpus} ${dataset}.tar.bz2 | tar -x --directory ${dataset}
    """
}

process kangaAlign {
  label 'align'
  tag("${dataset}"+" VS "+"${ref}")
  module = 'biokanga/4.3.9'

  input:
    set file(ref), val(dataset), file(dataDir) from kangaRefs.combine(datasetsForKanga) //cartesian product i.e. all input sets of reads vs all dbs - easy way of repeating ref for each dataset

  output:
    set val(meta), file(dataDir), file(outfile) into kangaAlignedDatasets

  script:
    meta = [tool: 'biokanga', id: dataset.replaceFirst("human_","")]
    outfile='Aligned.out.sam'
    CMD = "biokanga align --sfx ${ref} \
      --log kanga_align.log \
      --mode 0 \
      --format 5 \
      --maxns 2 \
      --pemode 2 \
      --pairmaxlen 50000 \
      --in ${dataDir}/*.forward.fa \
      --pair ${dataDir}/*.reverse.fa  \
	    --out ${outfile} \
	    --threads ${task.cpus} "
    CMD += "--substitutions 5 \
      --minchimeric 50 \
      --samplenthrawread 1000"
    """
    ${CMD}
    """
}
// process kangaAlign {
//   label 'align'
//   storeDir "${workflow.workDir}/tool_results/biokanga/alignment/dataset_human_hg19_RefSeq_${ds}"
//   tag("${dataset}"+" VS "+"${ref}")

//   input:
//     set val(dataset), file(ref) from datasetsForKanga.combine(kangaRefs) //cartesian product i.e. all input sets of reads vs all dbs

//   output:
//     set val(tool), val(ds) into kangaAlignedDatasets
//     file('Aligned.out.sam')

//   script:
//     tool='biokanga'
//     ds = dataset.replaceFirst("human_","")
//   """
//     mkdir -p ${workflow.workDir}/tool_results
//     SINGULARITY_CACHEDIR=${workflow.workDir}/singularity
//     singularity exec --writable \
//     --bind ${workflow.workDir}/dataset/human/:/project/itmatlab/aligner_benchmark/dataset/human/ \
//     --bind ${workflow.workDir}/tool_results/:/project/itmatlab/aligner_benchmark/tool_results/ \
//     ${params.container} /bin/bash -c \
//     "/bin/bash /project/itmatlab/aligner_benchmark/jobs/biokanga/biokanga-align-PE.sh ${task.cpus}   \
//     /project/itmatlab/aligner_benchmark/jobs/settings/dataset_human_hg19_${ds}.sh -#100 "
//   """
// }


process benchmark {
  //storeDir "${workflow.workDir}/statistics/human_${dataset}/${tool}"
  echo true
  tag("${meta}")

  input:
    set val(meta), file(dataDir), file(sam) from kangaAlignedDatasets.first() //.mix(hisat2AlignedDatasets)

  // output:
  //   set val(meta), file("dataset_${dataset}") into benchmarkedStats

//work/tool_results/biokanga/alignment/dataset_human_hg19_RefSeq_t1r1/Aligned.out.sam
  script:
  """
  mkdir -p statistics biokanga/alignment/dataset_human_hg19_RefSeq_${meta.id}
  #place SAM where aligner_benchmark expects it
  cp --preserve=links ${sam} biokanga/alignment/dataset_human_hg19_RefSeq_${meta.id}/
  SINGULARITY_CACHEDIR=${workflow.workDir}/singularity
  echo singularity exec --writable \
    --bind \${PWD}/${dataDir}:/project/itmatlab/aligner_benchmark/dataset/human/ \
    --bind \${PWD}:/project/itmatlab/aligner_benchmark/tool_results/ \
    --bind \${PWD}/statistics:/project/itmatlab/aligner_benchmark/statistics/ \
    ${params.container} /bin/bash -c \
    "cd /project/itmatlab/aligner_benchmark && ruby master.rb -v ${meta.id} ${meta.id} /project/itmatlab/aligner_benchmark -a${meta.tool}"
  """

}

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