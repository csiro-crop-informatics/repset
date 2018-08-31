#!/usr/bin/env nextflow

datasets = ['human_t1r1','human_t1r2','human_t1r3']

process downloadReference {
  storeDir '${workflow.workDir}/dataset/human/genome'

  output:
    file('ucsc.hg19.fa') into kangaRef, hisat2Ref

  shell:
    """
      wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
      tar xzvf chromFa.tar.gz
      cat $(ls | grep -E 'chr([0-9]{1,2}|X|Y)\.fa') > ucsc.hg19.fa
    """
}

process downloadDatasets {
  storeDir '${workflow.workDir}/dataset/human/${dataset}'

  input:
    each dataset from datasets

  output:
    val('${dataset}') into downloadedDatasets


  script:
    """
    wget http://bp1.s3.amazonaws.com/${dataset}.tar.bz2
    tar xjvf ${dataset}.tar.bz2
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


process kangaIndex {
  echo true
  // storeDir '${workflow.workDir}/dataset/human/${dataset}'
  //stageInMode 'link'
  runOptions = '--volume ${workflow.workDir}/dataset/human/:/project/itmatlab/aligner_benchmark/dataset/human/'
  container = 'rsuchecki/biokanga_benchmark:0.1.1'

  input:
    file('ucsc.hg19.fa') from kangaRef //not used from workdir

   script:
   // # indexing, generating alignments as PE only from the biokanga jobs directory and lastly processing for alignment statistics
  // cd /project/itmatlab/aligner_benchmark/jobs/biokanga
    """
    /project/itmatlab/aligner_benchmark/jobs/biokanga/biokanga-index.sh ${process.cpus} /project/itmatlab/aligner_benchmark/jobs/settings/dataset_human_hg19_t1r1.sh
    """
}
// process execute {
//   echo true
//   stageInMode 'link'
//   runOptions = '--volume ${workflow.workDir}/dataset/human/:/project/itmatlab/aligner_benchmark/dataset/human/'
//   container = 'rsuchecki/biokanga_benchmark:0.1.1'

//   input:
//     set filefile('hg19.chrom.sizes') from ref
//     val() from downloadDatasets
//     // file('ucsc.hg19.fa')

//    script:
//     """
//     ls -l
//     head *
//     """
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