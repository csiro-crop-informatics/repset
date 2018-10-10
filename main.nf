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
    file(downloadedRef) from refs

  output:
    file ref
    // file('ucsc.hg19.fa') into reference

  script:
  ref='ucsc.hg19.fa'
  COMMON="tar xzvf ${downloadedRef}"
  if(params.debug) {
   """
    ${COMMON}
    cat chr1.fa > ${ref}
    """
  } else {
    """
    ${COMMON}
    cat \$(ls | grep -E 'chr([0-9]{1,2}|X|Y)\\.fa' | sort -V)  > ${ref}
    """
  }

}

process kangaIndex {
  label 'index'
  label 'biokanga'
  tag("${ref}")

  input:
    file(ref)

  output:
    file("*.sfx") into kangaRefs

  script:
    """
    biokanga index --threads ${task.cpus} -i ${ref} -o ${ref}.sfx --ref human
    """
}


// tools = ['biokanga','dart','hisat2']
// tools = ['biokanga','hisat2']

// process indexGenerator {
//   label 'index'
//   tag("${tool}")


//   input:
//     file ref
//     each tool from tools

//   output:
//     set val(meta), file("${ref}*") into indices



//   script:
//     label("${tool}")
//     meta = [tool: "${tool}", ref: "${ref}"]
//     //EITHER THIS:
//     switch(tool) {
//       case 'biokanga':
//         """
//         biokanga index --threads ${task.cpus} -i ${ref} -o ${ref}.sfx --ref ${ref}
//         """
//         break
//       case 'dart':
//         """
//         bwt_index ${ref} ${ref}
//         """
//         break
//       case 'hisat2':
//         """
//         hisat2-build ${ref} ${ref} -p ${task.cpus}
//         """
//         break
//       default:
//         break
//     }
//     //OR USE TEMPLATES:
//     // template "${tool}_index.sh"
// }

process dartIndex {
  label 'index'
  label 'dart'
  tag("${ref}")

  input:
    file(ref)

  output:
    set val("${ref}"), file("${ref}.*") into dartRefs

  script:
    """
    bwt_index ${ref} ${ref}
    """
}

process hisat2Index {
  label 'index'
  label 'hisat2'
  tag("${ref}")

  input:
    file(ref)

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
    set val(dataset), file("${ds}")  into datasetsChannel //datasetsForKanga, datasetsForHisat2
    // file('*') into extractedDatasets

  script:
    ds = dataset.replaceFirst("human","dataset")
    """
    mkdir -p ${ds}
    pbzip2 --decompress --stdout -p${task.cpus} ${dataset}.tar.bz2 | tar -x --directory ${ds}
    """
}

process addAdapters {

  tag("${dataset}")

  input:
    set val(dataset), file(dataDir) from datasetsChannel
  output:
    set val(dataset), file(dataDir)  into datasetsForKanga, datasetsForHisat2, datasetsForDart

  script:
  """
  add_adapter2fasta_V3.pl ${dataDir}/*.forward.fa ${dataDir}/*.reverse.fa ${dataDir}/forward.adapters.fa ${dataDir}/reverse.adapters.fa
  """
}



// process align {
//   label 'align'
//   // label("${idxmeta.tool}")
//   // tag("${idxmeta.tool}"+" VS "+"${idxmeta.ref}")
//   echo true

//   input:
//     set val(idxmeta), file("${ref}.*"), val(dataset), file(dataDir) from indices.combine(datasetsChannel) //cartesian product i.e. all input sets of reads vs all dbs - easy way of repeating ref for each dataset

//   script:
//   """
//   ls -l
//   """


// }

process kangaAlign {
  label 'align'
  label 'biokanga'
  tag("${dataset}"+" VS "+"${ref}")

  input:
    set file(ref), val(dataset), file(dataDir) from kangaRefs.combine(datasetsForKanga) //cartesian product i.e. all input sets of reads vs all dbs - easy way of repeating ref for each dataset
    //each pemode from [2,3]

  output:
    set val(meta), file(dataDir), file(samfile) into kangaAlignedDatasets

  script:
    meta = [tool: 'biokanga', id: dataset.replaceFirst("human_","")]
    samfile='aligned.sam'
    CMD = "biokanga align --sfx ${ref} \
      --mode 0 \
      --format 5 \
      --maxns 2 \
      --pemode 2 \
      --pairmaxlen 50000 \
      --in ${dataDir}/forward.adapters.fa \
      --pair ${dataDir}/reverse.adapters.fa  \
	    --out ${samfile} \
	    --threads ${task.cpus} "
    CMD += "--substitutions 5 \
      --minchimeric 50"
    CMD += params.debug ? ' -# 1000' : '' //every thousandth read/pair
    """
    ${CMD}
    """
}

process dartAlign {
  label 'align'
  label 'dart'
  tag("${dataset}"+" VS "+"${ref}")

  input:
    set val(ref), file("${ref}.*"), val(dataset), file(dataDir) from dartRefs.combine(datasetsForDart) //cartesian product i.e. all input sets of reads vs all dbs - easy way of repeating ref for each dataset

  output:
    set val(meta), file(dataDir), file(samfile) into dartAlignedDatasets

  script:
    meta = [tool: 'dart', id: dataset.replaceFirst("human_","")]
    samfile='aligned.sam'
    CMDSFX = params.debug ? '| head -10000' : ''
    """
    dart -i ${ref} -f ${dataDir}/forward.adapters.fa -f2 ${dataDir}/reverse.adapters.fa \
      --threads ${task.cpus} ${CMDSFX} > ${samfile}
    """
}

process hisat2Align {
  label 'align'
  label 'hisat2'
  tag("${dataset}"+" VS "+"${ref}")

  input:
    set val(ref), file("${ref}.*.ht2"), val(dataset), file(dataDir) from hisat2Refs.combine(datasetsForHisat2) //cartesian product i.e. all input sets of reads vs all dbs - easy way of repeating ref for each dataset

  output:
    set val(meta), file(dataDir), file(samfile) into hisat2AlignedDatasets

  //READ-level-optimized: default-1-20-0.5-25-5-20-1-0-3-0
// MODE=default (end-to-end, alt local)
// NUM_MISMATCH=1
// SEED_LENGTH=20
// SEED_INTERVAL=0.5
// SEED_EXTENSION=25
// RE_SEED=5
// PENALITY_NONCANONICAL=20
// MAX_MISMATCH_PENALITY=1
// MIN_MISMATCH_PENALITY=0
// MAX_SOFTCLIPPING_PENALITY=3
// MIN_SOFTCLIPPING_PENALITY=0
  script:
    meta = [tool: 'hisat2', id: dataset.replaceFirst("human_","")]
    samfile='aligned.sam'
    CMDSFX = params.debug ? '| head -10000' : ''
    """
    hisat2 -x ${ref} -1 ${dataDir}/forward.adapters.fa -2 ${dataDir}/reverse.adapters.fa \
      --time \
      --threads ${task.cpus} \
      --reorder \
      -f \
      ${CMDSFX} > ${samfile}
    """
}

process nameSortSAM {
  label 'samtools'
   tag("${meta}")
   input:
    set val(meta), file(dataDir), file(samfile) from hisat2AlignedDatasets.mix(kangaAlignedDatasets).mix(dartAlignedDatasets)

  output:
    set val(meta), file(dataDir), file(sortedsam) into sortedSAMs

  script:
  sortedsam = 'sorted.sam'
  """
  samtools sort -n --threads ${task.cpus} -o ${sortedsam} ${samfile}
  """
}

process fixSAM {
  label 'benchmark'
  tag("${meta}")

  input:
    set val(meta), file(dataDir), file(samfile) from sortedSAMs //kangaAlignedDatasets.first() //hisat2AlignedDatasets.first() //kangaAlignedDatasets.mix(hisat2AlignedDatasets)

  output:
    set val(meta), file(dataDir), file(fixedsam) into fixedSAMs

  script:
  """
  fix_sam.rb ${samfile} > fixedsam
  """

  //tu run when erubis not available, e.g.using modules on cluster and .rb scripts in bin/
  // """
  // gem install erubis --install-dir \$PWD
  // echo "source \'https://rubygems.org\'" > Gemfile
  // echo "gem \'erubis\'" >> Gemfile
  // bundle exec ruby ${baseDir}/bin/fix_sam.rb ${samfile} > fixedsam
  // """
}

process compareToTruth {
  label 'benchmark'
  tag("${meta}")

  input:
    set val(meta), file(dataDir), file(fixedsam) from fixedSAMs

  script:
  """
  compare2truth.rb ${dataDir}/*.cig ${fixedsam} > comp_res.txt
  """
}
