tools = ['biokanga','dart','hisat2']
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
    cat ${params.debugChromosome}.fa > ${ref}
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
    set val("${ref}"), file("${ref}*") into dartRefs

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
  label 'download'

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
    set val(dataset), file("${dataset}") into extractedDatasets

  script:
    """
    mkdir -p ${dataset}
    pbzip2 --decompress --stdout -p${task.cpus} ${dataset}.tar.bz2 | tar -x --directory ${dataset}
    """
}

process prepareDatasets {
  tag("${dataset}")

  input:
    set val(dataset), file(dataDir) from extractedDatasets

  output:
    set val(meta), file(r1), file(r2), file(cig) into preparedDatasets, prepareDatasetsForAdapters

  script:
  meta = [dataset: dataset, adapters: false]
  if(params.debug) {
    """
    awk '\$2~/${params.debugChromosome}\$/' ${dataDir}/*.cig | sort -k1,1V | head -n ${params.debugReads} > cig
    fgrep --no-group-separator -A1 -f <(cut -f1 cig) ${dataDir}/*forward.fa > r1
    fgrep --no-group-separator -A1 -f <(cut -f1 cig) ${dataDir}/*reverse.fa > r2
    """
  } else {
    """
    cp --no-dereference ${dataDir}/*forward.fa r1
    cp --no-dereference ${dataDir}/*reverse.fa r2
    cp --no-dereference ${dataDir}/*.cig cig
    """
  }
}

//Each tool to get each dataset - multiplicate the channel and put these on Queue so that each tool can take one
preparedDatasetsMultiChannel = preparedDatasets.into(tools.size()) as Queue

process addAdapters {
  tag("${meta.dataset}")

  input:
    set val(inmeta), file(r1), file(r2), file(cig) from prepareDatasetsForAdapters //(preparedDatasetsMultiChannel.poll())
  output:
    set val(meta), file(a1), file(a2), file(cig)  into datasetsWithAdapters //datasetsForKanga, datasetsForHisat2, datasetsForDart

  script:
  meta = inmeta.clone()
  meta.adapters = true
    """
    add_adapter2fasta_V3.pl ${r1} ${r2} a1 a2
    """
}

//Each tool to get each dataset - multiplicate the channel and put these on Queue so that each tool can take one
datasetsWithAdaptersMultiChannel = datasetsWithAdapters.into(tools.size) as Queue


// // process align {
// //   label 'align'
// //   // label("${idxmeta.tool}")
// //   // tag("${idxmeta.tool}"+" VS "+"${idxmeta.ref}")
// //   echo true

// //   input:
// //     set val(idxmeta), file("${ref}.*"), val(dataset), file(dataDir) from indices.combine(datasetsChannel) //cartesian product i.e. all input sets of reads vs all dbs - easy way of repeating ref for each dataset

// //   script:
// //   """
// //   ls -l
// //   """
// // }

process kangaAlign {
  label 'align'
  label 'biokanga'
  tag("${inmeta}"+" VS "+"${ref}")

  input:
    set file(ref), val(inmeta), file(r1), file(r2), file(cig) from kangaRefs.combine(preparedDatasetsMultiChannel.poll().mix(datasetsWithAdaptersMultiChannel.poll())) //cartesian product i.e. all input sets of reads vs all dbs - easy way of repeating ref for each dataset

    //each pemode from [2,3]

  output:
    set val(meta), file(sam), file(cig) into kangaAlignedDatasets

  script:
    meta = inmeta.clone() + [tool: 'biokanga']
    CMD = "biokanga align --sfx ${ref} \
      --mode 0 \
      --format 5 \
      --maxns 2 \
      --pemode 2 \
      --pairmaxlen 50000 \
      --in ${r1} \
      --pair ${r2}  \
	    --out sam \
	    --threads ${task.cpus} "
    CMD += "--substitutions 5 \
      --minchimeric 50"
    """
    ${CMD}
    """
}

process dartAlign {
  label 'align'
  label 'dart'
  tag("${inmeta}"+" VS "+"${ref}")

  input:
    set val(ref), file("*"), val(inmeta), file(r1), file(r2), file(cig) from dartRefs.combine(preparedDatasetsMultiChannel.poll().mix(datasetsWithAdaptersMultiChannel.poll())) //cartesian product i.e. all input sets of reads vs all dbs - easy way of repeating ref for each dataset

  output:
    set val(meta), file(sam), file(cig) into dartAlignedDatasets

  script:
    meta = inmeta.clone() + [tool: 'dart']
    """
    dart -i ${ref} -f ${r1} -f2 ${r2} -t ${task.cpus} > sam
    """
}

process hisat2Align {
  label 'align'
  label 'hisat2'
  tag("${inmeta}"+" VS "+"${ref}")

  input:
    set val(ref), file("${ref}.*.ht2"), val(inmeta), file(r1), file(r2), file(cig) from hisat2Refs.combine(preparedDatasetsMultiChannel.poll().mix(datasetsWithAdaptersMultiChannel.poll())) //cartesian product i.e. all input sets of reads vs all dbs - easy way of repeating ref for each dataset

  output:
    set val(meta), file(sam), file(cig) into hisat2AlignedDatasets

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
    meta = inmeta.clone() + [tool: 'hisat2']
    """
    hisat2 -x ${ref} -1 ${r1} -2 ${r2} \
      --time \
      --threads ${task.cpus} \
      --reorder \
      -f \
      > sam
    """
}

// process nameSortSAM {
//   label 'samtools'
//    tag("${meta}")
//    input:
//     set val(meta), file(dataDir), file(samfile) from hisat2AlignedDatasets.mix(kangaAlignedDatasets).mix(dartAlignedDatasets)

//   output:
//     set val(meta), file(dataDir), file(sortedsam) into sortedSAMs

//   script:
//   sortedsam = 'sorted.sam'
//   """
//   samtools sort -n --threads ${task.cpus} -o ${sortedsam} ${samfile}
//   """
// }

// process fixSAM {
//   label 'benchmark'
//   tag("${meta}")

//   input:
//     set val(meta), file(dataDir), file(samfile) from sortedSAMs //kangaAlignedDatasets.first() //hisat2AlignedDatasets.first() //kangaAlignedDatasets.mix(hisat2AlignedDatasets)

//   output:
//     set val(meta), file(dataDir), file(fixedsam) into fixedSAMs

//   script:
//   """
//   fix_sam.rb ${samfile} > fixedsam
//   """

//   //tu run when erubis not available, e.g.using modules on cluster and .rb scripts in bin/
//   // """
//   // gem install erubis --install-dir \$PWD
//   // echo "source \'https://rubygems.org\'" > Gemfile
//   // echo "gem \'erubis\'" >> Gemfile
//   // bundle exec ruby ${baseDir}/bin/fix_sam.rb ${samfile} > fixedsam
//   // """
// }

// // actions = ['compare2truth', 'compare2truth_multi_mappers'] //
// actions = ['compare2truth']
// process compareToTruth {
//   label 'benchmark'
//   label 'stats'
//   tag("${meta}")

//   input:
//     set val(meta), file(dataDir), file(fixedsam) from fixedSAMs
//     each action from actions

//   output:
//     file '*' into stats

//   script:
//   """
//   ${action}.rb ${dataDir}/*.cig ${fixedsam} > "${meta.tool}_${meta.dataset}.${action}"
//   """
// }

// process aggregateStats {
//   label 'benchmark'
//   label 'stats'

//   input:
//     file '*' from stats.collect()

//   output:
//     file '*_combined.txt'

//   script:
//   """
//   for action in ${actions.join(" ")}; do
//     ls *.\${action} | sort | xargs read_stats.rb | sed "s/\${action}//g" > \${action}_combined.txt
//   done
//   """
// }