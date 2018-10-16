tools = ['biokanga','dart','hisat2']
datasets = ['human_t1r1','human_t1r2','human_t1r3']
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
    each tool from tools

  output:
    set val(meta), file("${ref}*") into indices

  script:
    meta = [tool: "${tool}", ref: "${ref}"]
    //EITHER THIS:
    switch(tool) {
      case 'biokanga':
        """
        biokanga index --threads ${task.cpus} -i ${ref} -o ${ref}.sfx --ref ${ref}
        """
        break
      case 'dart':
        """
        bwt_index ${ref} ${ref}
        """
        break
      case 'hisat2':
        """
        hisat2-build ${ref} ${ref} -p ${task.cpus}
        """
        break
      default:
        break
    }
    //OR USE TEMPLATES:
    // template "${tool}_index.sh"
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
    awk '\$2~/${params.debugChromosome}\$/' ${dataDir}/*.cig | cut -f1 > debug.ids
    fgrep --no-filename --no-group-separator -A1 -wf debug.ids ${dataDir}/*.forward.fa > r1
    fgrep --no-filename --no-group-separator -A1 -wf debug.ids ${dataDir}/*.reverse.fa > r2
    ln -sf "\$(readlink -f ${dataDir}/*.cig)" cig
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

process align {
  label 'align'
  // label("${idxmeta.tool}") // it is currently not possible to set dynamic process labels in NF, see https://github.com/nextflow-io/nextflow/issues/894
  container { this.config.process.get("withLabel:${idxmeta.tool}" as String).get("container") }

  tag("${idxmeta} << ${readsmeta}")
  // echo true

  input:
    set val(idxmeta), file("*"), val(readsmeta), file(r1), file(r2), file(cig) from indices.combine(datasetsWithAdapters.mix(preparedDatasets))
    // set val(idxmeta), file("*"), val(readsmeta), file(r1), file(r2), file(cig) from indices.combine(preparedDatasetsMultiChannel.poll().mix(datasetsWithAdaptersMultiChannel.poll()))

  output:
    set val(meta), file(sam), file(cig) into alignedDatasets

  script:
  meta = idxmeta.clone() + readsmeta.clone()
  //EITHER THIS:
    switch(idxmeta.tool) {
      case 'biokanga':
        """
        biokanga align --sfx ${idxmeta.ref}.sfx \
         --mode 0 \
         --format 5 \
         --maxns 2 \
         --pemode 2 \
         --pairmaxlen 50000 \
         --in ${r1} \
         --pair ${r2}  \
         --out sam \
         --threads ${task.cpus} \
         --substitutions 5 \
         --minchimeric 50
        """
        break
      case 'dart':
        """
        dart -i ${idxmeta.ref} -f ${r1} -f2 ${r2} -t ${task.cpus} > sam
        """
        break
      case 'hisat2':
        """
        hisat2 -x ${idxmeta.ref} -1 ${r1} -2 ${r2} \
        --time \
        --threads ${task.cpus} \
        --reorder \
        -f \
        > sam
        """
        break
      default: //If tool not specified above expecting a template in templates/
        template "${tool}_align.sh"
        break
    }
    //OR USE TEMPLATES:
    // template "${tool}_align.sh"
}

process nameSortSAM {
   tag("${meta}")
   input:
    set val(meta), file(sam), file(cig) from alignedDatasets //hisat2AlignedDatasets.mix(kangaAlignedDatasets, dartAlignedDatasets)

  output:
    set val(meta), file(sortedsam), file(cig) into sortedSAMs

  script:
    """
    grep -v '^@' ${sam} | sort -t'.' -k2,2n --parallel ${task.cpus} > sortedsam
    """
}

process fixSAM {
  label 'benchmark'
  tag("${meta}")

  input:
    set val(meta), file(sortedsam), file(cig) from sortedSAMs

  output:
    set val(meta), file(fixedsam), file(cig) into fixedSAMs

  script:
  """
  fix_sam.rb ${sortedsam} > fixedsam
  """
}

// actions = ['compare2truth', 'compare2truth_multi_mappers'] //
actions = ['compare2truth']
process compareToTruth {
  label 'benchmark'
  label 'stats'
  tag("${meta}")

  input:
    set val(meta), file(fixedsam), file(cig) from fixedSAMs
    each action from actions

  output:
    file '*' into stats

  script:
  outname = meta.dataset+""+(meta.adapters ? "_adapters" : "")+"_"+meta.tool+"."+action
  """
  ${action}.rb ${cig} ${fixedsam} > ${outname}
  """
}

process aggregateStats {
  label 'benchmark'
  label 'stats'

  input:
    file '*' from stats.collect()

  output:
    file '*_combined.txt'

  script:
  """
  for action in ${actions.join(" ")}; do
    ls *.\${action} | sort -V | xargs read_stats.rb | sed "s/\${action}//g" > \${action}_combined.txt
  done
  """
}