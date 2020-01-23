#!/usr/bin/env nextflow

import nextflow.util.Escape
import nextflow.container.SingularityCache


//Auxiliary container images
def containers = []
session.getConfig().process.each {k, v ->
  if((k.startsWith('withLabel:') || k.startsWith('withName:')) && v.containsKey('container') && !(k.endsWith('rrender'))) { //skipping rrender for CI
    // println "$k -> $v.container"
    containers << v.container
  }
}

//Mappers container images
Channel.from(params.mappersDefinitions)
.filter{ params.mappers == 'all' || it.tool.matches(params.mappers) }
.map {
  it.container
}
.set { mapperContainers }

SingularityCache scache = new SingularityCache() //to get NF-consitent image file names


process pull_container {
  tag { remote }
  maxForks 1
  storeDir "${params.singularitydir}"
  echo true

input:
  val(remote) from mapperContainers.mix(Channel.from(containers)).unique()

output:
  file(img)

script:
img = scache.simpleName(remote)
"""
singularity pull --name ${img} docker://${Escape.path(remote)}
"""
}

// String simpleName(String imageUrl) {
//     def p = imageUrl.indexOf('://')
//     def name = p != -1 ? imageUrl.substring(p+3) : imageUrl
//     String extension = '.img'
//     if( name.contains('.sif:') ) {
//         extension = '.sif'
//         name = name.replace('.sif:','-')
//     }
//     else if( name.endsWith('.sif') ) {
//         extension = '.sif'
//         name = name.substring(0,name.length()-4)
//     }
//     name = name.replace(':','-').replace('/','-')
//     return name + extension
// }