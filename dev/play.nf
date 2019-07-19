import static groovy.json.JsonGenerator.*
jsonGenerator = new groovy.json.JsonGenerator.Options()
                .addConverter(java.nio.file.Path) { java.nio.file.Path p, String key -> p.toUriString() }
                .build()


// x = ['containers':['bbmap':['38.44':'rsuchecki/bbmap:38.44_fae5e1e07240e69896dbf7095872fb6fea43d045', '38.49':'rsuchecki/bbmap:38.49_9e975d9bc6a657bc4306f4475be393b9fbe8e3fb'], 'minimap2':['2.17':'rsuchecki/minimap2:2.17_1d3f326820696496f025a95632979cd4ea4140cb']]]
// Channel.from(x.containers).flatMap().view()

// Channel.from(params.containers)
// .flatMap()
// .view {
//   println it
//   groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))
// }

// println params.mappers
// println groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(params.mappers))

def allModes = 'dna2dna|rna2rna|rna2dna'
mappers = Channel.from(params.mappers)
// .view()
.filter { record ->
   ['tool','version','container','index'].every { it in record.keySet() } \
   && allModes.split('\\|').any { it in record.keySet() }
}

def String bindTemplate(template, binding) {
  def engine = new groovy.text.SimpleTemplateEngine()
  engine.createTemplate(template).make(binding).toString()
}

process index {
  memory '100.MB'
  echo true

  input:
    val mapper from mappers

  output:
    set val(mapper), val(ref) into indices

  exec:
  // script:
  ref = 'reference/path.fa'
  binding = [ref: ref, task: task]
  println bindTemplate(mapper.index, binding) // "${bindTemplate(mapper.index, binding)}"
  // bindTemplate(mapper.index, binding)
}

/*
Read, sanitize and validate alignment/mapping param sets
*/
def validationMaps = [dna2dna: [:], rna2dna: [:], rna2rna: [:]]
def versionValidationMaps = [dna2dna: [:], rna2dna: [:], rna2rna: [:]]
// println groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(params.mapperParams))
mapperParams = Channel.from params.mapperParams
.each { //PUT IN DEFAULT VALUES
  if(!it.containsKey('label')) it.put('label', 'default');
  if(!it.containsKey('version')) it.put('version', 'ALL_AVAILABLE');
  if(!it.containsKey('mode')) it.put('mode', allModes);
}.each { rec -> //VALIDATE
  ['dna2dna', 'rna2dna', 'rna2rna']. each { mode ->
    validationMap = validationMaps."${mode}"
    if(mode.matches(rec.mode)) {
      key = [rec.tool, rec.version, rec.label].join("_")
      stored = validationMap.putIfAbsent(key, rec)
      if(stored != null) {
        log.error """Validation error: non-unique label for aligner params set for ${mode}
        Previously encountered: ${stored}
        Offending record: ${rec}
        If tool/version/mode overlap then the labels must be unique"""
        System.exit 1
      }
      versionValidationMap = versionValidationMaps."${mode}"
      noVersionKey = [rec.tool, rec.label].join(" label:")
      if(versionValidationMap.containsKey(noVersionKey)) {
        storedVersions = versionValidationMap.get(noVersionKey)
        storedVersions.add(rec.version)
        versionValidationMap.put(noVersionKey,storedVersions)
      } else {
        versionValidationMap.put(noVersionKey, [rec.version])
      }
    }
  }
}
versionValidationMaps.each { mode -> //VALIDATE SOME MORE
  mode.value.each { tool_label, versions ->
    if(versions.size() > 1 && 'ALL_AVAILABLE' in versions) {
      log.error """Validation error: ambiguous assignment of alignment parameter set labels.
        Parameter set assigned to 'ALL_AVAILABLE' as well as a specific version of a tool.
        ${tool_label} : ${versions}
        If tool/version/mode overlap then the labels must be unique"""
        System.exit 1
    }
  }
}
// println groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(versionValidationMaps))
// mapperParams.view()

process mapReads {
  input:

  exec:

}