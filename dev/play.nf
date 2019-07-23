import static groovy.json.JsonGenerator.*
jsonGenerator = new groovy.json.JsonGenerator.Options()
                .addConverter(java.nio.file.Path) { java.nio.file.Path p, String key -> p.toUriString() }
                .build()

//Input validation specified elswhere
def validators = new GroovyShell().parse(new File("$baseDir/Validators.groovy"))

//Read, parse, validate and sanitize alignment/mapping tools config
def allRequired = ['tool','version','container','index'] //Fields required for each tool in config
def allModes = 'dna2dna|rna2rna|rna2dna' //At leas one mode has to be defined as supported by each tool
def allVersions = validators.validateMappersDefinitions(params.mappers, allRequired, allModes)

//Read, sanitize and validate alignment/mapping param sets
validators.validateMapperParamsDefinitions(params.mapperParams, allVersions)

//Validated now, so gobble up mappers and their params definitions
mappersChannel = Channel.from(params.mappers)
.view {groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))}
mappersParamsChannel = Channel.from(params.mapperParams)


// println groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(allVersions))
// println groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(params.mapperParams))


// def String templateToScript(template, binding) {
//   def engine = new groovy.text.SimpleTemplateEngine()
//   engine.createTemplate(template).make(binding).toString()
// }

// process index {
//   memory '100.MB'
//   tag { "${mapper.container}" }
//   container { "${mapper.container}" }
//   //echo true

//   input:
//     val mapper from mappers

//   output:
//     val mapper into indices

//   exec:
//   // println (groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(mapper)))
//   def ref = 'reference/path.fa'
//   def binding = [ref: ref, task: task.clone()]
//   // println bindTemplate(mapper.index, binding) // "${bindTemplate(mapper.index, binding)}"
//   script:
//   templateToScript(mapper.index, binding)
// }

/*
Read, sanitize and validate alignment/mapping param sets
// */
// validateMapperParamsDefinitions(params.mapperParams, allVersions)

// println groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(params.mapperParams))







// mapperParams.view { groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it)) }
// // println groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(versionValidationMaps))
// // println groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(versionValidationMaps))
// mapperParams.view { groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it)) }
// // mapperParams.view()

// readsChn = Channel.from(['some_reads', 'some_more_reads'])

// // process mapReads {
// //   tag { "${mapper.tool} ${reads} ${par}" }
// //   container { "${mapper.container}" }
// //   // echo true
// //   input:
// //      set val(mapper), val(reads), val(par) from  indices.combine(readsChn).combine(mapperParams)
// //     //  set val(mapper), val(ref) from  indices

// //   when:
// //     mapper.tool == par.tool && mapper.version.matches()
// //   // script:
// //   exec:
// //     def binding = [task: task.clone(), reads: reads]
// //     // println "${mapper} ${ref}"
// //     // println bindTemplate(mapper.index, binding) // "${bindTemplate(mapper.index, binding)}"
// //     //templateToScript(mapper.dna2dna, binding)
// // }

// // indices.view()






// def addToListInMap (map, key, value) {
//   if(map.containsKey(key)) {
//     stored = map.get(key)
//     stored.add(value)
//     map.put(key,stored)
//   } else if(value instanceof Collection) {
//     map.put(key, value)
//   } else {
//     map.put(key, [value])
//   }
// }

// def validateMappersDefinitions (mappers, allRequired, allModes) {
//   def allVersions = [:] //Keep track of tool versions declared in config
//   mappers.each { rec ->
//     addToListInMap(allVersions, rec.tool, rec.version)
//     if(!allRequired.every { it in rec.keySet() } ) {
//       log.error """Validation error: required field not set.
//         Required fields: ${allRequired}
//         Offending record: ${rec}"""
//       System.exit 1
//     }
//     if(!allModes.split('\\|').any { it in rec.keySet() }){
//       log.error """Validation error: required field not set.
//         At least one of '${allModes}' is required.
//         Offending record: ${rec}"""
//       System.exit 1
//     }
//     if(!rec.container.contains(rec.version))
//       log.warn "Decalred tool version string ${it.version} not found in container image spec ${rec.container}."
//   }
//   allVersions.each {k, v ->
//     if(v.size()==1)
//       allVersions.put(k,v[0])
//   }
//   return allVersions
// }

// def validateMapperParamsDefinitions (mapperParams, allVersions) {
//   def validationMaps = [dna2dna: [:], rna2dna: [:], rna2rna: [:]]
//   mapperParams.each { //PUT IN DEFAULT VALUES
//     if(!it.containsKey('label')) it.put('label', 'default');
//     if(!it.containsKey('mode')) it.put('mode', allModes);
//     // if(!it.containsKey('version')) it.put('version', 'ALL_AVAILABLE');
//     if(!allVersions.containsKey(it.tool)) {
//       log.error """Validation error:  params defined for undefined mapper, please define in mapperParams
//       Offending record: ${it}"""
//       System.exit 1
//     }
//     if(it.containsKey('version')){
//       if(!it.version.split('\\||,').every {ver -> ver in allVersions."${it.tool}"}) {
//         log.error """Validation error:  params defined for undefined mapper version, please define in mapperParams
//           Offending record: ${it}"""
//         System.exit 1
//       }
//     } else {
//       it.put('version', allVersions."${it.tool}") //NOT SPECIFIED SO PARAM SET APPLIES TO ALL VERSIONS
//     }
//   }.each { rec -> //VALIDATE
//     ['dna2dna', 'rna2dna', 'rna2rna']. each { mode ->
//       validationMap = validationMaps."${mode}"
//       if(mode.matches(rec.mode)) {
//         versions = rec.version instanceof Collection ? rec.version : [rec.version]
//         versions.each{ ver ->
//           key = [rec.tool, ver, rec.label].join("_")
//           stored = validationMap.putIfAbsent(key, rec)
//           if(stored != null) {
//             log.error """Validation error: non-unique label for aligner params set for ${mode}
//             Previously encountered: ${stored}
//             Offending record: ${rec}
//             If tool/version/mode overlap then the labels must be unique"""
//             System.exit 1
//           }
//         }
//       }
//     }
//   }
// }