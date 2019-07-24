import static groovy.json.JsonGenerator.*
jsonGenerator = new groovy.json.JsonGenerator.Options()
                .addConverter(java.nio.file.Path) { java.nio.file.Path p, String key -> p.toUriString() }
                .build()

//Input validation specified elswhere
def validators = new GroovyShell().parse(new File("$baseDir/Validators.groovy"))

//Read, parse, validate and sanitize alignment/mapping tools config
def allRequired = ['tool','version','container','index'] //Fields required for each tool in config
def allModes = 'dna2dna|rna2rna|rna2dna' //At leas one mode has to be defined as supported by each tool
def allVersions = validators.validateMappersDefinitions(params.mappersDefinitions, allRequired, allModes)

//Check if specified template files exist
validators.validateTemplatesAndScripts(params.mappersDefinitions, (['index']+(allModes.split('\\|') as List)), '../templates')

//Read, sanitize and validate alignment/mapping param sets
validators.validateMapperParamsDefinitions(params.mapperParamsDefinitions, allVersions)

//Validated now, so gobble up mappers and their params definitions
mappersChannel = Channel.from(params.mappersDefinitions)
  .filter{ params.mappers == 'all' || it.tool.matches(params.mappers) } //TODO Could allow :version
  // .view {groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(it))}
mappersParamsChannel = Channel.from(params.mapperParamsDefinitions)

mapModesChannel = Channel.from(params.mapmode.split('\\||,'))

// println groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(allVersions))
// println groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(params.mapperParams))

/*
Resolve variables emebeded in single-quoted strings
*/
def String resolveScriptVariables(template, binding) {
  def engine = new groovy.text.SimpleTemplateEngine()
  engine.createTemplate(template).make(binding).toString()
}

process index {
  memory '100.MB'
  tag { "${mapper.container}" }
  container { "${mapper.container}" }
  echo true

  input:
    val mapper from mappersChannel

  output:
    val mapper into indices

  exec:
  // println (groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(mapper)))
    def ref = 'reference/path.fa'
    def binding = [ref: ref, task: task.clone()]
  // println bindTemplate(mapper.index, binding) // "${bindTemplate(mapper.index, binding)}"
  script:
    if('index' in mapper.templates) { //Indexing template expected
      template mapper.index == true ? "index/${mapper.tool}.sh" : "index/${mapper.index}"
    } else { //indexing script defined in config
      resolveScriptVariables(mapper.index, binding)
    }
}

// // readsChn = Channel.from(['some_reads', 'some_more_reads'])
// readsChn = Channel.from(['some_reads'])

// process mapReads {
//   maxForks 1
//   tag { "${mapper.tool}@${mapper.version} ${reads} params:${par.label}" }
//   container { "${mapper.container}" }
//   // echo true
//   input:
//      set val(mapper), val(reads), val(mode), val(par) from indices.combine(readsChn).combine(mapModesChannel).combine(mappersParamsChannel)
//     //  set val(mapper), val(ref) from  indices

//   when:
//     mapper.tool == par.tool && mapper.version.matches(par.version.join('|')) && mapper.containsKey(mode) && mode.matches(par.mode)
//   // script:
//   exec:
//     println (groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(mapper)))
//     println mode
//     println (groovy.json.JsonOutput.prettyPrint(jsonGenerator.toJson(par)))
//     def binding = [task: task.clone(), reads: reads]
//     // println "${mapper} ${ref}"
//     // println bindTemplate(mapper.index, binding) // "${bindTemplate(mapper.index, binding)}"
//     //resolveScriptVariables(mapper.dna2dna, binding)
// }

// // // indices.view()





