def addToListInMap (map, key, value) {
  if(map.containsKey(key)) {
    stored = map.get(key)
    if(value in stored) {
      System.err.println """Adding a previously stored value
        Current: ${value}
        Stored: ${stored}"""
      System.exit 1
    }
    stored.add(value)
    map.put(key,stored)
  } else if(value instanceof Collection) {
    map.put(key, value)
  } else {
    map.put(key, [value])
  }
}

def validateMappersDefinitions (mappers, allRequired, allModes) {
  def allVersions = [:] //Keep track of tool versions declared in config
  mappers.each { rec ->
    addToListInMap(allVersions, rec.tool, rec.version)
    rec.each {k, v ->
      if(!(k in (allModes.split('\\|')+allRequired))) {
        System.err.println """Validation error: unexpected field in mapper definition:
          Offending field: ${k}
          Offending record: ${rec}"""
        System.exit 1
      }
    }
    if(!allRequired.every { it in rec.keySet() } ) {
      System.err.println """Validation error: required field not set.
        Required fields: ${allRequired}
        Offending record: ${rec}"""
      System.exit 1
    }
    if(!allModes.split('\\|').any { it in rec.keySet() }){
      System.err.println """Validation error: required field not set.
        At least one of '${allModes}' is required.
        Offending record: ${rec}"""
      System.exit 1
    }
    if(!rec.container.contains(rec.version))
      System.err.println "Warning: decalred tool version string ${rec.version} not found in container image spec ${rec.container}."
  }
  allVersions.each {k, v ->
    if(v.size()==1)
      allVersions.put(k,v[0])
  }
  return allVersions
}

def fileExists(path, rec) {
  // println "validating ${path}"
  if(!new File(path).exists()) {
    System.err.println """Validation error: specified template file does not exist!
    Expected file path: ${path}
    Offending record: ${rec}"""
    System.exit 1
  }
}

def validateTemplatesAndScripts (mappers, keys, path) {
  // println keys
  mappers.each { rec ->
    rec.templates = []
    rec.each { k, v ->
      if(k in keys) { //one of index, dna2dna etc...
        if(v == true) {
          fileExists([path,k,rec.tool+'.sh'].join('/'), rec)
          rec.templates << k
        } else if(v ==~ /^[\w\-_]+.sh$/ ) { //assuming template filename, must end with .sh and no /
          fileExists([path,k,v].join('/'), rec)
          rec.templates << k
        }
      }
    }
  }
}

def validateMapperParamsDefinitions (mapperParams, allVersions) {
  def validationMaps = [dna2dna: [:], rna2dna: [:], rna2rna: [:]]
  mapperParams.each { //PUT IN DEFAULT VALUES
    if(!it.containsKey('label')) it.put('label', 'default');
    if(!it.containsKey('mode')) it.put('mode', allModes);
    // if(!it.containsKey('version')) it.put('version', 'ALL_AVAILABLE');
    if(!allVersions.containsKey(it.tool)) {
      System.err.println """Validation error:  params defined for undefined mapper, please define in mapperParams
      Offending record: ${it}"""
      System.exit 1
    }
    if(it.containsKey('version')){
      it.version = it.version instanceof String ? it.version.split('\\||,') as List : it.version
      if(!it.version.every {ver -> ver in allVersions."${it.tool}"}) {
        System.err.println """Validation error:  params defined for undefined mapper version, please define in mapperParams
          Offending record: ${it}"""
        System.exit 1
      }

    } else {
      it.put('version', allVersions."${it.tool}") //NOT SPECIFIED SO PARAM SET APPLIES TO ALL VERSIONS
    }
  }.each { rec -> //VALIDATE
    ['dna2dna', 'rna2dna', 'rna2rna']. each { mode ->
      validationMap = validationMaps."${mode}"
      if(mode.matches(rec.mode)) {
        versions = rec.version instanceof Collection ? rec.version : [rec.version]
        versions.each{ ver ->
          key = [rec.tool, ver, rec.label].join("_")
          // println "key "+key
          stored = validationMap.putIfAbsent(key, rec)
          if(stored != null) {
            System.err.println """Validation error: non-unique label for aligner params set for ${mode}
            Previously encountered: ${stored}
            Offending record: ${rec}
            If tool/version/mode overlap then the labels must be unique"""
            System.exit 1
          }
        }
      }
    }
  }
}