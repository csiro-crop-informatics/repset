import groovy.json.*
import ch.qos.logback.classic.Logger

def gitHubRelease(Logger log, Map args) {
  GH_TOKEN = System.getenv("GH_TOKEN")

  if(GH_TOKEN == null) {
    log.error "GH_TOKEN not set, unable to create GH release"
    return
  } else {
    log.info "GH_TOKEN found, attempting to create a GH release"
    log.info args.toString()
  }

  //CREATE RELEASE
  responseMap = gitHubApiCall(log, [
    GH_TOKEN: GH_TOKEN, method: 'POST', url:"https://api.github.com/repos/${args.REPO}/releases"  ,
    message: $/{"tag_name": "${args.RELEASE_TAG}", "target_commitish" : "${args.COMMIT}", "name":"${args.RELEASE_NAME}","body":"${args.RELEASE_BODY}", "draft": true }/$
  ])

  if(responseMap == null) {
    return
  }
  //println(JsonOutput.prettyPrint(JsonOutput.toJson(responseMap)))
  relaeseId = responseMap.id

  //UPLOAD ARTEFACT(s)
  uploadUrl = responseMap.upload_url.take(responseMap.upload_url.lastIndexOf('{'))
  //println uploadUrl

  for(fname in args.LOCAL_FILES) {
    File f = new File(fname)
    type = "application/".plus(fname[fname.lastIndexOf('.')+1..-1])

    //Displayed and remote file name could be customized using alternative '?label=' and '?name='
    payload = f.getBytes()
    uploadResponseMap = gitHubApiCall(log, [GH_TOKEN: GH_TOKEN, method: 'POST', url:"${uploadUrl}?name=${f.name}&label=${f.name}", payload: payload, type: type ])
  }

  //FINALIZE RELEASE
  finalResponseMap = gitHubApiCall(log, [GH_TOKEN: GH_TOKEN, method: 'POST', url:"https://api.github.com/repos/${args.REPO}/releases/${relaeseId}", message: $/{"draft": false }/$ ])
  //println(JsonOutput.prettyPrint(JsonOutput.toJson(finalResponseMap)))

}

def gitHubApiCall(Logger log, Map args) {
   HttpURLConnection connection =(HttpURLConnection)  new URL(args.url).openConnection();
    connection.setRequestProperty("Authorization","token ${args.GH_TOKEN}")
    connection.setRequestMethod(args.method)
    connection.setDoOutput(true)
    connection.setRequestProperty("Content-Type", args.type == null ? "application/json" : args.type)
    if(args.containsKey('message')) {
      connection.getOutputStream().write(args.message.getBytes("UTF-8"));
    } else if (args.containsKey('payload')) {
      connection.getOutputStream().write(args.payload);
    } else {
      log.error ('Either JSON control message or payload (file content) is required for the call, terminating')
      connection.disconnect();
      System.exit(1)
    }
    int postRC = connection.getResponseCode();
    log.info(postRC + " [${args.url}]");
    if(postRC >= 200 && postRC <= 202) { //201 == created
        response  = connection.getInputStream().getText()
        jsonSlurper = new JsonSlurper()
        respMap = jsonSlurper.parseText(response)
        return respMap
    } else {
      log.error connection.getErrorStream().getText()
      return null
    }
}
