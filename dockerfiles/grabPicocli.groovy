@Grab('info.picocli:picocli-groovy:4.5.1')
@Command(name = 'justGrabPicocli')
@picocli.groovy.PicocliScript
import groovy.transform.Field
import static picocli.CommandLine.*
println 'done'