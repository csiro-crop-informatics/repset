#!/usr/bin/env groovy

import static groovy.json.JsonOutput.*
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream


@Grab('info.picocli:picocli:4.0.0-alpha-3') //command line interface
@Command(header = [
       //Font Name: Calvin S (Caps)
       $/@|bold,blue  ╔═╗╦  ╦╔═╗╦    ╦═╗╔╗╔╔═╗   |@/$,
       $/@|bold,blue  ║╣ ╚╗╔╝╠═╣║    ╠╦╝║║║╠╣    |@/$,
       $/@|bold,blue  ╚═╝ ╚╝ ╩ ╩╩═╝  ╩╚═╝╚╝╚     |@/$
       ],
       description = "Evaluate SAM alignments against original positions RNF-encoded in read names. Paired End only.",
       showDefaultValues = true,
       footerHeading = "%nFootnote(s)%n",
       footer = ["[1] ASCII Art thanks to http://patorjk.com/software/taag/"]
)
@picocli.groovy.PicocliScript
import groovy.transform.Field
import java.security.MessageDigest
import static picocli.CommandLine.*

@Option(names = ["-f", "--faidx"], description = ["FAI index for the reference FASTA file"], required = true)
@Field private String fai

// @Option(names = ["-r", "--rname-match" ], description = "Ignore coordinates, only require reference sequnence name (RNAME) to match.")
// @Field private boolean rnameMatchOnly = false;

@Option(names = ["-d", "--allowed-delta"], description = "Allowed difference between true coordinates and aligned coordinates.")
@Field private int allowedDelta = 5;

@Option(names = ["-s", "--sam"], description = ["Input SAM file name"])
@Field private String sam = '/dev/stdin'

@Option(names = ["-O", "--output"], description = ["Summary output file name"])
@Field private String output = '/dev/stdout'

@Option(names = ["-E", "--es-output"], description = ["ES output file name"])
@Field private String outES

@Option(names= ["-h", "--help"], usageHelp=true, description="Show this help message and exit.")
@Field private boolean helpRequested




// boolean rnameMatchOnly = false


//final int PAD4BASES = 9
//final int PAD4CHROMOSOMES = 5
final int BUFFER_SIZE = 8192
final String NEWLINE = System.lineSeparator();
final String SEP = '\t';

File faiFile = new File(fai)
File samFile = new File(sam)
// File outFile = new File(output)
// File faiFile = new File('Arabidopsis_thaliana_TAIR10.fasta.fai')
// File samFile = new File('filtered.sam')


//append = false
//writer1 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile1, append)), "UTF-8"), BUFFER_SIZE);
//writer2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile2, append)), "UTF-8"), BUFFER_SIZE);

//Parse and store infor generated by gffread when extracting transcripts from genome
def refsList = []
def refsReverseLookup = [:]
//def refsLengthsMap = [:]
faiFile.withReader { source ->
  String line
  int refCount = 0
  while( line = source.readLine()) {
    toks = line.split('\t')
    ref = toks[0]
    // int len = toks[1].toInteger()
    refsList << ref
    refsReverseLookup.put(ref, (++refCount as String))
    //refsLengthsMap.put(ref, len)
  }
}



BufferedWriter esWriter = null;

if(outES != null) {
  File outFileES = new File(outES);
  append = false;
  esWriter = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFileES, append)), "UTF-8"), BUFFER_SIZE);
}



//Process alignments
try {
  if(esWriter != null) {
    esWriter.write """# RN:   read name
# Q:    is mapped with quality
# Chr:  chr id
# D:    direction
# L:    leftmost nucleotide
# R:    rightmost nucleotide
# Cat:  category of alignment assigned by LAVEnder
#         M_i    i-th segment is correctly mapped
#         m      segment should be unmapped but it is mapped
#         w      segment is mapped to a wrong location
#         U      segment is unmapped and should be unmapped
#         u      segment is unmapped and should be mapped
# Segs: number of segments
#
# RN\tQ\tChr\tD\tL\tR\tCat\tSegs"""
  }

  samContentBuffer = new BufferedReader(new InputStreamReader(new FileInputStream(samFile), "UTF-8"), BUFFER_SIZE);

  // def long alignedToCorrectReferenceCount = 0;
  // def wrongReference = new long[256];
  // def unalignedCount = new long[256];
  // def mateAlignedToDifferentReferenceCount = new long[256];

  // def bothMatch = new int[256]
  // def match = new int[256]
  // def mateMatch = new int[256]
  // def mismatch = new int[256]

  def stats = [:]
  ['primary','secondary','supplementary'].each { level ->
    stats."${level}" = [match: [:], mateMatch: [:], mismatch: [:], unaligned: [:], wrongReference: [:], mateReferenceMismatch: [:]]
  }


  while ((line = samContentBuffer.readLine()) != null && !line.isEmpty()) {
    (QNAME, FLAG, RNAME, POS, MAPQ, CIGAR, RNEXT, PNEXT, TLEN, SEQ, QUAL) = line.split('\t')
    def coordsSplit = QNAME.split('__')[2].replaceAll('\\(','').replaceAll('\\)','').split(',')
    int mapq = MAPQ.toInteger()
    int flag = FLAG.toInteger()

    //Reference ID check
    ref = coordsSplit[1].toInteger()
    refRecord = refsList[ref-1]
    boolean isUnaligned = false
    boolean isWrongRef = false

    //MULTI-MAPPER? we want to separate those in reporting
    boolean isSecondary = ((flag & 256) != 0)
    boolean isSupplementary = ((flag & 2048) != 0)
    boolean isPrimary = (!isSecondary && !isSupplementary)
    statsKey = isPrimary ? 'primary' : isSecondary ? 'secondary' : isSupplementary ? 'supplementary' : "ERROR - record should be primary or secondary or supplementary"

    /**
    Treating missing CIGAR as indication of unaligned. This is for consistency with how RNF deals with e.g. kallisto pseudobam.
    We could go back to the idea of only checking if a read aligned to the correct reference in rna2rna mode which could be a sufficient requirement
    e.g. when only interested in quantification.
    **/
    if(RNAME == '*' || CIGAR  == '*') {
      // unalignedCount[mapq]++
      addStat(stats, statsKey, 'unaligned', mapq)
      isUnaligned = true
    } else if(refRecord != RNAME) { //RNAME MATCH
      // wrongReference[mapq]++;
      addStat(stats, statsKey, 'wrongReference', mapq)
      isWrongRef = true
    } else {
      // alignedToCorrectReferenceCount++;
      // mateAlignedToDifferentReferenceCount[mapq] += (RNEXT == '=' || RNEXT == refRecord) ? 0 : 1
      if(RNEXT != '=' && RNEXT != refRecord) {
         addStat(stats, statsKey, 'mateReferenceMismatch', mapq)
      }
    }

    //ORIENTATION
    boolean reverse = ((flag & 16) != 0)

    //Alignment details
    def cigar = CIGAR
    int alnStart = isUnaligned ? 0 : POS.toInteger() //- getStartClip(cigar)
    int alnEnd = isUnaligned ? -1 : alnStart + getAlignmentLength(cigar)-1 //Get aln len from CIGAR string

    category = ''
    if(isUnaligned) {
      category = 'u'
    } else if(isWrongRef) {
      category = 'w'
    } else {
      //RNF encoded coordinates in read name
      int simStart = coordsSplit[3].toInteger()
      int simEnd = coordsSplit[4].toInteger()
      int simStartMate = coordsSplit[8].toInteger()
      int simEndMate = coordsSplit[9].toInteger()

      //Expected v actual
      int startDelta = (alnStart-simStart).abs()
      int endDelta = (alnEnd-simEnd).abs()
      int startDeltaMate = (alnStart-simStartMate).abs()
      int endDeltaMate = (alnEnd-simEndMate).abs()

      if(startDelta <= allowedDelta  && endDelta <= allowedDelta) {
        // match[mapq]++
        addStat(stats, statsKey, 'match', mapq)
        category = 'M_1'
      } else if(startDeltaMate <= allowedDelta && endDeltaMate <= allowedDelta) {
        // mateMatch[mapq]++
        addStat(stats, statsKey, 'mateMatch', mapq)
        category = 'M_2'
      } else {
        // mismatch[mapq]++
        addStat(stats, statsKey, 'mismatch', mapq)
            // matches = stats.get(statsKey).get('mismatch')
            // matches.put(mapq, matches.containsKey(mapq) ? matches.get(mapq)+1 : 1)
        category = 'w'
      }
    }

  //     if(isMatch && isMateMatch) { //THIS DOES HAPPEn BUT RNF JUST RETURNS AFTER FIRST (R1) MATCH
  //       //both reads aligned to same location
  //       bothMatch[mapq]++
  // //      println simStart+ ' '+simEnd+ ' '+simStartMate+ ' '+simEndMate+ ' '+flag+' R1='+isRead1+' R2='+isRead2+' '+ (flag & 64)+ ' ' + (flag & 128) + ' isReverse='+reverse
  // //      println alnStart+ ' '+alnEnd
  // //      println line
  // //      println startDelta+ ' '+endDelta+ ' '+startDeltaMate+' '+endDeltaMate
  // //      println ''
  //       category = 'both-segments-match'
  //     }

      //END RECORD PROCESSING, NOW REPORT
      if(esWriter != null) {
        esWriter.write NEWLINE
        esWriter.write QNAME
        esWriter.write SEP
        esWriter.write isUnaligned ? "unmapped" : "mapped_"+MAPQ
        esWriter.write SEP
        esWriter.write isUnaligned ? ref.toString() : refsReverseLookup.get(RNAME) //ref.toString()  //RNAME //RNF puts the expected ref for unmapped, 'None' would make more sense
        esWriter.write SEP
        esWriter.write reverse ? 'R' : 'F' //should be 'N' if not known
        esWriter.write SEP
        esWriter.write POS
        esWriter.write SEP
        esWriter.write alnEnd < 1 ? "None" : alnEnd.toString()
        esWriter.write SEP
        esWriter.write category
        esWriter.write SEP
        esWriter.write '2' //Currently we're fixed on 2 segs anyway
      }


  }

  //RE-FORMAT DEEPEST LEVEL FROM MAP(MAPQ -> count) TO LIST WITH ["MAPQ" -> value, "count" -> value]
  stats.each { naryKey, naryValue ->
    naryValue.each { category, mapqMap ->
      mapqList = []
      mapqMap.each {mapq, count ->
        mapqList << [MAPQ: mapq, count: count]
      }
      naryValue.put(category,mapqList.sort())
    }
  }

  new File(output).withWriter { out ->
    // out.println "[EDGE-CASES] Both match: "+bothMatch.sum()+" "+bothMatch // This happens but RNF assignes all these to M_1
//     out.println "Match: "+match.sum()+" "+match
//     out.println "Mate_match: "+mateMatch.sum()+" "+mateMatch
//     out.println "Mismatch: "+mismatch.sum()+" "+mismatch
// //  println alignedToCorrectReferenceCount
//     out.println "Wrong_reference: "+wrongReference.sum()+" "+wrongReference
//     out.println "Unaligned: "+unalignedCount.sum()+" "+unalignedCount
//     out.println "Total_reported: "+(unalignedCount.sum()+wrongReference.sum()+match.sum()+mateMatch.sum()+mismatch.sum())
//     out.println "Mate_reference_mismatch: "+mateAlignedToDifferentReferenceCount.count()+" "+mateAlignedToDifferentReferenceCount

    out.println(prettyPrint(toJson(stats)))
    // out.println "Total: "+(unalignedCount.sum()+wrongReference.sum()+bothMatch.sum()+match.sum()+mateMatch.sum()+mismatch.sum())
  }
} catch (FileNotFoundException ex) {
  ex.printStackTrace();
} catch (InterruptedException ex) {
  ex.printStackTrace();
} catch (IOException ex) {
  ex.printStackTrace();
} finally {
  try {
    if (esWriter != null) {
      esWriter.close();
    }
  } catch (IOException ex) {
    ex.printStackTrace();
  }
}


def int getAlignmentLength(String cigar) {
  int len = 0
  StringBuilder segment = new StringBuilder();
  for(int i=0; i < cigar.length(); i++){
    char c=cigar.charAt(i)
    if(Character.isDigit(c)){
      segment.append(c)
    } else {
      // if(c == 'M' || c == '=' || c == 'X' || c == 'D' || c == 'N' || c == 'S' || c =='H') {
      if(c == 'M' || c == '=' || c == 'X' || c == 'D' || c == 'N') {
        len+=segment.toInteger();
      // }else if(c=='I'){
      } else if(c=='I' || c == 'S' || c =='H'){
        //ignore insertions in reference?
        //also ignore soft- and hard-clipping?
      } else if(c=='*'){
        //CIGAR UNAVAILABLE
        return 1 //In some cases kallisto does not provide a CIGAR in pseudobam. Consequently, for a given record we may have alnStart but no alnEnd
      } else{
        System.err.println("Unexpected character ${c} in CIGAR string: ${cigar}")
        System.exit(1)
      }
      segment.setLength(0) //clear buffer
    }
  }
  return len;
}

def int getStartClip(String cigar) {
  def m = (cigar =~ /^[0-9]+[HS]/)
  return m.count == 0 ? 0 : (m[0][0..-2]).toInteger()
}

def int getEndClip(String cigar) {
  def m = (cigar =~ /[0-9]+[HS]$/)
  return m.count == 0 ? 0 : (m[0][0..-2]).toInteger()
}

def addStat(Map stats, String statsKey, String statsCategory, int mapq) {
  matches = stats.get(statsKey).get(statsCategory)
  matches.put(mapq, matches.containsKey(mapq) ? matches.get(mapq)+1 : 1)
}