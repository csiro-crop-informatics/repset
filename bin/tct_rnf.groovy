#!/usr/bin/env groovy

import static groovy.json.JsonOutput.*
import java.util.zip.GZIPInputStream
import java.util.zip.GZIPOutputStream




@Grab('info.picocli:picocli:4.0.0-alpha-3')
@Command(header = [

        $/@|bold,blue  ╔╦╗╔═╗╔╦╗  ╦═╗╔╗╔╔═╗   |@/$,
        $/@|bold,blue   ║ ║   ║   ╠╦╝║║║╠╣    |@/$,
        $/@|bold,blue   ╩ ╚═╝ ╩   ╩╚═╝╚╝╚     |@/$
        ],
        description = "Conver RNF reads coordinates from transcriptome to genome space.",
        version = 'tct_rnf v?.?.?', showDefaultValues = true,
        footerHeading = "%nFor more details, see:%n",
        footer = ["[1] ASCII Art thanks to http://patorjk.com/software/taag/"]
)
@picocli.groovy.PicocliScript
import groovy.transform.Field
import java.security.MessageDigest
import static picocli.CommandLine.*


// @Parameters(arity="1", paramLabel="FILE", description="The file(s) whose checksum to calculate.")
// @Field private File[] files

@Option(names = ["-t", "--transcriptome"], description = ["Transcriptome from `gffread -W` used to simulate input reads."], required = true)
@Field private String transcriptome

@Option(names = ["-f", "--in-forward"], description = ["R1 input file name"], required = true)
@Field private String inforward
@Option(names = ["-r", "--in-reverse"], description = ["R2 input file name"], required = true)
@Field private String inreverse

@Option(names = ["-F", "--out-forward"], description = ["R1 output file name"], required = true)
@Field private String outforward
@Option(names = ["-R", "--out-reverse"], description = ["R2 output file name"], required = true)
@Field private String outreverse

@Option(names= ["-h", "--help"], usageHelp=true, description="Show this help message and exit.")
@Field private boolean helpRequested

@Option(names= ["-V", "--version"], versionHelp=true, description="Show version info and exit.")
@Field private boolean versionInfoRequested

// System.setProperty(picocli.usage.width,160)
// files.each {
//   println MessageDigest.getInstance(algorithm).digest(it.bytes).encodeHex().toString() + "\t" + it
// }



final int PAD4BASES = 9
final int PAD4CHROMOSOMES = 5
final int BUFFER_SIZE = 8192
final String NEWLINE = System.lineSeparator();

File refFile = new File(transcriptome)
File fastqFile1 = new File(inforward)
File fastqFile2 = new File(inreverse)
//File fastqFile1 = new File('A_thaliana_TAIR10_chr1_with_gff_ArtIllumina_reads.1.fq.gz')
//File fastqFile2 = new File('A_thaliana_TAIR10_chr1_with_gff_ArtIllumina_reads.2.fq.gz')
File outFile1 = new File(outforward);
File outFile2 = new File(outreverse);
append = false
writer1 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile1, append)), "UTF-8"), BUFFER_SIZE);
writer2 = new BufferedWriter(new OutputStreamWriter(new GZIPOutputStream(new FileOutputStream(outFile2, append)), "UTF-8"), BUFFER_SIZE);

refs = []
refFile.withReader { source ->
  String line
  while( line = source.readLine()) {
    if(line =~ /^>/) {
      toks = line.split(' ')
      record = [:]
      toks.each { tok ->
        subtoks = tok.split(':|=')
//        println subtoks[1].split('-|,|\\|')
        record.(subtoks[0]) = subtoks[0] == 'loc' ? subtoks[1].split('\\|') : subtoks[1].split(',|\\|')*.split('-')
        if(subtoks[0] =~ /(exons)|(segs)/ ) {
          record.(subtoks[0]) = record.(subtoks[0]).collect{ it*.toInteger() }
        }
//        record.(subtoks[0]) = subtoks[1].split('-|,|\\|')
        //          record.(subtoks[0]) = subtoks[1].split(',')
      }
      //      record.loc = record.loc[0].split("\\|")

//      record.exons = record.exons*.toInteger()
      //          record.exons = record.exons*.split('-')
      //        record.exons = record.exons*.replace('-','..') //.split('-') //*.toInteger()
//      record.segs = record.segs*.toInteger()
      refs << record
    }
  }
}

//println(prettyPrint(toJson(refs[1])))
//
//def record = refs[1]
//
//println record.segs
//println record.exons
//for(int i=0; i<record.segs.size(); i+=2) {
//  slen = record.segs[i+1] - record.segs[i]
//  println record.segs[i] + " - " + record.segs[i+1] + " slen = "+slen
//  elen = record.exons[i+1] - record.exons[i]
//  println record.exons[i] + " - " + record.exons[i+1] + " elen = "+elen
//}



//System.exit 0


try {
  gzipStream1 = new GZIPInputStream(new FileInputStream(fastqFile1), BUFFER_SIZE);
  content1 = new BufferedReader(new InputStreamReader(gzipStream1, "UTF-8"), BUFFER_SIZE);
  gzipStream2 = new GZIPInputStream(new FileInputStream(fastqFile2), BUFFER_SIZE);
  content2 = new BufferedReader(new InputStreamReader(gzipStream2, "UTF-8"), BUFFER_SIZE);
  i = 4
  while ((line1 = content1.readLine()) != null && !line1.isEmpty() && (line2 = content2.readLine()) != null && !line2.isEmpty() ) {
    if(i++ % 4 != 0 ) {
      //JUST OUTPUT/STORE LINE AS IS
      writer1.write(line1);
      writer1.write(NEWLINE);
      writer2.write(line2);
      writer2.write(NEWLINE);
      continue;
    }
//    println line1
//    println line2

    split1 = line1.split('__')

    ////  split2 = line2.split('__')
    //
    //
    coords = split1[2].replaceAll('\\(','').replaceAll('\\)','')
    coordsSplit = coords.split(',')

    //FIELDS TO BE USED
    ref = coordsSplit[1].toInteger()
    //FILEDS TO BE MODIFIED
    start = coordsSplit[3].toInteger()
    end = coordsSplit[4].toInteger()
    startMate = coordsSplit[8].toInteger()
    endMate = coordsSplit[9].toInteger()
    //
    //  println(prettyPrint(toJson(coordsSplit)))
    ////  println "Current ref: "+(ref-1)
    refRecord = refs[ref-1]
//      println(prettyPrint(toJson(refRecord)))
    //
    //CONVERT REF ID
    coordsSplit[1] = coordsSplit[6] = refRecord.loc[0].padLeft(PAD4CHROMOSOMES,'0')

    //CONVERT COORDINATES
    translatedStart = translateAndPad(start, refRecord, PAD4BASES)
    translatedEnd = translateAndPad(end, refRecord, PAD4BASES)
    translatedStartMate = translateAndPad(startMate, refRecord, PAD4BASES)
    translatedEndMate = translateAndPad(endMate, refRecord, PAD4BASES)
    coordsSplit[3] = translatedStart < translatedEnd ? translatedStart : translatedEnd
    coordsSplit[4] = translatedStart >= translatedEnd ? translatedStart : translatedEnd
    coordsSplit[8] = translatedStartMate < translatedEndMate ? translatedStartMate : translatedEndMate
    coordsSplit[9] = translatedStartMate >= translatedEndMate ? translatedStartMate : translatedEndMate
//    coordsSplit[3] = translateAndPad(start, refRecord, PAD4BASES)
//    coordsSplit[4] = translateAndPad(end, refRecord, PAD4BASES)
//    coordsSplit[8] = translateAndPad(startMate, refRecord, PAD4BASES)
//    coordsSplit[9] = translateAndPad(endMate, refRecord, PAD4BASES)

    //    println(prettyPrint(toJson(coordsSplit)))

    //RE-CONSTITUTE READ HEADER
    //  println split1
    //  println(prettyPrint(toJson(split1)))
    StringBuilder sb = new StringBuilder()
    sb.append('(')
    sb.append(coordsSplit[0..4].join(',') )
    sb.append('),(')
    sb.append(coordsSplit[5..9].join(',') )
    sb.append(')')
    split1[2] = sb.toString()
    //  println(prettyPrint(toJson(split1)))
    //  println split1
    //  println line1
    outline = split1.join('__')
    //  println line2
    //  println split1.join('__')[0..-2]+'2'

    println outline

//
//    println '''
//    SHOULD BE: 7807-8260
//    7807-7835
//    intron
//    7942-7987 = 45
//    intron
//    8236-8260
//    '''
//    System.exit 0


    writer1.write(outline);
    writer1.write(NEWLINE);
    writer2.write(outline[0..-2]+'2');
    writer2.write(NEWLINE);
  }
} catch (FileNotFoundException ex) {
  ex.printStackTrace();
} catch (InterruptedException ex) {
  ex.printStackTrace();
} catch (IOException ex) {
  ex.printStackTrace();
} finally {
  try {
    if (writer1 != null) {
      writer1.close();
    }
    if (writer2 != null) {
      writer2.close();
    }
  } catch (IOException ex) {
    ex.printStackTrace();
  }
}


String translateAndPad(int position, Map record, int padding) {
  return translate(position, record).toString().padLeft(padding,'0')
}

int translate (int position, Map record) {
  boolean forward = record.loc[2].equals('+')
  //  println record.exons
  //  assert record.exons.class == java.util.ArrayList
  //  println position
  //  try {

  int numrecords = record.segs.size()

  for(int i=0; i<numrecords; i++) {
//    println i+" "+record.segs
    seg = record.segs[i]
    if(position >= seg[0] && position <= seg[1]) {
      if(forward) {
        return record.exons[i][0]+position-seg[0]
      } else {
//        println 'S '+i+ ', E '+(numrecords-i-1)
        return record.exons[numrecords-i-1][0]+(seg[1]-position)
      }
    }
  }
  return Integer.MIN_VALUE;


//  record.segs.eachWithIndex { seg, i ->
//    println "current segment "+seg
//    println "current position "+position
////    println "current offset1 "+(seg[0]-position)
////    println "current offset2 "+(seg[1]-position)
////    println "current exon "+record.exons[numrecords-i-1]
////    println (position >= seg[0] && position <= seg[1] ? "within" : "without")
//    if(position >= seg[0] && position <= seg[1]) {
//      if(forward) {
//        translated = record.exons[i][0]+position-seg[0]
//      } else {
////        println 'S '+i+ ', E '+(numrecords-i-1)
//        translated = record.exons[numrecords-i-1][0]+(seg[1]-position)
////        return record.exons[i]+position-record.segs[i]
//      }
//    }
//  }
//  return translated

    //  } catch (Exception e) {
  //    println(prettyPrint(toJson(record)))
  //    e.printStackTrace()
  //  }
}

//'''
//SHOULD BE: 7807-8260
//exons:
//0 6788-7069,
//1 7157-7450,
//2 7564-7649,
//3 7762-7835, <<--
//4 7942-7987,
//5 8236-8325, <-
//6 8417-8464,
//7 8571-8737
//
//segs:
//0 1-167,
//1 168-215,
//2 216-305, <-
//3 306-351,
//4 352-425, <<--
//5 426-511,
//6 512-805,
//7 806-1087
//'''