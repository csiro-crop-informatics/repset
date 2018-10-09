require 'optparse'
require "erubis"
require 'logger'
path = File.expand_path(File.dirname(__FILE__))
require "#{path}/logging"
include Logging

#####
#
#   Runs the statistics for a given dataset
#   IN: dataset_name source_of_tree
#   out:
#   1) Sorted and appropriate sam files
#   2) runs compare2truth
#   3) runs compare junctions
#
####

# 2015/8/10 Katharina Hayer

$logger = Logger.new(STDERR)

# Initialize logger
def setup_logger(loglevel)
  case loglevel
  when "debug"
    $logger.level = Logger::DEBUG
  when "warn"
    $logger.level = Logger::WARN
  when "info"
    $logger.level = Logger::INFO
  else
    $logger.level = Logger::ERROR
  end
end

def setup_options(args)
  options = {
    :loglevel => "error",
    :debug => false,
    :read_length => nil,
    :cut_bases => 0,
    :single_end => false
  }

  opt_parser = OptionParser.new do |opts|
    opts.banner = "\nUsage: ruby compare2truth.rb [options] truth.cig sorted.sam"
    opts.separator ""
    opts.separator "truth.cig:"
    opts.separator "seq.1a  chr10 123502684 123502783 100M  123502684-123502783 - CTTAAGTATGGGGAAGGTAGAAAGTTCATTTCATTACTTATAAAATATGTCTTCTCAAGAACAAAACTGTGCTGTTACAACTCAGTGTTCAATGTGAAAT"
    opts.separator "seq.1b  chr10 123502530 123502629 100M  123502530-123502629 + AGGGAGTACTATATTCTAGGGGAAAAAACTATGCCAAGACAACAGACATGAACAGGACTGTCCTGAACAAATGGATTCCTGATGCTAACACAAGCTCCAT"
    opts.separator ""
    opts.separator "sorted.sam:"
    opts.separator "seq.1a  83  chr10 123502684 255 100M  = 123502530 -254  ATTTCACATTGAACACTGAGTTGTAACAGCACAGTTTTGTTCTTGAGAAGACATATTTTATAAGTAATGAAATGAACTTTCTACCTTCCCCATACTTAAG  * NH:i:1  HI:i:1  AS:i:192  nM:i:3"
    opts.separator "seq.1b  163 chr10 123502530 255 100M  = 123502684 254 AGGGAGTACTATATTCTAGGGGAAAAAACTATGCCAAGACAACAGACATGAACAGGACTGTCCTGAACAAATGGATTCCTGATGCTAACACAAGCTCCAT  * NH:i:1  HI:i:1  AS:i:192  nM:i:3"
    opts.separator ""
    # enumeration
    #opts.on('-a', '--algorithm ENUM', [:all, :contextmap2,
    #  :crac, :gsnap, :hisat, :mapsplice2, :olego, :rum,
    #  :star,:soap,:soapsplice, :subread, :tophat2],'Choose from below:','all: DEFAULT',
    #  'contextmap2','crac','gsnap','hisat', 'mapsplice2',
    #  'olego','rum','star','soap','soapsplice','subread','tophat2') do |v|
    #  options[:algorithm] = v
    #end

    opts.on("-d", "--debug", "Run in debug mode") do |v|
      options[:log_level] = "debug"
      options[:debug] = true
    end

    #opts.on("-o", "--out_file [OUT_FILE]",
    #  :REQUIRED,String,
    #  "File for the output, Default: overview_table.xls") do |anno_file|
    #  options[:out_file] = anno_file
    #end

    opts.on("-r", "--read_length [INT]",
      :REQUIRED,Integer,
      "read length, if not specified it will be taken from cig file") do |s|
      options[:read_length] = s
    end

    opts.on("-c", "--cut_bases [INT]",
      :REQUIRED,Integer,
      "cut bases of cig file") do |s|
      options[:cut_bases] = s
    end

    opts.on("-s", "--single_end", "Run in single_end mode") do |v|
      options[:single_end] = true
    end

    opts.on("-v", "--verbose", "Run verbosely") do |v|
      options[:log_level] = "info"
    end

    opts.separator ""
  end

  args = ["-h"] if args.length == 0
  opt_parser.parse!(args)
  setup_logger(options[:log_level])
  if args.length != 2
    $logger.error("You only provided #{args.length} fields, but 2 required!")
    raise "Please specify the input (truth.cig sorted.sam)"
  end
  options
end

class Stats
  def initialize()
    @total_number_of_bases_of_reads = 0
    @total_number_of_reads = 0
    @total_number_of_bases_aligned_correctly = 0
    @total_number_of_reads_aligned_correctly = 0
    @total_number_of_bases_aligned_incorrectly = 0
    @total_number_of_reads_aligned_incorrectly = 0
    @total_number_of_bases_aligned_ambiguously = 0
    @total_number_of_reads_aligned_ambiguously = 0
    @total_number_of_bases_unaligned = 0
    @total_number_of_reads_unaligned = 0
    @total_number_of_bases_in_true_insertions = 0
    @total_number_of_bases_in_true_deletions = 0
    @total_number_of_bases_in_true_skipping = 0
    @total_number_of_bases_in_true_skipping_binary = 0
    @total_number_of_reads_in_true_skipping_binary = 0
    @total_number_of_bases_called_insertions = 0
    @total_number_of_bases_called_deletions = 0
    @total_number_of_bases_called_skipped = 0
    @total_number_of_bases_called_skipped_binary = 0
    @insertions_called_correctly = 0
    @deletions_called_correctly = 0
    @skipping_called_correctly = 0
    @skipping_called_correctly_binary = 0
    # Sides can be one of "none", "left", "right", "ambiguous" or "both"
    @skipping_sides = [0,0,0,0]
  end

  attr_accessor :total_number_of_bases_of_reads,
    :total_number_of_reads,
    :total_number_of_bases_aligned_correctly,
    :total_number_of_reads_aligned_correctly,
    :total_number_of_bases_aligned_incorrectly,
    :total_number_of_reads_aligned_incorrectly,
    :total_number_of_bases_aligned_ambiguously,
    :total_number_of_reads_aligned_ambiguously,
    :total_number_of_bases_unaligned,
    :total_number_of_reads_unaligned,
    :total_number_of_bases_in_true_insertions,
    :total_number_of_bases_in_true_deletions,
    :total_number_of_bases_in_true_skipping,
    :total_number_of_bases_in_true_skipping_binary,
    :total_number_of_reads_in_true_skipping_binary,
    :total_number_of_bases_called_insertions,
    :total_number_of_bases_called_deletions,
    :total_number_of_bases_called_skipped,
    :total_number_of_bases_called_skipped_binary,
    :insertions_called_correctly,
    :deletions_called_correctly,
    :skipping_called_correctly,
    :skipping_called_correctly_binary,
    :skipping_sides

  def to_s
    %{total_number_of_bases_of_reads: #{@total_number_of_bases_of_reads}
total_number_of_reads: #{@total_number_of_reads}
total_number_of_bases_aligned_correctly: #{@total_number_of_bases_aligned_correctly}
total_number_of_reads_aligned_correctly: #{@total_number_of_reads_aligned_correctly}
total_number_of_bases_aligned_incorrectly: #{@total_number_of_bases_aligned_incorrectly}
total_number_of_reads_aligned_incorrectly: #{@total_number_of_reads_aligned_incorrectly}
total_number_of_bases_aligned_ambiguously: #{@total_number_of_bases_aligned_ambiguously}
total_number_of_reads_aligned_ambiguously: #{@total_number_of_reads_aligned_ambiguously}
total_number_of_bases_unaligned: #{@total_number_of_bases_unaligned}
total_number_of_reads_unaligned: #{@total_number_of_reads_unaligned}
total_number_of_bases_in_true_insertions: #{@total_number_of_bases_in_true_insertions}
total_number_of_bases_in_true_deletions: #{@total_number_of_bases_in_true_deletions}
total_number_of_bases_in_true_skipping: #{@total_number_of_bases_in_true_skipping}
total_number_of_bases_in_true_skipping_binary: #{@total_number_of_bases_in_true_skipping_binary}
total_number_of_bases_called_insertions: #{@total_number_of_bases_called_insertions}
total_number_of_bases_called_deletions: #{@total_number_of_bases_called_deletions}
total_number_of_bases_called_skipped: #{@total_number_of_bases_called_skipped}
total_number_of_bases_called_skipped_binary: #{@total_number_of_bases_called_skipped_binary}
insertions_called_correctly: #{@insertions_called_correctly}
deletions_called_correctly: #{@deletions_called_correctly}
skipping_called_correctly: #{@skipping_called_correctly}
skipping_called_correctly_binary: #{@skipping_called_correctly_binary}
skipping_sides: #{@skipping_sides.join(":")}}
  end

  def fill_skipping_sides(word, num=1)
    #Sides can be one of "none", "left", "right", "ambiguous" or "both"
    case word
    when "none"
      @skipping_sides[0] += num
    when "left"
      @skipping_sides[1] += num
    when "right"
      @skipping_sides[2] += num
    else
      @skipping_sides[3] += num
    end

  end

  def process
    # READ LEVEL
    out = "--------------------------------------\n"
    out += "total_number_of_reads:\t#{@total_number_of_reads}\n"
    percent_reads_aligned_correctly = (@total_number_of_reads_aligned_correctly.to_f / @total_number_of_reads.to_f * 10000).to_i / 100.0
    out += "accuracy over all reads:\t#{percent_reads_aligned_correctly}%\n"
    total_num_unique_aligners = @total_number_of_reads_aligned_correctly + @total_number_of_reads_aligned_incorrectly
    #$logger.debug("total_num_unique_aligned_reads=#{total_num_unique_aligners}")
    if total_num_unique_aligners == 0
      accuracy_on_unique_aligners = 0
    else
      accuracy_on_unique_aligners = (@total_number_of_reads_aligned_correctly.to_f / total_num_unique_aligners.to_f * 10000).to_i / 100.0
    end
    ##print "% unique aligners correct:\t$accuracy_on_unique_aligners%\n";
    out += "accuracy over uniquely aligned reads:\t#{accuracy_on_unique_aligners}%\n"
    percent_reads_aligned_incorrectly = (@total_number_of_reads_aligned_incorrectly.to_f / @total_number_of_reads.to_f * 10000.0).to_i / 100.0
    ##print "total_number_of_bases_aligned_incorrectly = $total_number_of_bases_aligned_incorrectly\n";
    out += "% reads aligned incorrectly:\t#{percent_reads_aligned_incorrectly}%\n"
    percent_reads_aligned_ambiguously = (@total_number_of_reads_aligned_ambiguously.to_f / @total_number_of_reads.to_f * 10000).to_i / 100.0
    ##print "total_number_of_bases_aligned_ambiguously = $total_number_of_bases_aligned_ambiguously\n";
    out += "% reads aligned ambiguously:\t#{percent_reads_aligned_ambiguously}%\n"
    percent_reads_unaligned = (@total_number_of_reads_unaligned.to_f / @total_number_of_reads.to_f * 10000).to_i / 100.0
    ##print "total_number_of_bases_unaligned = $total_number_of_bases_unaligned\n";
    out += "% reads unaligned:\t#{percent_reads_unaligned}%\n"
    percent_reads_aligned = 100 - percent_reads_unaligned
    out += "% reads aligned:\t#{percent_reads_aligned}%\n"

    intron_rate = (@total_number_of_reads_in_true_skipping_binary.to_f / @total_number_of_reads.to_f * 1000000).to_i / 10000.0
    out += "% of reads with true introns:\t#{intron_rate}%\n"


    # BASE LEVEL
    out += "--------------------------------------\n"
    out += "total_number_of_bases_of_reads:\t#{@total_number_of_bases_of_reads}\n"
    percent_bases_aligned_correctly = (@total_number_of_bases_aligned_correctly.to_f / @total_number_of_bases_of_reads.to_f * 10000).to_i / 100.0
    out += "accuracy over all bases:\t#{percent_bases_aligned_correctly}%\n"
    total_num_unique_aligners = @total_number_of_bases_aligned_correctly + @total_number_of_bases_aligned_incorrectly
    $logger.debug("total_num_unique_aligners=#{total_num_unique_aligners}")
    if total_num_unique_aligners == 0
      accuracy_on_unique_aligners = 0
    else
      accuracy_on_unique_aligners = (@total_number_of_bases_aligned_correctly.to_f / total_num_unique_aligners.to_f * 10000).to_i / 100.0
    end
    ##print "% unique aligners correct:\t$accuracy_on_unique_aligners%\n";
    out += "accuracy over uniquely aligned bases:\t#{accuracy_on_unique_aligners}%\n"
    percent_bases_aligned_incorrectly = (@total_number_of_bases_aligned_incorrectly.to_f / @total_number_of_bases_of_reads.to_f * 10000.0).to_i / 100.0
    ##print "total_number_of_bases_aligned_incorrectly = $total_number_of_bases_aligned_incorrectly\n";
    out += "% bases aligned incorrectly:\t#{percent_bases_aligned_incorrectly}%\n"
    percent_bases_aligned_ambiguously = (@total_number_of_bases_aligned_ambiguously.to_f / @total_number_of_bases_of_reads.to_f * 10000).to_i / 100.0
    ##print "total_number_of_bases_aligned_ambiguously = $total_number_of_bases_aligned_ambiguously\n";
    out += "% bases aligned ambiguously:\t#{percent_bases_aligned_ambiguously}%\n"
    percent_bases_unaligned = (@total_number_of_bases_unaligned.to_f / @total_number_of_bases_of_reads.to_f * 10000).to_i / 100.0
    ##print "total_number_of_bases_unaligned = $total_number_of_bases_unaligned\n";
    out += "% bases unaligned:\t#{percent_bases_unaligned}%\n"
    percent_bases_aligned = 100 - percent_bases_unaligned
    out += "% bases aligned:\t#{percent_bases_aligned}%\n"
    #puts "number of bases in true insertions = #{@total_number_of_bases_in_true_insertions}"
    insertion_rate = (@total_number_of_bases_in_true_insertions.to_f / @total_number_of_bases_of_reads.to_f * 1000000).to_i / 10000.0
    out += "% of bases in true insertions:\t#{insertion_rate}%\n"
    deletion_rate = (@total_number_of_bases_in_true_deletions.to_f / @total_number_of_bases_of_reads.to_f * 1000000).to_i / 10000.0
    out += "% of bases in true deletions:\t#{deletion_rate}%\n"


    # INSERTIONS DELETIONS SKIPPING

    if(@total_number_of_bases_in_true_insertions==0)
      out += "insertions FN/FD rate:\tNo insertions exist in true data.\n"
    else
      if(@total_number_of_bases_called_insertions>0)
        #false_discovery_rate
        insertions_false_discovery_rate = ((1 - (@insertions_called_correctly.to_f / @total_number_of_bases_called_insertions.to_f * 10000).to_i / 10000.0) * 100 * 10000).to_i/10000.0
        out += "insertions FD rate:\t#{insertions_false_discovery_rate}%\n"
      else
        out += "insertions FD rate:\t0% (no insertions called)\n"
      end
      #false_negative_rate
      insertions_false_negative_rate = ((1 - (@insertions_called_correctly.to_f / @total_number_of_bases_in_true_insertions.to_f * 10000).to_i / 10000.0) * 100* 10000).to_i/10000.0
      out += "insertions FN rate:\t#{insertions_false_negative_rate}%\n"
    end

    if(@total_number_of_bases_in_true_deletions==0)
      out += "deletions FN/FD rate:\tNo deletions exist in true data.\n"
    else
      if(@total_number_of_bases_called_deletions>0)
        #false_discovery_rate
        deletions_false_discovery_rate = ((1 - (@deletions_called_correctly.to_f / @total_number_of_bases_called_deletions.to_f * 10000).to_i / 10000.0) * 100 * 10000).to_i/10000.0
        out += "deletions FD rate:\t#{deletions_false_discovery_rate}%\n"
      else
        out += "deletions FD rate:\t0% (no deletions called)\n"
      end
      #false_negative_rate
      deletions_false_negative_rate = ((1 - (@deletions_called_correctly.to_f / @total_number_of_bases_in_true_deletions.to_f * 10000).to_i / 10000.0) * 100 * 10000).to_i/10000.0
      out += "deletions FN rate:\t#{deletions_false_negative_rate}%\n"
    end

    if(@total_number_of_bases_in_true_skipping==0)
      out += "skipping FN/FD rate:\tNo skipping exist in true data.\n"
    else
      if(@total_number_of_bases_called_skipped>0)
        #false_discovery_rate
        skipping_false_discovery_rate = ((1 - (@skipping_called_correctly.to_f / @total_number_of_bases_called_skipped.to_f * 10000).to_i / 10000.0) * 100 * 10000).to_i/10000.0
        out += "skipping FD rate:\t#{skipping_false_discovery_rate}%\n"
      else
        out += "skipping FD rate:\t0% (no skipping called)\n"
      end
      #false_negative_rate
      skipping_false_negative_rate = ((1 - (@skipping_called_correctly.to_f / @total_number_of_bases_in_true_skipping.to_f * 10000).to_i / 10000.0) * 100 * 10000).to_i/10000.0
      out += "skipping FN rate:\t#{skipping_false_negative_rate}%\n"
    end

    #if(@total_number_of_bases_in_true_skipping_binary==0)
    #  out += "skipping FN/FD rate:\tNo skipping exist in true data.\n"
    #else
    #  if(@total_number_of_bases_called_skipped_binary>0)
    #    #false_discovery_rate
    #    skipping_false_discovery_rate = ((1 - (@skipping_called_correctly_binary.to_f / @total_number_of_bases_called_skipped_binary.to_f * 10000).to_i / 10000.0) * 100 * 10000).to_i/10000.0
    #    out += "skipping FD rate:\t#{skipping_false_discovery_rate}%\n"
    #  else
    #    out += "skipping FD rate:\t0% (no skipping called)\n"
    #  end
    #  #false_negative_rate
    #  skipping_false_negative_rate = ((1 - (@skipping_called_correctly_binary.to_f / @total_number_of_bases_in_true_skipping_binary.to_f * 10000).to_i / 10000.0) * 100 * 10000).to_i/10000.0
    #  out += "skipping FN rate:\t#{skipping_false_negative_rate}%\n"
    #end
    # JUNCTIONS LEVEL
    out += "--------------------------------------\n"


    if(@total_number_of_bases_in_true_skipping_binary==0)
      out += "junctions FN/FD rate:\tNo skipping exist in true data.\n"
    else
      if(@total_number_of_bases_called_skipped_binary>0)
        #false_discovery_rate
        skipping_false_discovery_rate = ((1 - (@skipping_called_correctly_binary.to_f / @total_number_of_bases_called_skipped_binary.to_f * 10000).to_i / 10000.0) * 100 * 10000).to_i/10000.0
        out += "junctions FD rate:\t#{skipping_false_discovery_rate}%\n"
      else
        out += "junctions FD rate:\t0% (no junctions called)\n"
      end
      #false_negative_rate
      skipping_false_negative_rate = ((1 - (@skipping_called_correctly_binary.to_f / @total_number_of_bases_in_true_skipping_binary.to_f * 10000).to_i / 10000.0) * 100 * 10000).to_i/10000.0
      out += "junctions FN rate:\t#{skipping_false_negative_rate}%\n"
    end
    out += "Junctions Sides (none|left|right|both):\t#{@skipping_sides.join("|")}\n"
    if total_number_of_bases_called_skipped_binary > 0
      out += "Junctions Sides (none|left|right|both)% of all called:\t#{@skipping_sides.map { |e| "#{(((e.to_f/total_number_of_bases_called_skipped_binary.to_f * 10000).to_i / 10000.0) * 100 * 10000).to_i/10000.0}%"}.join("|")}\n"
    else
      out += "Junctions Sides (none|left|right|both)% of all called:\tNaN|NaN|NaN|NaN"
    end
    out
  end
end

class MappingObject
  def initialize()
    # Current Pos 100
    # 95M [100,195]
    @matches = []
    # 3I [100,3]
    @insertions = []
    # 4D [100,104]
    @deletions = []
    # 123N [100,223]
    @skipped = []
    # 30S/H [100,30]
    @unaligned = []
  end

  attr_accessor :matches,
    :insertions,
    :deletions,
    :skipped,
    :unaligned


  def melt_ranges()
    $logger.debug("HERE #{matches.length} #{matches}")
    if matches.length > 1
      stops = []
      starts = []
      matches.each do |m|
        starts << m[0]
        stops << m[1]
      end
      $logger.debug("matches: #{starts.join("|")}; #{stops.join("|")}")
      starts.each_with_index do |s,i|
        next if i == 0
        if s-stops[i-1] == 1
          $logger.debug("GOT HERE")
          starts.delete_at(i)
          stops.delete_at(i-1)
        end
      end
      $logger.debug("matches: #{starts.join("|")}; #{stops.join("|")}")
      
      matches_new = []
      starts.each_with_index do |s,i|
        matches_new << [s,stops[i]]
      end
      self.matches = matches_new
    end
  end

  def to_s
    %{Matches: #{matches.join(":")},
Insertions: #{insertions.join(":")},
Deletions: #{deletions.join(":")},
Skipped: #{skipped.join(":")},
Unaligned: #{unaligned.join(":")}
}
  end


  #def fix
  #  @matches = comp(@matches)
  #  @insertions = comp(@insertions)
  #  @deletions = comp(@deletions)
  #  @skipped = comp(@skipped)
  #  @unaligned = comp(@unaligned)
  #end
#
  #private
#
  #def comp(some)
  #  out = some.dup
  #  out.flatten!
  #  out.each_with_index do |t1, i|
  #    next unless i.odd?
  #    next unless out[i+1]
  #    puts "MOMA"
  #    if t1 == out[i+1]
  #      out.delete_at(i)
  #      out.delete_at(i)
  #    end
  #  end
  #  out
  #end

end

def cut_adapters(cig_group,num_cut_bases)
  $logger.debug("NiNA #{num_cut_bases}")
  cig_group_new = []

  cig_group.each do |line|
    $logger.debug(line)
    fields = line.split("\t")

    cig_cigar_nums = fields[4].split(/\D/).map { |e|  e.to_i }
    cig_cigar_letters = fields[4].split(/\d+/).reject { |c| c.empty? }
    starts = fields[5].split(", ").map { |chr| chr.split("-")[0].to_i  }
    ends = fields[5].split(", ").map { |chr| chr.split("-")[1].to_i  }
    $logger.debug starts.length
    $logger.debug ends
    cig_cigar_nums_new = cig_cigar_nums.dup
    cig_cigar_letters_new = cig_cigar_letters.dup
    starts_dup = starts.dup
    ends_dup = ends.dup

    if fields[0] =~ /a$/ 
      #forward = true
      count = 0
      (cig_cigar_nums.length-1).downto(0) do |i|
        e = cig_cigar_nums[i]
        if ["D","I","S","H","N"].include?(cig_cigar_letters[i])
          cig_cigar_nums_new.delete_at(i)
          cig_cigar_letters_new.delete_at(i)
          count -= 1
          next
        end
        if e < num_cut_bases
          num_cut_bases = num_cut_bases - e
          cig_cigar_letters_new.delete_at(i)
          starts_dup.delete_at(i)
          ends_dup.delete_at(i)
        else
          cig_cigar_nums_new[i] = cig_cigar_nums[i]-num_cut_bases
          ends_dup[count] = ends_dup[count]-num_cut_bases
          break
        end
        count += 1
      end 
    else
      count = 0
      cig_cigar_nums.each_with_index do |e,i|
        if ["D","I","S","H","N"].include?(cig_cigar_letters[i])
          cig_cigar_nums_new.delete_at(i)
          cig_cigar_letters_new.delete_at(i)
          count -= 1
          #count -= 1 if ["D","I","S","H"].include?(cig_cigar_letters[i])
          next
        end
        if e < num_cut_bases
          num_cut_bases = num_cut_bases - e
          cig_cigar_letters_new.delete_at(i)
          starts_dup.delete_at(i) unless ["D","I","S","H"].include?(cig_cigar_letters[i+1])
          ends_dup.delete_at(i) unless ["D","I","S","H"].include?(cig_cigar_letters[i+1])
        else
          cig_cigar_nums_new[i] = cig_cigar_nums[i]-num_cut_bases
          $logger.debug "count #{count}"
          starts_dup[count] = starts_dup[count]+num_cut_bases
          break
        end
        count += 1
      end

    end
    new_cig = ""
    cig_cigar_nums_new.each_with_index do |e,i|
      new_cig += "#{e}#{cig_cigar_letters_new[i]}"
    end
    fields[4] = new_cig
    new_cig_range = []
    starts_dup.each_with_index do |e,i|
      new_cig_range << "#{e}-#{ends_dup[i]}"
    end
    fields[5] = new_cig_range.join(", ")
    fields[2] = starts_dup[0]
    fields[3] = ends_dup[-1]
    cig_group_new << fields.join("\t")
  end
  cig_group_new
end

def files_valid?(truth_cig,sam_file,options)
  l = `grep ^seq #{truth_cig} | head -1`
  l.chomp!
  l =~ /seq.(\d+)/
  first_truth = $1
  l =~ /\t([^\t]+)$/;
  options[:read_length] ||= $1.length;
  l = `tail -1 #{truth_cig}`
  l.chomp!
  l =~ /seq.(\d+)/
  last_truth = $1
  l = `grep ^seq #{sam_file} | head -1`
  l.chomp!
  l =~ /seq.(\d+)/
  first_sam = $1
  l = `tail -1 #{sam_file}`
  l.chomp!
  l =~ /seq.(\d+)/
  last_sam = $1
  unless last_sam == last_truth && first_sam == first_truth && last_sam &&  first_sam
    $logger.error("Sam file and cig file don't start and end in the same sequence!")
    $logger.debug("last_sam #{last_sam}, last_truth #{last_truth}")
    $logger.debug("first_sam #{first_sam}, first_truth #{first_truth}")
    raise "both files must start and end with the same sequence number and must have an entry for every sequence number in between."
  end
end

def fill_mapping_object(mo, start, cigar_nums, cigar_letters)
  current_pos = start
  add = 0
  first = true
  cigar_nums.each_with_index do |num,i|
    first = false if i > 0
    case cigar_letters[i]
    when "M"
      mo.matches << [current_pos, current_pos + num + add]
      current_pos += num #+ add
      add = 0
    when "I"
      mo.insertions << [current_pos, current_pos + num]
      add = 1 #num
    when "D"
      mo.deletions << [current_pos, current_pos + num]
      current_pos += num
    when "N"
      mo.skipped << [current_pos, current_pos + num]
      current_pos += num
    when "H","S"
      if first
        $logger.debug("MEEEA")
        mo.unaligned << [start-num, start]
      else
        mo.unaligned << [current_pos, current_pos + num]
      end
      #current_pos += num
    end
  end
end

# Returns [ #matches, #misaligned]
def compare_ranges(true_ranges, inferred_ranges, insertion_mode = false)
  matches = 0
  misaligned = 0
  # Sides can be one of "none", "left", "right", "ambiguous" or "both"
  sides = []
  true_ranges.each_with_index do |t1, i|
    next unless i.even?
    t2 = true_ranges[i+1]
    inferred_ranges.each_with_index do |i1, k|
      next unless k.even?
      old_matches = matches
      old_misaligned = misaligned
      i2 = inferred_ranges[k+1]
      if t1 <= i1 && t2 >= i2
        matches += (i2 - i1)
        side = "none"
        side = "left" if t1 == i1
        side = "right" if t2 == i2
        side = "both" if t2 == i2 && t1 == i1
        sides << side
      elsif !insertion_mode && t1 <= i1 && i1 < t2 && t2 <= i2
        matches += (t2 - i1)
        misaligned += i2 - t2
        if t1 == i1
          sides << "left"
        else
          sides << "none"
        end
      elsif insertion_mode && t1 <= i1 && i1 <= t2 && t2 <= i2
        #puts "YOUNK"
        matches += (t2 - i1)
        misaligned += i2 - t2
      #elsif t1 >= i1  && t2 <= i2
      #  matches += (t2 - t1)
      #  misaligned += (i2 - t2) + (t1 - i1)
      #  $logger.debug "BUBBLES"
      elsif !insertion_mode && t1 >= i1  && t2 >= i2 && t1 < i2
        matches += (i2 - t1)
        misaligned += (t1 - i1)
        if t2 == i2
          sides << "right"
        else
          sides << "none"
        end
      elsif insertion_mode && t1 >= i1  && t2 >= i2 && t1 <= i2
        #puts "YOUNS"
        matches += (i2 - t1)
        misaligned += (t1 - i1)

      end
      $logger.debug "Matches #{matches}"
      $logger.debug "Misaligned #{misaligned}"
      if matches != old_matches || misaligned != old_misaligned
        inferred_ranges.delete_at(k)
        inferred_ranges.delete_at(k)
        break
      end
      #puts misaligned
    end
  end

  inferred_ranges.each_with_index do |i1, k|
    next unless k.even?
    i2 = inferred_ranges[k+1]
    misaligned += i2-i1
    sides << "none"
  end
  if matches < 0 || misaligned < 0
    puts matches
    puts misaligned
    exit
  end
  $logger.debug("SIDES #{sides}")
  [matches, misaligned, sides]
end

def fix_cigar(t_nums,t_letters,i_nums,i_letters)
  #$logger.debug t_nums.join("T")
  #$logger.debug t_letters.join("T")
  #$logger.debug i_nums.join("I")
  #$logger.debug i_letters.join("I")
  t_nums.each_with_index do |t_num, i|
    next if t_num == i_nums[i]
    case t_letters[i]
    when 'M'
      if ['N','I','D'].include?(t_letters[i+1])
        if t_nums[i+1] == i_nums[i+1] && i_nums[i+2] && (i_nums[i]-t_num).abs == (i_nums[i+2]-t_nums[i+2]).abs
          i_nums[i] = t_num
          i_nums[i+2] = t_nums[i+2]
        end
      end
    end
  end
  #$logger.debug "LALA"
  #$logger.debug t_nums.join("T")
  #$logger.debug t_letters.join("T")
  #$logger.debug i_nums.join("I")
  #$logger.debug i_letters.join("I")
end

# Returns [#matches, #misaligned]
#def compare_pos_count(true_pos_count, inferred_pos_count)
#  correct = 0
#  incorrect = 0
#  true_pos_count.each_with_index do |t_pos, i|
#    next unless i.even?
#    t_count = true_pos_count[i+1]
#    inferred_pos_count.each_with_index do |i_pos, k|
#      next unless k.even?
#      old_correct = correct
#      i_count = inferred_pos_count[k+1]
#      if t_pos <= i_pos &&
#      if old_correct != correct
#        inferred_pos_count.delete_at(k)
#        inferred_pos_count.delete_at(k)
#      end
#      #puts misaligned
#    end
#  end
#  inferred_pos_count.each_with_index do |i1, k|
#    next unless k.even?
#    i2 = inferred_pos_count[k+1]
#    incorrect += i2-i1
#  end
#
#  [correct, incorrect]
#end

def exists?(file)
  out = File.exist?(file) &&
  !File.zero?(file)
  raise "File #{file} does not exist or is empty!" unless out
  out
end

def comp_base_by_base(s_sam,c_cig,stats,skipping_length,skipping_binary,options)
  $logger.debug(s_sam.join("::"))
  $logger.debug(c_cig.join("::"))
  cig_cigar_nums = c_cig[4].split(/\D/).map { |e|  e.to_i }
  cig_cigar_letters = c_cig[4].split(/\d+/).reject { |c| c.empty? }
  sam_cigar_nums = s_sam[5].split(/\D/).map { |e|  e.to_i }
  sam_cigar_letters = s_sam[5].split(/\d+/).reject { |c| c.empty? }

  c_cig_mo = MappingObject.new()
  fill_mapping_object(c_cig_mo, c_cig[2].to_i, cig_cigar_nums, cig_cigar_letters)
  c_cig_mo.melt_ranges()
  $logger.debug(c_cig_mo)

  #c_cig_mo.fix
  s_sam_mo = MappingObject.new()
  if (cig_cigar_letters & ["I","D","N"]).length > 0 && (sam_cigar_letters & ["I","D","N"]).length > 0 &&
    cig_cigar_letters == sam_cigar_letters
    # In case I, D or N is ambigous
    $logger.debug("I,D,N ambigous")
    fix_cigar(cig_cigar_nums,cig_cigar_letters,sam_cigar_nums,sam_cigar_letters)
  end
  fill_mapping_object(s_sam_mo, s_sam[3].to_i, sam_cigar_nums, sam_cigar_letters)
  #s_sam_mo.fix
  s_sam_mo.melt_ranges()
  $logger.debug(s_sam_mo)
  # How many matches?
  $logger.debug("MATCHES")
  matches_misaligned = compare_ranges(c_cig_mo.matches.flatten, s_sam_mo.matches.flatten)
  $logger.debug("matches_misaligned #{matches_misaligned.join("|")}")
  $logger.debug("DELETIONS")
  deletions_incorrect = compare_ranges(c_cig_mo.deletions.flatten, s_sam_mo.deletions.flatten,true)
  stats.deletions_called_correctly += deletions_incorrect[0]
  stats.total_number_of_bases_called_deletions += deletions_incorrect[1] + deletions_incorrect[0]
  stats.total_number_of_bases_aligned_incorrectly += deletions_incorrect[1]
  
  if matches_misaligned[0]  > options[:read_length]
    matches_misaligned[0] = options[:read_length]  
  end
  if matches_misaligned[0] + matches_misaligned[1] > options[:read_length]
    
    matches_misaligned[1] = options[:read_length] - matches_misaligned[0]
    #puts matches_misaligned
    #STDIN.gets
  end
  stats.total_number_of_bases_aligned_correctly += matches_misaligned[0] - deletions_incorrect[1]
  stats.total_number_of_bases_aligned_incorrectly += matches_misaligned[1]

  if matches_misaligned[0] > 0
    stats.total_number_of_reads_aligned_correctly += 1
    if matches_misaligned[0] != options[:read_length]
      #cur_sum = options[:read_length] - matches_misaligned[1]- matches_misaligned[0]
      stats.total_number_of_bases_unaligned += options[:read_length] - matches_misaligned[1]- matches_misaligned[0] #if cur_sum > 0
    end
  else
    stats.total_number_of_reads_aligned_incorrectly += 1

    if matches_misaligned[1] > 0
      #cur_sum = options[:read_length] - matches_misaligned[1]
      stats.total_number_of_bases_unaligned += options[:read_length] - matches_misaligned[1] #if cur_sum > 0
    end

  end


  # Insertions
  $logger.debug("INSERTIONS")
  insertions_incorrect = compare_ranges(c_cig_mo.insertions.flatten, s_sam_mo.insertions.flatten,true)
  $logger.debug("insertions_incorrect #{insertions_incorrect.join('|')}" )
  stats.insertions_called_correctly += insertions_incorrect[0]
  stats.total_number_of_bases_called_insertions += insertions_incorrect[1] + insertions_incorrect[0]
  #stats.total_number_of_bases_aligned_incorrectly += insertions_incorrect[1]
  # Deletions
  # Skipping
  $logger.debug("SKIPPING")
  skipping_incorrect = compare_ranges(c_cig_mo.skipped.flatten, s_sam_mo.skipped.flatten)
  stats.skipping_called_correctly += skipping_incorrect[0]
  stats.total_number_of_bases_called_skipped += skipping_incorrect[1] + skipping_incorrect[0]
  $logger.debug(skipping_incorrect)
  #if skipping_incorrect[0]-skipping_length == 0 &&  skipping_incorrect[1] == 0 && skipping_incorrect[0] > 0
  #  stats.skipping_called_correctly_binary += skipping_incorrect[2].length
  #  #skipping_incorrect[2] = "both"
  #end
  #stats.total_number_of_bases_called_skipped_binary += skipping_incorrect[2].length if skipping_incorrect[0] > 0 || skipping_incorrect[1] > 0 #|| (skipping_incorrect[0]-skipping_length).abs > 10
  if  skipping_incorrect[1] > 0 || skipping_incorrect[0]  > 0
    skipping_incorrect[2].each do |e|
      stats.fill_skipping_sides(e)
      stats.skipping_called_correctly_binary += 1 if e == "both"
      stats.total_number_of_bases_called_skipped_binary += 1
    end
  end
  # How many clippings?
  #$logger.debug("CLIPPING")
  #unaligned = compare_ranges(c_cig_mo.unaligned.flatten, s_sam_mo.unaligned.flatten)
  #$logger.debug(unaligned.join(",,,"))
  #stats.total_number_of_bases_unaligned += unaligned[1]
  #stats.total_number_of_bases_aligned_incorrectly += matches_misaligned[1]
  #puts unaligned
end

def process(current_group, cig_group, stats,options)
  $logger.debug("PROCESSING")
  if options[:cut_bases] > 0
    cig_group = cut_adapters(cig_group,options[:cut_bases])
  end
  if options[:single_end]
    stats.total_number_of_reads += 1
  else
    stats.total_number_of_reads += 2
  end
  cig_group.each do |l|
    l = l.split("\t")
    k = l[4].dup
    inserts = 0
    while k =~ /(\d+)I/
      inserts = inserts+$1.to_i
      k.sub!(/(\d+)I/,"")
    end
    stats.total_number_of_bases_in_true_insertions += inserts
    k = l[4].dup
    deletions = 0
    while k =~ /(\d+)D/
      deletions = deletions+$1.to_i
      k.sub!(/(\d+)D/,"")
    end
    stats.total_number_of_bases_in_true_deletions += deletions
    k = l[4].dup
    skipping = 0
    skipping_binary = 0
    while k =~ /(\d+)N/
      skipping = skipping+$1.to_i
      skipping_binary += 1
      k.sub!(/(\d+)N/,"")
    end
    stats.total_number_of_bases_in_true_skipping += skipping
    if skipping > 0
      stats.total_number_of_bases_in_true_skipping_binary += skipping_binary
      stats.total_number_of_reads_in_true_skipping_binary += 1
    end
    stats.total_number_of_bases_of_reads += options[:read_length]
    if current_group.length > 2 || (current_group.length > 1 && options[:single_end])
      ##### HERE MULTIMAPPER ROUTINE!
      all_unaligned = true
      current_group.each do |s|
        s = s.split("\t")
        all_unaligned = false unless s[2] == "*" || s[5] == "*"
      end
      stats.total_number_of_bases_aligned_ambiguously += options[:read_length] unless all_unaligned
      stats.total_number_of_reads_aligned_ambiguously += 1 unless all_unaligned
      if all_unaligned
        stats.total_number_of_bases_unaligned += options[:read_length]
        stats.total_number_of_reads_unaligned += 1
      end
    else
      current_group.each do |s|
        s = s.split("\t")
        next unless l[0] == s[0]
        if s[2] == "*" || s[5] == "*"
          stats.total_number_of_bases_unaligned += options[:read_length]
          stats.total_number_of_reads_unaligned += 1
        else
          if s[2] != l[1]
            stats.total_number_of_bases_aligned_incorrectly += options[:read_length]
            stats.total_number_of_reads_aligned_incorrectly += 1
          else
            if s[3] == l[2] && s[5] == l[4]
              stats.total_number_of_bases_aligned_correctly += options[:read_length]
              stats.insertions_called_correctly += inserts
              stats.total_number_of_bases_called_insertions += inserts
              stats.deletions_called_correctly += deletions
              stats.total_number_of_bases_called_deletions += deletions
              stats.skipping_called_correctly += skipping
              stats.total_number_of_bases_called_skipped += skipping
              #if skipping > 0
              stats.skipping_called_correctly_binary += skipping_binary
              stats.total_number_of_bases_called_skipped_binary += skipping_binary
              stats.fill_skipping_sides("both",skipping_binary)
              #end
              stats.total_number_of_reads_aligned_correctly += 1
            else
              $logger.debug("SKIPPING_LENGTH #{skipping}")
              comp_base_by_base(s,l,stats,skipping,skipping_binary,options)
            end
          end
        end
      end
    end
  end
  $logger.debug(current_group.length)
  $logger.debug(cig_group[0])
end

def compare(truth_cig, sam_file, options)
  stats = Stats.new()
  $logger.debug(stats)
  truth_cig_handler = File.open(truth_cig)
  sam_file_handler = File.open(sam_file)
  current_group = []
  cig_group = []
  current_num = nil
  current_letter = nil
  count = 0
  while !sam_file_handler.eof?
    # process one sequence name at a time
    line = sam_file_handler.readline.chomp
    next unless line =~ /^seq/
    line =~ /seq.(\d+)/
    current_num ||= $1
    now_num = $1
    line =~ /seq.\d+(a|b)/
    current_letter ||= $1
    now_letter = $1
    $logger.debug("current_num: #{current_num}, current_letter #{current_letter}")
    $logger.debug("current_num: #{current_num}, current_letter #{current_letter}")
    if current_num == now_num && !options[:single_end]
      current_group << line
    elsif current_num == now_num && options[:single_end] && current_letter == now_letter 
      current_group << line
    else
      cig_group << truth_cig_handler.readline.chomp
      cig_group << truth_cig_handler.readline.chomp unless options[:single_end]
      count += 1
      $logger.debug(cig_group)
      if (count % 50000 == 0)
        STDERR.puts "finished #{count} reads"
      end
      #STDIN.gets
      process(current_group, cig_group,stats,options)
      current_num = now_num
      current_letter = now_letter
      current_group = []
      cig_group = []
      current_group << line
    end
  end

  cig_group << truth_cig_handler.readline.chomp
  cig_group << truth_cig_handler.readline.chomp unless options[:single_end]
  process(current_group, cig_group,stats,options)
  stats
end

def run(argv)
  options = setup_options(argv)
  truth_cig = argv[0] if exists?(argv[0])
  sam_file = argv[1] if exists?(argv[1])
  $logger.info("Options are #{options}")

  files_valid?(truth_cig,sam_file,options)
  stats = compare(truth_cig, sam_file, options)
  $logger.info(stats)
  puts stats.process()

  $logger.info("All done!")
end

if __FILE__ == $0
  run(ARGV)
end
