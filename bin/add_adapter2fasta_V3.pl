#!/usr/bin/env perl

# Script to add retained adapter tails into read sequences.
# 75% of all reads will have retained adapters with lengths ranging from 1 to a maximum of 30bp
# These retained adaptors will be inserted into the read sequence at the 5' end with the read 3' trimmed back to it's original length
# The actual length of retained adapters will be randomly chosen with 50% in length range 1bp to 5bp, 25% in length range 6bp to 10bp, 12.5% in length range 11bp to 15bp, and remander 16bp to 30bp
# Will be using the Illumina Universal Adapter as the retained adaptor template sequence
# Usage:
# add_adapter2fasta_V3.pl read_forward.fa read_reverse.fa read_forward_adapters.fa read_reverse_adapters.fa

# Illumina Universal Adapter
$adapter = "AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT";
chomp($adapter);
$adapter = uc($adapter);

# reverse complement of the adapter sequence (for reverse reads)
$adapter_rev = reversecomplement($adapter);

# for reproducabilty always explicity seed the random number generator
srand(int(5637));	# haven't given much thought to this seed, it's a prime so should be good enough for this application!

$forward_out = $ARGV[2];
$reverse_out = $ARGV[3];
open(FORWARD, ">$forward_out");

open(INFILE, $ARGV[0]); # the forward reads fasta file
$flag = 0;
$count = 0;
while($line = <INFILE>) {
    if($flag == 0) {
	print FORWARD $line;
	$flag = 1;
	next;
    }
    chomp($line);
	$line = uc($line);
	$SeqLen=length($line);
    $flag = 0;
	$RetainLen=int(rand(4));	# 75% of reads to contain at least 1bp of retained adapter
	if($RetainLen > 0) {		# > 0 if read to have retained adapter
		$RetainLen=int(rand(10))+1;	# 50% of reads with retained adapters to retain 1..5bp
		if($RetainLen > 5) {
			$RetainLen = int(rand(10))+6; # 25% of reads with retained adapters to retain 6..10bp
			if($RetainLen > 10) {
				$RetainLen = int(rand(11))+11; # 12.5% of reads with retained adapters to retain 11..15bp
				if($RetainLen > 15) {
					$RetainLen = int(rand(15) + 16);	# remainder (12.5%) to retain 16..30bp
					}
				}

			}
		$line = substr($adapter,-1 * $RetainLen) . substr($line,0,$SeqLen-$RetainLen);
		}
    print FORWARD "$line\n";
    $count++;
	}
close(INFILE);
close(FORWARD);

open(REVERSE, ">$reverse_out");
open(INFILE, $ARGV[1]); # the reverse reads fa file
$flag = 0;
$count = 0;
while($line = <INFILE>) {
    if($flag == 0) {
	print REVERSE $line;
	$flag = 1;
	next;
    }
    chomp($line);
	$line = uc($line);
	$SeqLen=length($line);
    $flag = 0;
	$RetainLen=int(rand(4));	# 75% of reads to contain at least 1bp of retained adapter
	if($RetainLen > 0) {		# > 0 if read to have retained adapter
		$RetainLen=int(rand(10))+1;	# 50% of reads with retained adapters to retain 1..5bp
		if($RetainLen > 5) {
			$RetainLen = int(rand(10))+6; # 25% of reads with retained adapters to retain 6..10bp
			if($RetainLen > 10) {
				$RetainLen = int(rand(11))+11; # 12.5% of reads with retained adapters to retain 11..15bp
				if($RetainLen > 15) {
					$RetainLen = int(rand(15) + 16);	# remainder (12.5%) to retain 16..30bp
					}
				}

			}
		$line = substr($adapter_rev,-1 * $RetainLen) . substr($line,0,$SeqLen-$RetainLen);
		}
    print REVERSE "$line\n";
    $count++;
	}


sub getrandbase () {
    $x = int(rand(4));
    if($x == 0) {
        return "A";
    }
    if($x == 1) {
        return "C";
    }
    if($x == 2) {
        return "G";
    }
    if($x == 3) {
        return "T";
    }
}

sub reversecomplement () {
    ($sq) = @_;
    @A = split(//,$sq);
    $rev = "";
    for($i=@A-1; $i>=0; $i--) {
        $flag = 0;
        if($A[$i] eq 'A') {
            $rev = $rev . "T";
            $flag = 1;
        }
        if($A[$i] eq 'T') {
            $rev = $rev . "A";
            $flag = 1;
        }
        if($A[$i] eq 'C') {
            $rev = $rev . "G";
            $flag = 1;
        }
        if($A[$i] eq 'G') {
            $rev = $rev . "C";
            $flag = 1;
        }
        if($flag == 0) {
            $rev = $rev . $A[$i];
        }
    }

    return $rev;
}

