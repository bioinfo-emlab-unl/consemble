#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd;

my $version = 0.1;

# Set global variables to be used in the program
my $q = 20;
my $ref = "";
my $time = 0;
my $force = 0;
my $noNorm = 0;
my $normCov = 50;
my $species = "";
my $single = "";
my $reads1 = "";
my $reads2 = "";
my $baseDir = getcwd;
my $verFlag;
my $help = 0;
my $h = 0;
my $pe = 0;
# End: Set global variables to be used in the program

## Get options from command line input and link them to variables
GetOptions(	"assemblyName=s" => \$species,
			"ref=s" => \$ref,
			"q=i" => \$q,
			"time" => \$time,
			"force" => \$force,
			"noNorm" => \$noNorm,
			"normCov" => \$normCov,
			"single=s" => \$single,
			"reads1=s" => \$reads1,
			"reads2=s" => \$reads2,
			"outputDir=s" => \$baseDir,
			"help" => \$help,
			"h" => \$h,
			"version" => \$verFlag
);
# End: Get options from command line input and link them to variables

## Print help message when help option is passed or no single/reads opitions are passed
if ($h or $help or not ($single or ($reads1 and $reads2))) {
	print "Please specify either a single merged file of single end reads (using --single) or two merged files for paired-end reads (using --reads1 and --read2)\n";
	print "\nOther options:\n--species\t\tSpecify the species being assembled to help organize output files\n";
	print "--q\t\t\tChange q-value threshold for reads filtering (default q20)\n";
	print "--baseDir\t\tChange base directory for output files (default is current working directory)\n";
	print "--ref\t\tPath to fasta for reference genome\n";
	print "--force\t\tIgnore the log file and always run entire pipeline\n";
	print "--noNorm\t\tSkip the digital normalization of filtered reads\n";
	print "--normCov\t\tChange the max kmer coverage for digital normalization. Lowering coverage speeds up assembly time, but increases risk of artifacts (default 50)\n";
	print "--time\t\tRecords the amount of time each stage of the pipeline takes\n";
	print "--version\t\tPrints pipeline version and dies\n";
	print "--h or --help\tPrints this message and quits\n\n";
	die;
}
# End: Print help message when help option is passed or no single/reads opitions are passed

## Print consemble version when version option is passed
if ($verFlag) {
	print "Pipeline version $version\n";
	die;
}
# End: Print consemble version when version option is passed

### Begin processing reads
print "Checking read length\n";
my $readLength = 0;
## Process single reads
if ($single) {
	print "Checking single end reads\n";
	$pe = 0;
	my $S;
	open ($S, $single) or die "Can't open single end reads at $single, please check the file location: $!\n"; # Open single read file for reading
	while (my $line = readline($S)) {
		next unless $line =~ /^@/; # Skip unless line begins with "@"
		my $r = readline($S); # Read each line in the single read file
		chomp $r; # Remove trailing spaces
		$readLength = length($r); # Get length of readline
		print "Detected read length of $readLength\n";
		last;
	}
	close $S;
}
# End: Process single reads

## Processing paired-end reads
if ($reads1 and $reads2) {
	print "Checking paired-end reads\n";
	$pe = 1;
	my ($P1, $P2);
	open ($P1, $reads1) or die "Can't open paired end reads 1 at $reads1, please check the file location: $!\n"; # Open first read file for reading

	open ($P2, $reads2) or die "Can't open paired end reads 2 at $reads2, please check the file location: $!\n"; # Open second read file for reading
	my ($rl1, $rl2);
	while (my $line = readline($P1)) {
		next unless $line =~ /^@/; # Skip unless line begins with "@"
		my $r = readline($P1); # Read each line in first paired-end read file
		chomp $r; # Remove trailing spaces
		$rl1 = length($r); # Get length of the readline
		print "Detected read length of $rl1 for reads1\n";
		last;
	}
	while (my $line = readline($P2)) {
		next unless $line =~ /^@/; # Skip unless line begins with "@"
		my $r = readline($P2); # Read each line in second paired-end read file
		chomp $r; # Remove trailing spaces
		$rl2 = length($r); # Get length of the readline
		print "Detected read length of $rl2 for reads2\n";
		last;
	}
	die "Can't determine read length, please check the input files" unless ($rl1 and $rl2);
	if ($rl1 != $rl2) {
		warn "Read lengths in $reads1 and $reads2 do not match.  Please check the files\n";
	}
	$readLength = $rl1; # Set read length from the readline length of one of the files
	close $P1;
	close $P2;
}
# End: Processing paired-end reads

my @q = split /,/, $q; # Split q value from input into an array @q

## Set up output and log directories
my $filReads = "$baseDir/FilteredReads/$species";
my $normReads = "$baseDir/NormalizedReads/$species";
my $asReads = "$baseDir/AssembledReads/$species";
my $runDir = getcwd;
my $log = "$asReads/log.txt";
my ($LOG, $TIME);
my %status;

print "Creating output directories\n";

if (-d "$baseDir") { # Check if working directory exists
	print "Found $baseDir\n";
} else {
	mkdir ("$baseDir"); # Create working directory if it does not exist
	print "Made $baseDir\n";
}

if (-d "$baseDir/AssembledReads/") { # Check if AssembledReads directory exists
	print "Found $baseDir/AssembledReads/\n";
} else {
	mkdir ("$baseDir/AssembledReads/"); # Create AssembledReads directory if it does not exist
	print "Made $baseDir/AssembledReads/\n";
}

if (-d "$baseDir/AssembledReads/$species/") { # Check if AssembledReads/<species> directory exists
	print "Found $baseDir/AssembledReads/$species/\n";
} else {
	mkdir ("$baseDir/AssembledReads/$species/"); # Create AssembledReads/<specie> directory if it does not exist
	print "Made $baseDir/AssembledReads/$species/\n";
}
# End: Set up output and log directories

## Open log file for appending
print "Opening log\n";
if (-e $log and $force != 1) {
	open ($LOG, "$log") or die "$log exists, but can't be read: $!\n";
	while (my $line = readline($LOG)) {
		chomp $line;
		next unless $line;
		
		my ($pr, $st) = split /\t/, $line;
		$status{$pr} = $st;
	}
	close ($LOG) or die "Can't close $log: $!\n";
	open ($LOG, ">>$log") or die "Can't append $log: $!\n";
} else {
	open ($LOG, ">$log") or die "Can't write to $log: $!\n";
}
# End: Open log file for appending

## Check previously completed jobs
my %finished;
print "Checking previously completed jobs\n";
foreach my $pr (sort {$a cmp $b} keys %status) {
	my @info = split /\-/, $pr;
	my $con;
	if ($info[2]) {
		$con = "$info[1]-$info[2]";
	} else {
		$con = $info[1];
	}
	
	if ($finished{$info[0]}) {
		$finished{$info[0]} .= ",$con";
	} else {
		$finished{$info[0]} = $con;
	}
}
# End: Check previously completed jobs

# Open time file
if ($time) {
	open ($TIME, ">$asReads/time.txt");
}
# End: Open time file

## Filter reads using ErneFilter
my $efStartTime = time;
foreach my $q (@q) {
	if ($status{"erneFilter-q$q"}) {
		print "Previously finished erne-filter-q$q, skipping this run\n";
	} else {	
		if (-d "$filReads") { # Check if FilteredReads directory exists
			print "Found $filReads\n";
		} else {
			mkdir ("$filReads"); # Create FilteredReads directory if it does not exist
			print "Made $filReads\n";
		}
		if (-d "$filReads/q$q") { # Check if FilteredReads/q<quality> directory exists
			print "Found $filReads/q$q\n";
		} else {
			mkdir ("$filReads/q$q"); # Create FilteredReads/q<quality> directory if it does not exist
			print "Made $filReads/q$q\n";
		}
		if ($pe) {
			# Run erne-flter on paired-end reads using the given options
			system("erne-filter --query1 $reads1 --query2 $reads2 --output-prefix $filReads/q$q/merged --ultra-sensitive --force-standard");
		} else {
			# Run erne-flter on single reads using the given options
			system("erne-filter --query1 $single --output-prefix $filReads/q$q/merged --ultra-sensitive --force-standard");
		}
		print $LOG "erneFilter-q$q\t1\n"; # Write to log file
	}
}
my $efEndTime = time; # Set end time
my $efRunTime = $efEndTime - $efStartTime; # Calculate runtime
print $TIME "ErneFilter:\t$efRunTime\n" if $time;
# End: Filter reads using ErneFilter

## Normalize unfiltered reads using khmer
my $khStartTime = time; # Set start time
my $finK = $finished{"khmer"}; #Set key for finished job entry
$finK = "none" unless $finK;
foreach my $q (@q) {		
	if (-d "$normReads") { # Check if NormalizedReads directory exists
		print "Found $normReads\n";
	} else {
		mkdir ("$normReads"); # Create NormalizedReads directory if it does not exist
		print "Made $normReads\n";
	}
}

if ($finK =~ /q$q/) {
	print "Previously finished khmer, skipping this run\n";
} else {
	if ($pe) {
		print "Interleaving reads\n";
		# Interleave reads using python library
		system("interleave-reads.py  -o $filReads/q$q/mergedI.fastq --force $filReads/q$q/merged_1.fastq $filReads/q$q/merged_2.fastq");	
		print "Normalizing reads\n";
		if ($noNorm) {
			system("cp $filReads/q$q/mergedI.fastq $normReads/mergedI.fq");	
		} else {
			# Normalize reads
			system("normalize-by-median.py -o $normReads/mergedI.fq -x 1000000000 -C $normCov --n_tables 99 --force --ksize 32 $filReads/q$q/mergedI.fastq");	
		}
		print "Splitting reads\n";
		chdir("$normReads/");
		system("extract-paired-reads.py $normReads/mergedI.fq"); # Extract interleaved paired-end reads into separate files
		system("split-paired-reads.py --force $normReads/mergedI.fq.pe"); # Split interleaved paired-end reads into separate files
		system("mv $normReads/mergedI.fq.pe.1 $normReads/merged.fq.1");
		system("mv $normReads/mergedI.fq.pe.2 $normReads/merged.fq.2");
		system("mv $normReads/mergedI.fq.pe $normReads/mergedI.fq");
	} else {
		if ($noNorm) {
			system("cp $filReads/q$q/merged_1.fastq $normReads/merged.fq");
		} else {
			# Normalize reads
			system("normalize-by-median.py -o $normReads/merged.fq -x 1000000000 -C 50 --n_tables 99 --force $filReads/q$q/merged_1.fastq");
		}
		system("fastq-to-fasta.py -o $normReads/merged.fa $normReads/merged.fq"); # Convert fastq to fasta file
	}
	print $LOG "khmer-q$q\t1\n"; # Write to log file
}
my $khEndTime = time; # Set end time
my $khRunTime = $khEndTime - $khStartTime; # Calculate runtime
print $TIME "Khmer:\t$khRunTime\n" if $time;
chdir("$runDir/");
# End: Normalize unfiltered reads using khmer

my $soStartTime = time;
#Prepare reference
system("bowtie2-build $ref $ref");
# Assembly using TOphat 2
system("tophat2 -o $baseDir/AssembledReads/$species/ -p 8 $ref $normReads/merged.fq.1 $normReads/merged.fq.2");

my %nucl;
# Assembly using Bayesembler
system("bayesembler -b $baseDir/AssembledReads/$species/accepted_hits.bam -o $baseDir/AssembledReads/$species/bayesembler-");
readAssembly("$baseDir/AssembledReads/$species/bayesembler-assembly", "bayesembler"); # Read Bayesembler assembly
# Assembly using Cufflinks
system("cufflinks -o $baseDir/AssembledReads/$species/cufflinks $baseDir/AssembledReads/$species/accepted_hits.bam");
readAssembly("$baseDir/AssembledReads/$species/cufflinks/transcripts", "cufflinks"); # Read Cufflinks assembly
# Assembly using Scallop
system("scallop -i $baseDir/AssembledReads/$species/accepted_hits.bam -o $baseDir/AssembledReads/$species/scallop.gtf");
readAssembly("$baseDir/AssembledReads/$species/scallop", "scallop"); # Read Scallop assembly
# Assembly using Stringtie
system("stringtie -o $baseDir/AssembledReads/$species/stringtie.gtf $baseDir/AssembledReads/$species/accepted_hits.bam");
readAssembly("$baseDir/AssembledReads/$species/stringtie", "stringtie"); # Read Stringtie assembly

my ($OUT, $CONNUCL, $CONNUCLONG);
# Open consemble3+g.faa for writing
open ($OUT, ">$baseDir/AssembledReads/$species/consemble3+g.faa") or die "Can't open $baseDir/AssembledReads/$species/consemble3+g.faa: $!";
# Open consemble3+g.fa for writing
open ($CONNUCL, ">$baseDir/AssembledReads/$species/consemble3+g.fa") or die "Can't open $baseDir/AssembledReads/$species/consemble3+g.faa: $!";
# Open consemble3+g_longest.fa for writing
open ($CONNUCLONG, ">$baseDir/AssembledReads/$species/consemble3+g_longest.fa") or die "Can't open $baseDir/AssembledReads/$species/consemble3+g.faa: $!";
my $i = 1;
my %seqNames;
my %seqCounts;

## Counting and matching consensus sequences
foreach my $seq (keys %seqCounts) {
	next unless $seq;
	if ($seqCounts{$seq} > 2) { # Check if sequence at $seqCounts{$seq} appears more than twice
		my $shortest = 999999; # Initiate $shortest sequence variable
		my $longest = 0; # Initiate longest sequence variable
		my @assem = split /,/, $seqNames{$seq}; # Split sequece ids and store in @assem array
		my ($nuclSeq, $nuclSeqLong);
		foreach my $s (@assem) { # Loop throught @assem array of ids belonging to $seq sequence being processed in current iteration
			#print "$s\n";
			if (length($nucl{$s}) < $shortest) { # Check if sequence at $nucl{$s} is shorter than $shortest
				$nuclSeq = $nucl{$s}; # Set $nuclSeq to sequence at $nucl{$s}
				$shortest = length($nucl{$s}); # Set value of $shortest to length of $nucl{$s}
			}
			if (length($nucl{$s}) > $longest) { # Check if sequence at $nucl{$s} is longer than $longest
				$nuclSeqLong = $nucl{$s}; # Set $nuclSeqLong to sequence at $nucl{$s}
				$longest = length($nucl{$s}); # Set value of $longest to length of $nucl{$s}
			}
		}
		print $OUT ">Contig_$i\n$seq\n"; # Wtite to consemble3+g.faa file
		print $CONNUCL ">Contig_".$i."\n$nuclSeq\n"; # Write to consemble3+g.fa file
		print $CONNUCLONG ">Contig_".$i."\n$nuclSeqLong\n"; # Write to consemble3+g_longest.fa file
		$i++;
	}
}
# End: Counting and matching consensus sequences

## Subroutine readAssembly
sub readAssembly {
	my ($f, $a) = @_; # Set assembly file $f and assembly method $a parameters
	system("gtf_to_fasta $f.gtf $ref $f.fa"); # Convert fasta reference to gtf
	readNucAssembly($f, $a);
	system("ORFfinder -in $f.fa -out $f-ORFs.fa"); # Run ORFinder on assembly file
	getLongestORF($f);
	my $F;
	open ($F, "$f-ORFs.fa.faa") or die "Can't open $f-ORFs.fa.faa: $!\n"; # Open $f-ORFs.fa.faa, product of getLongestORF($f), for reading
	my %seen;
	my $seq = "";
	my $id;
	while (my $line = readline($F)) {
		chomp $line; # Remove trailing spaces
		next unless $line; # skip if read line is empty
		
		if ($line =~ /^>/) { # Check if read line is a label
			if (!$seen{$seq}) { # Check if value at $seen{$seq} is false
				$seqCounts{$seq}++; # Increment value at $seqCounts{$seq}
				$seen{$seq}; # Set $seen{$seq} to true
			}
			$line =~ s/^>//; # Remove ">" from label
			my @info = split /\s+/, $line; # Split label into an array @info
			if ($seqNames{$seq}) { # Check if value at $seqNames{$seq} is set
				$seqNames{$seq} .= ",$id" if $seq; # Concatenate value at $seqNames{$seq} with ",$id", from previous iteration, if $seq is set
			} else {
				$seqNames{$seq} = "$id" if $seq; # Set $seqNames{$seq} to $id if $seq is set
			}
			$seq = ""; # Set $seq to empty string
			$info[0] =~ s/\_\d+//; # Replace "_<any-digit>" in first item of @info
			$id = "$a-$info[0]"; # Set $id to a concatenation of $a and first element of $info. To be used in next iteration
			print "$id\n";
		} else {
			$seq .= $line; # Concatenate read line to sequence to create one long string 
		}
	}
	if (!$seen{$seq}) { # Check if value at $seen{$seq} is false
		$seqCounts{$seq}++; # Increment value at $seqCounts{$seq}
		$seen{$seq}; # Set $seen{$seq} to true
	}
	if ($seqNames{$seq}) { # Check if value at $seqNames{$seq} is set
		$seqNames{$seq} .= ",$id"; # Concatenate value at $seqNames{$seq} with ",$id", from last iteration in while loop
	} else {
		$seqNames{$seq} = "$id"; # Set $seqNames{$seq} to $id from last iteration in while loop
	}
	$seq = ""; # Set $seq to empty string
}
# End: Subroutine readAssembly

## Subroutine readNucAssembly
sub readNucAssembly {
        my ($f, $a) = @_; # Set assembly file $f and assembly method $a parameters

        my $seq = "";
        my $id;
		my $F;
		open ($F, "$f.fa") or die "Can't open $f.fa: $!\n"; # Open assembly file for reading
        while (my $line = readline($F)) {
                chomp $line; # Remove trailing spaces
                next unless $line; # Skip if read line is empty

                if ($line =~ /^>/) { # Check if read line is label
                    $line =~ s/^>//; # Remove ">" at teh beginning of label
                    my @info = split /\s+/, $line; # Split label in array @info
                    $nucl{$id} = $seq if $seq; # If $seq is set, from previous iteration, set $nucl to that $seq where id = $id 
                    $seq = ""; # Set $seq to empty string
                    $id = "$a-$info[0]"; # Set $id to a concatenation of $a and first element of $info. To be used in next iteration
					print "$id\n";
                } else {
                    $seq .= $line; # Concatenate read line to sequence to create one long string
                }
        }
        $nucl{$id} = $seq if $seq; # If $seq is set, from last iteration in while loop, set $nucl to that $seq where id = $id 
        $seq = ""; # Set $seq to empty string
}
# End: Subroutine readNucAssembly

## Subroutine getLongestORF
sub getLongestORF {
	my $f = shift @_; # Set $f to first item in passed arguments
	my %tmpGenes;
	my ($tmpId, $tmpSeq);
	my ($TMP, $TMPOUT);
	open ($TMP, "$f-ORFs.fa") or die "Can't open$f-ORFs.fa\n"; # Open $f-ORFs.fa for reading
	open ($TMPOUT, ">$f-ORFs.fa.faa") or die "Can't open $f-ORFs.fa.faa\n"; # Open $f-ORFs.fa.faa for writing

	while (my $line = readline($TMP)) {
		chomp $line; # Remove trailing spaces
		next unless $line; # Skip if read line is not set
		if ($line =~ /^>/) { # check if read line is a label
			if ($tmpSeq) { # Check if $tmpSeq is set
				$tmpGenes{$tmpId} = "" unless $tmpGenes{$tmpId}; # Set $tmpGenes{$tmpId} to empty unless it already has a value
				# Set $tmpGenes{$tmpId} to $tmpSeq, from previous iteration, if it's longer than the current value at $tmpGenes{$tmpId}
				$tmpGenes{$tmpId} = $tmpSeq if (length($tmpSeq) > length($tmpGenes{$tmpId}));
			}
			my @info = split /ORF\d+_/, $line; # Split label into an array @info
			my @gn = split /:/, $info[1]; # Split first item of @info into array @gn
			$tmpId = $gn[0]; # Set $tempId to first item in @gn
			$tmpSeq = ""; # Set $tempSeq to empty string
		} else {
			$tmpSeq .= $line; # Concatenate read line to sequence to create one long string
		}
	}
	# Set $tmpGenes{$tmpId} to $tmpSeq, from last iteration in while loop, if it's longer than the current value at $tmpGenes{$tmpId}
	$tmpGenes{$tmpId} = $tmpSeq if (length($tmpSeq) > length($tmpGenes{$tmpId}));

	foreach my $g (sort {$a cmp $b} keys %tmpGenes) {
		print $TMPOUT ">$g\n$tmpGenes{$g}\n";
	}

	close ($TMP);
	close ($TMPOUT);
}
# End: Subroutine getLongestORF