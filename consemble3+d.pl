#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Cwd;
use File::Path qw(make_path remove_tree);

my $version = 0.1;

# Set global variables to be used in the program
my $q = 20;
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
	open ($S, $single) or die "Can't open single end reads at $single, please check the file location: $!\n";
	while (my $line = readline($S)) {
		next unless $line =~ /^@/;
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
	open ($P1, $reads1) or die "Can't open paired end reads 1 at $reads1, please check the file location: $!\n";
	open ($P2, $reads2) or die "Can't open paired end reads 2 at $reads2, please check the file location: $!\n";
	my ($rl1, $rl2);
	while (my $line = readline($P1)) {
		next unless $line =~ /^@/;
		my $r = readline($P1); # Read each line in first paired-end read file
		chomp $r; # Remove trailing spaces
		$rl1 = length($r); # Get length of the readline
		print "Detected read length of $rl1 for reads1\n";
		last;
	}
	while (my $line = readline($P2)) {
		next unless $line =~ /^@/;
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

my @q = split /,/, $q; # Split q value from input into an array @

## Set up output and log directories
my $filReads = "$baseDir/FilteredReads/$species";
my $normReads = "$baseDir/NormalizedReads/$species";
my $asReads = "$baseDir/AssembledReads/$species";
my $runDir = getcwd;
my $log = "$baseDir/log.txt";
my ($LOG, $TIME);
my %status;

print "Creating output directories\n";

if (-d "$baseDir") { # Check if working directory exists
	print "Found $baseDir\n";
} else {
	make_path ("$baseDir"); # Create working directory if it does not exist
	print "Made $baseDir\n";
}

if (-d "$baseDir/AssembledReads/") { # Check if AssembledReads directory exists
	print "Found $baseDir/AssembledReads/\n";
} else {
	make_path ("$baseDir/AssembledReads/"); # Create AssembledReads directory if it does not exist
	print "Made $baseDir/AssembledReads/\n";
}

if (-d "$baseDir/AssembledReads/$species/") { # Check if AssembledReads/<species> directory exists
	print "Found $baseDir/AssembledReads/$species/\n";
} else {
	make_path ("$baseDir/AssembledReads/$species/"); # Create AssembledReads/<specie> directory if it does not exist
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
	open ($TIME, ">$baseDir/time.txt");
}
# End: Open time file

## Filter reads using ErneFilter
my $efStartTime = time; # Set start time
foreach my $q (@q) {
	if ($status{"erneFilter-q$q"}) {
		print "Previously finished erne-filter-q$q, skipping this run\n";
	} else {	
		if (-d "$filReads") { # Check if FilteredReads directory exists
			print "Found $filReads\n";
		} else {
			make_path ("$filReads"); # Create FilteredReads directory if it does not exist
			print "Made $filReads\n";
		}
		if (-d "$filReads/q$q") {# Check if FilteredReads/q<quality> directory exists
			print "Found $filReads/q$q\n";
		} else {
			make_path ("$filReads/q$q"); # Create FilteredReads/q<quality> directory if it does not exist
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
		make_path ("$normReads"); # Create NormalizedReads directory if it does not exist
		print "Made $normReads\n";
	}
}

if ($finK =~ /q$q/) {
	print "Previously finished khmer, skipping this run\n";
} else {
	if ($pe) {
		print "Interleaving reads\n";
		system("interleave-reads.py  -o $filReads/q$q/mergedI.fastq --force $filReads/q$q/merged_1.fastq $filReads/q$q/merged_2.fastq");	# Interleave reads using python library
		print "Normalizing reads\n";
		if ($noNorm) {
			system("cp $filReads/q$q/mergedI.fastq $normReads/mergedI.fq");	 
		} else {
			system("normalize-by-median.py -o $normReads/mergedI.fq -x 1000000000 -C $normCov --n_tables 99 --force --ksize 32 $filReads/q$q/mergedI.fastq");	# Normalize reads
		}
		print "Splitting reads\n";
		# chdir("$normReads/");
		system("extract-paired-reads.py $normReads/mergedI.fq"); # Extract interleaved paired-end reads into separate files
		system("split-paired-reads.py --force $normReads/mergedI.fq.pe"); # Split interleaved paired-end reads into separate files
		system("mv $normReads/mergedI.fq.pe.1 $normReads/merged.fq.1");
		system("mv $normReads/mergedI.fq.pe.2 $normReads/merged.fq.2");
		system("mv $normReads/mergedI.fq.pe $normReads/mergedI.fq");
		system("fastq-to-fasta.py -o $normReads/mergedI.fa $normReads/mergedI.fq"); # Convert fastq to fasta file
		system("fastq-to-fasta.py -o $normReads/merged.fa.1 $normReads/merged.fq.1"); # Convert fastq to fasta file
		system("fastq-to-fasta.py -o $normReads/merged.fa.2 $normReads/merged.fq.2"); # Convert fastq to fasta file
	} else {
		if ($noNorm) {
			system("cp $filReads/q$q/merged_1.fastq $normReads/merged.fq");
		} else {
			system("normalize-by-median.py -o $normReads/merged.fq -x 1000000000 -C 50 --n_tables 99 --force $filReads/q$q/merged_1.fastq"); # Normalize reads
		}
		system("fastq-to-fasta.py -o $normReads/merged.fa $normReads/merged.fq"); # Convert fastq to fasta file
	}
	print $LOG "khmer-q$q\t1\n"; # Write to log file
}
my $khEndTime = time; # Set end time
my $khRunTime = $khEndTime - $khStartTime; # Calculate runtime
print $TIME "Khmer:\t$khRunTime\n" if $time;
# chdir("$runDir/");
# End: Normalize unfiltered reads using khmer

## Assembly using SOAPdenovo-Trans assembly
my $soStartTime = time; # Set start time

my $fq = "$normReads/merged.fq";
my $fqI = "$normReads/mergedI.fq";
my $kMin = 15; # Set minimum k-mer value
my $kMax;
if ($readLength < 127) { # Why 127????
	$kMax = $readLength; # Set maximum k-mer value
} else {
	$kMax = 127; # Set maximum k-mer value
}

my $CON;

if (-d "$asReads/SOAP/") { # Check if AssembledReads/SOAP directory exists
	print "Found $asReads/SOAP/\n";
} else {
	make_path ("$asReads/SOAP/"); # Make AssembledReads/SOAP if directory does not exist
	print "Made $asReads/SOAP/\n";
}

my $finS = $finished{"soap"}; # Set finished job key value
$finS = "none" unless $finS;
for (my $k = $kMin; $k <= $kMax; $k += 4) { # Run for k-mer in increaments of 4
	if ($finS =~ /q$q-k$k/) { # Skip run if already completed
		print "Previously finished soap-q$q-k$k, skipping this run\n";
	} else {
		if (-d "$asReads/SOAP/k$k/") { # Check if AssembledReads/SOAP/k<k-mer_value> directory exists
			print "Found $asReads/SOAP/$k/\n";
		} else {
			make_path ("$asReads/SOAP/k$k/"); # Make AssembledReads/SOAP/k<k-mer_value> if directory does not exist
			print "Made $asReads/SOAP/$k/\n";
		}

		open ($CON, ">SOAP_Conf.txt") or die "Can't open configuration file: $!\n"; # Open SOAPdenovo-Trans configuration file for writing
		# Write to configuration file
		print $CON "#maximal read length\nmax_rd_len=$readLength\n[LIB]\n#maximal read length in this lib\nrd_len_cutof=150\n#average insert size\navg_ins=150\n#if sequence needs to be reversed \nreverse_seq=0\n#in which part(s) the reads are used\nasm_flags=3\n#minimum aligned length to contigs for a reliable read location (at least 32 for short insert size)\nmap_len=32\n";
		if ($pe) {
			print $CON "#fastq file for read 1\nq1=$fq.1\nq2=$fq.2\n"; # Set reads to paired-end reads in configuration file
		} else {
			print $CON "#fastq file for read 1\nq=$fq\n"; # Set reads to singlereads in configuration file
		}
		close ($CON);

		system ("SOAPdenovo-Trans-127mer all -s SOAP_Conf.txt -K $k -o $asReads/SOAP/k$k -F -p 4"); # Run SOAPdenovo-Trans-127mer using assembly created config file and save output to AssembledReads/SOAP/k<k-mer_value>
		
		print $LOG "soap-q$q-k$k\t1\n"; # Write to log file
	}
}
my $soEndTime = time; # Set end time
my $soRunTime = $soEndTime - $soStartTime; # Calculate runtime
print $TIME "SOAP:\t$soRunTime\n" if $time;
# End: Assembly using SOAPdenovo-Trans assembly

## Assembly using SPAdes assembly
my $spStartTime = time; # Set start time

$kMin = 19; # Set minimum k-mer value
if ($readLength < 71) {
	$kMax = $readLength; # Set maximum k-mer value
} else {
	$kMax = 71; # Set maximum k-mer value
}

if (-d "$asReads/SPAdes/") { # Check if AssembledReads/SPAdes/ directory exists
	print "Found $asReads/SPAdes/\n";
} else {
	make_path ("$asReads/SPAdes/");  # Make AssembledReads/SPAdes/ if directory does not exist
	print "Made $asReads/SPAdes/\n";
}


my $finSp = $finished{"SPAdes"}; # Set finished job key value
$finSp = "none" unless $finSp;
for (my $k = $kMin; $k <= $kMax; $k += 4) { # Set k-mer value increment of 4
	if ($finSp =~ /q$q-k$k/) { # Skip if job already run
		print "Previously finished rnaspades-q$q-k$k, skipping this run\n";
	} else {
		if ($pe) {
			system("rnaspades.py -k $k -o $asReads/SPAdes/k$k --12 $fqI -m 500 -t 4"); # Run SPAdes on paired-end reads and save output to AssembledReads/SPAdes/k<k-mer_value>
		} else {
			system("rnaspades.py -k $k -o $asReads/SPAdes/k$k -s $fq -m 500 -t 4")  # Run SPAdes on single reads and save output to AssembledReads/SPAdes/k<k-mer_value>
		}
		
		print $LOG "SPAdes-q$q-k$k\t1\n"; # Write to log file
	}
}

my $spEndTime = time; # Set end time
my $spRunTime = $spEndTime - $spStartTime; # Calculate run time
print $TIME "SPAdes:\t$spRunTime\n" if $time;
# End: Assembly using SPAdes assembly

## Assembly using Trinity assembly
foreach my $q (@q) {
my $trStartTime = time; # Set start time
for (my $k = 15; $k < 32; $k += 4) { # Set k-mer increments of 4
		if ($status{"trinity-q$q-$k"}) { # Skip if job finished running
			print "Previously finished trinity-q$q, skipping this run\n";
		} else {
			if ($pe) {
				system("Trinity --KMER_SIZE $k --no_version_check --seqType fq --max_memory 250G --CPU 4 --output $asReads/Trinity/trinity-$k/ --full_cleanup --left $fq.1 --right $fq.2"); # Run Trinity on paired-end reads and save output to AssembledReads/Trinity/trinity-k<k-mer_value>
			} else {
				system("Trinity --KMER_SIZE $k --no_version_check --seqType fq --max_memory 250G --CPU 4 --output $asReads/Trinity/trinity-$k/ --full_cleanup --single $fq"); # Run Trinity on single reads and save output to AssembledReads/Trinity/trinity-k<k-mer_value>
			}
			print $LOG "trinity-q$q-$k\t1\n"; # Write to log file
		}
	}
	my $trEndTime = time; # Set end time
	my $trRunTime = $trEndTime - $trStartTime; # Calculate run timw
	print $TIME "Trinity:\t$trRunTime\n" if $time;
}
# End: Assembly using Trinity assembly
	
## Assembly using IDBA-Tran assembly	
my $idStartTime = time;
if ($status{"idba-q$q"}) { # Skip if job finished running
	print "Previously finished idba_trans, skipping this run\n";
} else {
	if ($pe) {
		system("idba_tran -o $asReads/idba/ -l $normReads/mergedI.fa --num_threads 4 --max_isoforms 50"); # Run IDBA-Tran on paired-end reads and save output to AssembledReads/idba/
	} else {
		system("idba_tran -o $asReads/idba/ -l $normReads/merged.fa --num_threads 4 --max_isoforms 50"); # Run IDBA-Tran on single reads and save output to AssembledReads/idba/
	}
	print $LOG "idba-q$q\t1\n"; # Write to log file
}
my $idEndTime = time; # Set end time
my $idRunTime = $idEndTime - $idStartTime; # Calculate run time
print $TIME "idbaTrans:\t$idRunTime\n" if $time;
# End: Assembly using IDBA-Tran assembly	

### Merging of transcripts from assembly of different k-mer lengths 
my ($IN, $OUT, $OUTAA, $TMP, $SPEC);

## Merge IDBA-Tran output transcripts
my $meStartTime = time; # Set start time
if ($status{"merged"}) { # Skip if transcripts already merged
	print "Previously finished merging transcrispts, skipping this run\n";
} else {

	open ($OUT, ">$asReads/mergedTranscripts.fa\n") or die "Can't open output: $!\n"; # Open AssembledReads/mergedTranscripts.fa for writing
	system("touch $asReads/mergedTranscripts.fa.faa"); # Create final contigs output file

	my $t;
	my %spec;

	my @idbaK = (20,30,40,50,60); # Set array for k-mer sizes
	foreach my $k (@idbaK) {
		if (-e "$asReads/idba/transcript-$k.fa") { # Check if AssembledReads/idba/transcript-<k-mer_value> exists
			open ($TMP, ">$asReads/mergingTmp.fa") or die "Can't open tmp file: $!\n"; # Open temp file for merging
			print "Found $asReads/idba/transcript-$k.fa, adding to merged transcripts\n";
			open ($IN, "$asReads/idba/transcript-$k.fa") or warn "Can't open $asReads/idba/transcripts-$k.fa"; # Open AssembledReads/idba/transcriopt-<k-mer_value>
			my $lt = 0;
			my ($id, $seq);
			while (my $line = readline($IN)) {
				if ($line =~ /^>/) { # Check if line is label
					$line =~ s/^>/>idbaK$k-q$q-/;
					if (length($line) > 50) {
						$line = substr($line, 0, 50);
						$line .= "\n";
					}
					$spec{"idba"}{$id} = $seq if $seq; # Add to general hash holding merged transcripts for all assembly methods
					$id = $line;
					$seq = "";
					print $OUT $line;
					print $TMP $line;
					$lt++;
					$t++;
				} else {
					$seq .= $line; # Set sequence to read line
					print $OUT $line; # Write to final merge transcripts output file
					print $TMP $line; # Write to temp merge transcripts output file
				}
			}
			$spec{"idba"}{$id} = $seq if $seq; # Add to general hash holding merged transcripts for all assembly methods
			close($IN);
			close($TMP);
			
		 	chdir($asReads);
			system("ORFfinder -in mergingTmp.fa -out mergingTmpORFs.fa"); # Run ORFfinder to find open reading frames in the temp merged transcripts
			 	
			getLongestORF();
				
		 	system("cat mergingTmpORFs.fa.faa >> mergedTranscripts.fa.faa"); # Append contigs to final contigs output file
		 	chdir($runDir);
		 	
			print "Successfully added $lt contigs\n";
		}
	}
	# End: Merge IDBA-Tran output transcripts
	
	## Merge Trinity output transcripts
	for (my $k = $kMin; $k <= $kMax; $k += 4) {
		if ($k <= 31) {
			if (-e "$asReads/Trinity/trinity-$k.Trinity.fasta") { # Check if AssembledReads/Trinity/trinity-<k-mer_value> exists
				open ($TMP, ">$asReads/mergingTmp.fa") or die "Can't open tmp file: $!\n";  # Open temp file for merging
				print "Found $asReads/Trinity/trinity-$k.Trinity.fasta, adding to merged transcripts\n";
				open ($IN, "$asReads/Trinity/trinity-$k.Trinity.fasta") or warn "Can't open $asReads/Trinity/trinity-$k.Trinity.fasta: $!\n"; # Open AssembledReads/Trinity/trinity-<k-mer_value>
				my $lt = 0;
				my ($id, $seq);
				while (my $line = readline($IN)) { # Check if line is label
					if ($line =~ /^>/) {
						$line =~ s/^>/>Trinity-q$q-$k-/;
						if (length($line) > 50) {
							$line = substr($line, 0, 50);
							$line .= "\n";
						}
						$spec{"Trinity"}{$id} = $seq if $seq; # Add to general hash holding merged transcripts for all assembly methods
						$id = $line;
						$seq = "";
						print $OUT $line;
						print $TMP $line;
						$lt++;
						$t++;
					} else {
						$seq .= $line;  # Set sequence to read line
						print $OUT $line; # Write to final merge transcripts output file 
						print $TMP $line; # Write to temp merge transcripts output file
					}
				}
				close $IN;
				close($TMP);
				$spec{"Trinity"}{$id} = $seq if $seq; # Add to general hash holding merged transcripts for all assembly methods
							
				chdir($asReads);
			 	system("ORFfinder -in mergingTmp.fa -out mergingTmpORFs.fa"); # Run ORFfinder to find open reading frames in the temp merged transcripts
			 	
			 	getLongestORF();
				
		 		system("cat mergingTmpORFs.fa.faa >> mergedTranscripts.fa.faa"); # Append contigs to final contigs output file
				chdir($runDir);
				
				print "Successfully added $lt contigs\n";
			}
		}
		# End: Merge Trinity output transcripts

		## Merge SOAPdenovo output transcripts
		if (-e "$asReads/SOAP/k$k.scafSeq") { # Check if AssembledReads/SOAP/k-<k-mer_value> exists
			open ($TMP, ">$asReads/mergingTmp.fa") or die "Can't open tmp file: $!\n"; # Open temp file for merging
			print "Found $asReads/SOAP/k$k.scafSeq, adding to merged transcripts\n";
			open ($IN, "$asReads/SOAP/k$k.scafSeq") or warn "Can't open $asReads/SOAP/k$k.scafSeq: $!\n"; # Open AssembledReads/SOAP/k-<k-mer_value>
			my $lt = 0;
			my ($id, $seq);
			while (my $line = readline($IN)) {
				if ($line =~ /^>/) { # Check if line is label
					$line =~ s/^>/>SOAPdenovoK$k-q$q-/;
					if (length($line) > 50) {
						$line = substr($line, 0, 50);
						$line .= "\n";
					}
					$spec{"SOAPdenovo"}{$id} = $seq if $seq; # Add to general hash holding merged transcripts for all assembly methods
					$id = $line;
					$seq = "";
					print $OUT $line;
					print $TMP $line;
					$lt++;
					$t++;
				} else {
					$seq .= $line; # Set sequence to read line
					print $OUT $line; # Write to final merge transcripts output file
					print $TMP $line; # Write to temp merge transcripts output file
				}
			}
			$spec{"SOAPdenovo"}{$id} = $seq if $seq; # Add to general hash holding merged transcripts for all assembly methods
			close $IN;
			close($TMP);
			
		 	chdir($asReads);
		 	system("ORFfinder -in mergingTmp.fa -out mergingTmpORFs.fa"); # Run ORFfinder to find open reading frames in the temp merged transcripts
		 	
			getLongestORF();
			
		 	system("cat mergingTmpORFs.fa.faa >> mergedTranscripts.fa.faa"); # Append contigs to final contigs output file
		 	chdir($runDir);
		 	
			print "Successfully added $lt contigs\n";
		}
		# End: Merge SOAPdenovo output transcripts
	
		## Merge SPAdes output transcripts
		if (-e "$asReads/SPAdes/k$k/transcripts.fasta") { # Check if AssembledReads/SPAdes/k<k-mer_value> exists
			open ($TMP, ">$asReads/mergingTmp.fa") or die "Can't open tmp file: $!\n"; # Open temp file for merging
			print "Found $asReads/SPAdes/k$k/transcripts.fasta, adding to merged transcripts\n";
			open ($IN, "$asReads/SPAdes/k$k/transcripts.fasta") or warn "Can't open $asReads/SPAdes/k$k/transcripts.fasta: $!\n"; # Open AssembledReads/SPAdes/k<k-mer_value>
			my $lt = 0;
			my ($id, $seq);
			while (my $line = readline($IN)) {
				if ($line =~ /^>/) { # Check if line is label
					$line =~ s/^>/>SPAdesK$k-q$q-/;
					if (length($line) > 50) {
						$line = substr($line, 0, 50);
						$line .= "\n";
					}
					$spec{"SPAdes"}{$id} = $seq if $seq; # Add to general hash holding merged transcripts for all assembly methods
					$id = $line;
					$seq = "";
					print $OUT $line;
					print $TMP $line;
					$lt++;
					$t++;
				} else {
					$seq .= $line; # Set sequence to read line
					print $OUT $line; # Write to final merge transcripts output file
					print $TMP $line; # Write to temp merge transcripts output file
				}
			}
			$spec{"SPAdes"}{$id} = $seq if $seq; # Add to general hash holding merged transcripts for all assembly methods
			close $IN;
			close($TMP);
			
		 	chdir($asReads);
			system("ORFfinder -in mergingTmp.fa -out mergingTmpORFs.fa"); # Run ORFfinder to find open reading frames in the temp merged transcripts
			
			getLongestORF();
				
		 	system("cat mergingTmpORFs.fa.faa >> mergedTranscripts.fa.faa"); # Append contigs to final contigs output file
		 	chdir($runDir);
		 	
			print "Successfully added $lt contigs\n";
		}
		# End: Merge SPAdes output transcripts
	}

	print "Printing assembler specific results:\n";
	foreach my $as (sort {$a cmp $b} keys %spec) {
		print "Printing merged results for $as:\n";
		open ($SPEC, ">$asReads/$as.fasta") or warn "Can't open $as output: $!\n";
		foreach my $g (sort {$a cmp $b} keys %{$spec{$as}}) {
			print $SPEC "$g\n$spec{$as}{$g}\n"
		}
	}

	print "Added a total of $t contigs to $asReads/mergedTranscripts.fa\n";
	print $LOG "merged\t1\n"; # Write to log file
	close $OUT;
}

my $meEndTime = time; # Set end time
my $meRunTime = $meEndTime - $meStartTime; # Calculate run time
print $TIME "Merging:\t$meRunTime\n" if $time;# 
# 
# if ($status{"gm-q$q"}) {
# 	print "Previously finished GeneMarkS, skipping this run\n";
# } else {
# 	print "Running GeneMarkS to predict amino acid sequences:\n";
# 	chdir($asReads);
# 	system("gmsn.pl --faa --euk mergedTranscripts.fa");
# 	chdir($runDir);
# 	print $LOG "gm-q$q\t1\n";
# }

# End: Merging of transcripts from assembly of different k-mer lengths 

### Analysis of merged contigs - Finding of consensus sequences
my ($FASTA, $PROT, $POUT, $LOUT, $SOUT, $SOUT4, $LOUT4,  $POUT4);

# my $asReads = "$asReads/";

open ($FASTA, "$asReads/mergedTranscripts.fa") or die "Can't open $asReads/mergedTranscripts.fa file: $!\n"; # Open merged transcripts for reading
open ($PROT, "$asReads/mergedTranscripts.fa.faa") or die "Can't open $asReads/mergedTranscripts.fa.faa file: $!\n"; # Open ORF found merged transcripts for reading
open ($SOUT, ">$asReads/ConSemble3shortest.fasta") or die "Can't open nt output: $!\n"; # Open consemble3+ file to hold short consensus sequences for writing
open ($LOUT, ">$asReads/ConSemble3longest.fasta") or die "Can't open nt output: $!\n"; # Open consemble3+ file to hold long consensus sequences for writing
open ($POUT, ">$asReads/ConSemble3.fa.faa") or die "Can't open aa output: $!\n"; # Open consemble3+ file to hold consensus sequences
open ($SOUT4, ">$asReads/ConSemble4shortest.fasta") or die "Can't open nt output: $!\n"; # Open consemble4+ file to hold shortest consensus sequences for writing
open ($LOUT4, ">$asReads/ConSemble4longest.fasta") or die "Can't open nt output: $!\n"; # Open consemble4+ file to hold longest consensus sequences for writing
open ($POUT4, ">$asReads/ConSemble4.fa.faa") or die "Can't open aa output: $!\n"; # Open consemble4+ file to hold consensus sequences

my %fasta;
my ($seq, $id);
print "Finding consensus between assemblers:\n";
while (my $line = readline($FASTA)) {
	chomp $line;
	next unless $line;
	
	if ($line =~ /^>/) { # Check if read line is label
		$fasta{$id} = $seq if $seq; # If $seq is set, mostly from previous iteration, set $fasta to that sequence where id = $id
		$seq = ""; # Set sequence to empty string
		$line =~ s/^>//; # Remove ">" from read label
		my @info = split /\s+/, $line; # Split read label
		$id = $info[0]; # Set id to first element in split array of label, to be used in the next iteration
	} else {
		$seq .= $line; # Concatenate read line to sequence to create one long string
	}
}
$fasta{$id} = $seq if $seq; # If $seq is set, from last iteration in while loop, set $fasta to that sequence where id = $id

my (%prot, %match);
($seq, $id) = ("", "");
while (my $line = readline($PROT)) {
	chomp $line;
	next unless $line;
	
	if ($line =~ /^>/) { # Check if read line is label
		$prot{$id} = $seq if $seq; # If $seq is set, mostly from previous iteration, set $prot to that sequence where id = $id
		$line =~ s/^>//; # Remove ">" from read label
		my @info = split /\s+/, $line; # Split read label
		$info[0] =~ s/^>//; # Remove ">" from read label
		if ($seq) { # Check if $seq has value
			if ($match{$seq}) { # Check if $match where seq = $seq has value
				$match{$seq} .= ",$id" unless $match{$seq} =~ /$id/; # Concantenate id, set in the previous iteration, to $match where seq = $seq, unless $id is excatly the same
			} else { 
				$match{$seq} = $id # Set $match where seq = $seq to $id, set in previous iteration
			}
		}
		$id = $info[0]; # Set id to first element in split array of label, to be used in the next iteration
		$seq = "";
	} else {
		$seq .= $line; # Concatenate read line to sequence to create one long string
	}
}

$prot{$id} = $seq if $seq; # If $seq is set, from last iteration in while loop, set $prot to that sequence where id = $id

if ($seq) { # Check if $seq has value
	if ($match{$seq} and $id) { # Check if $match where seq = $seq and $id have values
		$match{$seq} .= ",$id" unless $match{$seq} =~ /$id/; # Concantenate id, set in the last iteration of while loop, to $match where seq = $seq
	} else { 
		$match{$seq} = $id; # Set $match where seq = $seq to $id, set in last iteration of while loop
	}
}

## Counting and matching consensus sequences
my (%fm_count, %counted);
my $t = 0;
my $c = 1;
foreach my $hit (sort {$a cmp $b} keys %match) { # Loop through the $match hash holding the sequences and ids of those sequences from merged transcripts
	next if $counted{$hit}; # Skip if sequence has been counted
	$counted{$hit} = 1; # Set $counted attribute for sequence in $match to 1
	
	my $hits = $match{$hit}; # Get values of ids belonging to the sequence being processed in current iteration
	my @list = split /,/, $hits; # Split ids and store in @list array
	
	my %as;
	foreach my $h (@list) {	# Loop through ids array
		my @ask = split /-/, $h; # Split id being processed in current iteration and store in @ask array
		my @as = split /K/, $ask[0]; # Split first item from @ask array - This would contain the assembly method from which the id was extracted
		$as{$as[0]} = 1; # Set $as attribute for extracted assembly method to 1
	}
	my $as_list = join ",",  sort({$a cmp $b} keys %as); # Join extratcted assembly methods from ids in $hits into a single string $as_list
	$fm_count{$as_list}++; # Increment counter $fm_count($as_list) - Tracking the number of times assembly method overlap
	
	my $denovoCount = 0; # Variable to keep track of number of assembly methods
	foreach my $a (sort {$a cmp $b} keys %as) { # Loop through all extracted assembly methods
		$denovoCount++ if ($a eq "SOAPdenovo" or $a eq "idba" or $a eq "SPAdes" or $a eq "Trinity"); # Increment $denovoCount if assembly is one of the listed
	}
	if ($denovoCount >= 3) { # Check if $denovoCount is 3 or greater
		my $longestLength = 0; # Set attibute to store length of longest contig
		my $shortestLength = 999999; # Set attribute to store length of shortest contig
		my $longest = ""; # Set attibute to store longest contig
		my $shortest = ""; # Set attribute to store shortest contig
		foreach my $h (@list) { # Loop throught @list array of ids belonging to $hit sequence being processed in current iteration
			next if $counted{$h}; # Skip if id is counted
			$counted{$h} = 1; # Set $counted attribute for id in @list to 1
			my $clen = length($fasta{$h}); # Get length of sequnce belonging to id $h in $fasta
			if (length($fasta{$h}) > $longestLength) {
				$longestLength = length($fasta{$h}); # Set $longestLength to length of $clen
				$longest = $h; # Set @longest attribute to $h
			}
			if (length($fasta{$h}) < $shortestLength) {
				$shortestLength = length($fasta{$h}); # Set $shortestLength to length of $clen
				$shortest = $h; # Set @shortest attribute to $h
			}
		}
		next unless $longest; # Go to next sequence $hit if not longest
		$t++;
		print $POUT ">contig_$c\n$prot{$longest}\n"; # Wite to consemble3+ file
		print $POUT4 ">contig_$c\n$prot{$longest}\n" if $denovoCount >= 4; # Wite to consemble4+ file
		if ($longest) {
			print $LOUT ">contig_$c\n$fasta{$longest}\n"; # Wite to longest consemble3+ file
			print $LOUT4 ">contig_$c\n$fasta{$longest}\n" if $denovoCount >= 4; # Wite to longest consemble4+ file
		}
		if ($shortest) {
			print $SOUT ">contig_$c\n$fasta{$shortest}\n"; # Wite to shortest consemble3+ file
			print $SOUT4 ">contig_$c\n$fasta{$shortest}\n" if $denovoCount >= 4; # Wite to longest consemble4+ file
		}
		$c++;
	}
}
# End: Counting and matching consensus sequences

#Uncomment this block to see the number of contigs shared between each assembler
#######################
# foreach my $al (sort {$a cmp $b} keys %fm_count) {
# 	print "$al\t$fm_count{$al}\n";
# }
#######################

## Count contigs in final assembly
print "Number of contigs in the final assembly: $t\n";

## Subroutine to get logest reading frame
sub getLongestORF {
	my %tmpGenes;
	my ($tmpId, $tmpSeq);
	my ($TMP, $TMPOUT);
	open ($TMP, "mergingTmpORFs.fa") or die "Can't open mergingTmpORFs.fa\n"; # Open mergingTmpORFs.fa for reading
	open ($TMPOUT, ">mergingTmpORFs.fa.faa") or die "Can't open mergingTmpORFs.fa.faa\n"; # Open mergingTmpORFs.fa.faa for writing

	while (my $line = readline($TMP)) {
		chomp $line;
		next unless $line;
		if ($line =~ /^>/) { # CHeck if line is label
			if ($tmpSeq) { # Check if sequence is set
				# Set sequence for $tmpGenes{$tmpId} to empty string unless $tmpGenes{$tmpId} is set, usually in previous iteration
				$tmpGenes{$tmpId} = "" unless $tmpGenes{$tmpId};
				# Set sequence for $tmpGenes{$tmpId} to #tempseq to longest sequence, compared with the one set in previous iteration 
				$tmpGenes{$tmpId} = $tmpSeq if (length($tmpSeq) > length($tmpGenes{$tmpId}));
			}
			my @info = split /ORF\d+_/, $line; # Split read line and store in @info array
			my @gn = split /:/, $info[1]; # Split the first item in @info arrayand stor in @gn array
			$tmpId = $gn[0]; # Set $tempId to first item of @gn array
			$tmpSeq = ""; # Set $tmpSeq to empty string
		} else {
			$tmpSeq .= $line; # Concantenate read line to $tmpSeq
		}
	}
	# Set sequence for $tmpGenes{$tmpId} to #tempseq to longest sequence, compared with the one set in last iteration of while loop
	$tmpGenes{$tmpId} = $tmpSeq if (length($tmpSeq) > length($tmpGenes{$tmpId}));

	foreach my $g (sort {$a cmp $b} keys %tmpGenes) {
		print $TMPOUT ">$g\n$tmpGenes{$g}\n"; # Write to output file
	}

	close ($TMP);
	close ($TMPOUT);
}
## End: Subroutine to get logest reading frame
