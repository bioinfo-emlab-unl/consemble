(consemble.pl: last updated on January 8, 2019)
(README: last updated on March 30, 2021)

####################################################
TESTED:
Linux
#####################################################


####################################################
Example Usages:

Help:
perl consemble3+g.pl --help

Single-end read assembly
perl --assemblyName testData --single testData_1.fastq --ref reference.fasta

Paired-end read assembly (Produces the output shown below)
perl --assemblyName testData --reads1 testData_1.fastq --reads2 testData_2.fastq --ref reference.fasta
####################################################


####################################################
PREPARATION: (DO THIS BEFORE YOU RUN THE PIPELINE)

For the assembly pipeline to run, the following programs MUST be in the PATH:

	erne-filter v2.0 (http://erne.sourceforge.net/index.php)
	khmer v2.0 (https://github.com/dib-lab/khmer)
	Tophat2 v2.1.1 (https://ccb.jhu.edu/software/tophat/index.shtml)
	Bayesembler v1.2.0 (https://github.com/bioinformatics-centre/bayesembler)
	Cufflinks v2.2.1 (https://github.com/cole-trapnell-lab/cufflinks)
	Scallop v0.10.0 (https://github.com/Kingsford-Group/scallop)
	StringTie v2.0.0 (https://github.com/skovaka/stringtie2)
	ORFfinder (https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/)

Later versions of these programs should work, but have not been tested.

Create a directory and copy the input read sequence file(s) in FASTQ format in this directory.  All of the reads for the assembly should be combined into 1 single-end read file, or 2 paired-end read files.

ALL output files will be saved in the current working directory the pipeline is called from unless a directory is specified with --outputDir option. This directory and all of the subdirectories will be automatically created by the pipeline.

In the above "Paired-end read assembly" example, the following directories are generated:
	./FilteredReads/testData/q20/	Location of quality filtered reads before normalization
	./NormalizedReads/testData/	Location of normalized reads used for assemblies
	./AssembledReads/testData/	Location of all output files produced for the assemblies and final ConSemble results
If the directories cannot be created, the execution of the pipeline will fail.
####################################################


####################################################
INPUTS and OPTIONS:
	--single	Specifies the file including single-end reads
	--reads1	Specifies the first file for paired-end reads
	--reads2	Specifies second file for paired-end reads
	--assemblyName	Specifies the assembly name (e.g. testData) to organize output files
	--ref		Location for the reference Fasta file for read mapping
	--q		Changes q-value threshold for reads filtering by erne-filter (default 20)
	--maxMemory	Maximum memory for pipeline in GB (default 64)
	--maxThreads	Maxmimum number of threads (default 8)
	--outputDir	Changes the base directory for output files (default is current working directory)
	--force		Ignores the log file and runs entire pipeline.  See the description of log.txt below.
	--noNorm	Skips the digital normalization of filtered reads
	--normCov	Changes the max kmer coverage for digital normalization. Lowering coverage speeds up assembly time, but increases risk of artifacts (default 50)
	--time		Records the amount of time each stage of the pipeline takes
	--version	Prints pipeline version without executing the pipeline
	--h or --help	Prints these options without executing the pipeline
	

####################################################
OUTPUT FILES: all the files mentioned below will be written to the output directory chosen by the "--outputDir" option and the subdirectory "testData" is chosen by the "--assemblyName" option (see the example above).

SUBDIRECTORY - "FilteredReads/testData/q20"
1. merged_1.fastq
--------------------------
Reads that pass the quality filtering from the first read file.  Output by erne-filter.

2. merged_2.fastq
--------------------------
Reads that pass the quality filtering from the second read file.  Output by erne-filter.

3. mergedI.fastq
--------------------------
Interleaved reads from the filtered merged_1.fastq and merged_2.fastq, input file for read normalization.  Output by interleave-reads.py from khmer.

4. merged_unpaired.fastq
--------------------------
Reads orphaned by quality filtering.  Output by interleave-reads.py from khmer.


SUBDIRECTORY - "NormalizedReads/testData"
1. merged.fa.1
--------------------------
Fasta version of normalized reads from the first read file.  Output by split-paired-reads.py from khmer.

2. merged.fa.2
--------------------------
Fasta version of normalized reads from the second read file.  Output by split-paired-reads.py from khmer.

3. merged.fq.1
--------------------------
Normalized reads from the first read file.  Output by split-paired-reads.py from khmer.

4. merged.fq.2
--------------------------
Normalized reads from the second read file. Output by split-paired-reads.py from khmer.

5. mergedI.fq
--------------------------
Interleaved normalized reads. Output by normalize-by-median.py from khmer.

6. mergedI.fa
--------------------------
Interleaved reads in fasta format. Output by fastq-to-fasta.py from khmer.

7. mergedI.fq.se
--------------------------
Reads orphaned by the read normalization.  Output by extract-paired-reads.py from khmer.

SUBDIRECTORY - "AssembledReads/testData"
1. consemble3+g.faa
--------------------------
Protein sequences produced by the ConSemble pipeline (protein sequences produced by at least 3 of the 4 assemblers)

2. consemble3+g.fa
--------------------------
Nucleotide sequences produced by ConSemble pipeline.  If multiple contigs produce a protein sequence in consensus.aa, the shortest contig is chosen.

3. accepted_hits.bam
--------------------------
Read mapping produced by Tophat2

4. bayesembler-assembly.gtf
--------------------------
GTF file from the assembly for Bayesembler

5. bayesembler-assembly.fa
--------------------------
Fasta file based on bayesembler-assmebly.gtf and the reference genome from gtf_to_fasta (part of Tophat2)

6. bayesembler-assembly-ORFs.fa
--------------------------
Protein translations of from all open reading frames (ORFs) in bayesembler-assembly.fa from ORFfinder

7. bayesembler-assembly-ORFs.fa.faa
--------------------------
The longest ORF for each contig in bayesembler-assembly-ORFs.fa

8. cufflinks/transcripts.gtf
--------------------------
GTF file from the assembly Cufflinks

9. cufflinks/transcripts.fa
--------------------------
Fasta file based on cufflinks/transcripts.gtf and the reference genome from gtf_to_fasta (part of Tophat2)

10. cufflinks/transcripts-ORFs.fa
--------------------------
Protein translations of from all open reading frames (ORFs) in cufflinks/transcripts.fa from ORFfinder

11. cufflinks/transcripts-ORFs.fa.faa
--------------------------
The longest ORF for each contig in cufflinks/transcripts-ORFs.fa

12. scallop.gtf
--------------------------
GTF file for assembly for Scallop

13. scallop.fa
--------------------------
Fasta file based on scallop.gtf and the reference genome from gtf_to_fasta (part of Tophat2)

14. scallop-ORFs.fa
--------------------------
Protein translations of from all open reading frames (ORFs) in scallop.fa from ORFfinder

15. scallop-ORFs.fa.faa
--------------------------
The longest ORF for each contig in scallop-ORFs.fa

16. stringtie.gtf
--------------------------
GTF file for assembly for StringTie2

17. stringtie.fa
--------------------------
Fasta file based on stringtie.gtf and the reference genome from gtf_to_fasta (part of Tophat2)

18. stringtie-ORFs.fa
--------------------------
Protein translations of from all open reading frames (ORFs) in stringtie.fa from ORFfinder

19. stringtie-ORFs.fa.faa
--------------------------
The longest ORF for each contig in stringtie-ORFs.fa

20. log.txt
--------------------------
Log of previously completed steps for assembly.  If the assembly is interrupted, any steps listed in this file will be skipped when the assembly is restarted.  If the --force flag is used, the log file is ignored and the full pipeline is run reglardless of previously complete steps

Depending on the environment the pipeline is run in, additional files may be generated.

####################################################
[NOTES]
1.  Sequences produced by only 1 or 2 assemblers are not kept in ConSemble3+g.faa or ConSemble3+g.fa.

####################################################
