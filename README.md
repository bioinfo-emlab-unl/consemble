# ConSemble
A consensus-based ensemble approach to improve transcriptome assembly

(consemble.pl: last updated on January 8, 2019)
(README: last updated on March 23, 2021)

####################################################
TESTED:
Linux
#####################################################


####################################################
Example Usages:

Help:
perl consemble3+d.pl --help

Single-end read assembly
perl --assemblyName testData --single testData_1.fastq

Paired-end read assembly (Produces the output shown below)
perl --assemblyName testData --reads1 testData_1.fastq --reads2 testData_2.fastq
####################################################


####################################################
PREPARATION: (DO THIS BEFORE YOU RUN THE PIPELINE)

For the assembly pipeline to run, the following programs MUST be in the PATH:

	erne-filter v2.0 (http://erne.sourceforge.net/index.php)
	khmer v2.0 (https://github.com/dib-lab/khmer)
	IDBA-Tran v1.1.1 (http://i.cs.hku.hk/~alse/hkubrg/projects/idba_tran/)
	SOAPdenovo-Trans v1.0.3 (http://soap.genomics.org.cn/SOAPdenovo-Trans.html)
	rnaSPAdes v3.10.0 (http://bioinf.spbau.ru/en/rnaspades)
	Trinity v2.4.0 (https://github.com/trinityrnaseq/trinityrnaseq/wiki)
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
1. consensus.aa
--------------------------
Protein sequences produced by the ConSemble pipeline (protein sequences produced by at least 3 of the 4 assemblers)

2. consensus.fasta
--------------------------
Nucleotide sequences produced by ConSemble pipeline.  If multiple contigs produce a protein sequence in consensus.aa, the shortest contig is chosen.

3. idba.fasta
--------------------------
Merged nucleotide sequences from all IDBA-Tran assemblies

4. SOAPdenovo.fasta
--------------------------
Merged nucleotide sequences from all SOAPdenovo-trans assemblies

5. SPAdes.fasta
--------------------------
Merged nucleotide sequences from all rnaSPAdes assemblies

6. Trinity.fasta
--------------------------
Merged nucleotide sequences from all Trinity assemblies

7. mergedTranscripts.fa
--------------------------
Merged nucleotide sequences from all assemblies.  The sequence names are changed to identify which assembly produced the sequence.

8. mergedTranscripts.fa.faa
--------------------------
Merged protein sequences produced by ORFfinder from all assemblies.  The sequence names are changed to identify which assembly produced the sequence.  Used as the input for ConSemble.

9. idba/
--------------------------
Contains all of the files produced by the IDBA-Tran assemblies

10. SOAP/
--------------------------
Contains all of the files produced by the SOAPdenovo-Trans assemblies

11. SPAdes/
--------------------------
Contains all of the files produced by the rnaSPAdes assemblies

12. Trinity/
--------------------------
Contains all of the files produced by the Trinity assemblies

13. log.txt
--------------------------
Log of previously completed steps for assembly.  If the assembly is interrupted, any steps listed in this file will be skipped when the assembly is restarted.  If the --force flag is used, the log file is ignored and the entire pipeline is run.

14. mergingTmp.fa and mergingTmp.fa.faa
--------------------------
Temporary files for ORFfinder to translate individual assemblies

Depending on the environment the pipeline is run in, additional files may be generated.

####################################################
[NOTES]
1.  Sequences produced by only 1 or 2 assemblers are not kept in consensus.aa or consensus.fasta, but can still be found in the merged transcripts files.
2.  Merged assemblies in the AssembledReads subdirectory are not unique at the nucleotide or amino acid level, and may contain a large number of duplicates.
3.  All khmer lengths used in the assemblies are hard-coded as followed:
	IDBA-Tran: 20-60 every 10
	SOAPdenovo-Trans: 15-127 every 4 until the max read length (automatically determined by the script)
	rnaSPAdes: 19-71 every 4 until the max read length (automatically determined by the script)
	Trinity: 15-31 every 4.

####################################################
