Copyright (c) 2014 Y. William Yu. All rights reserved.

This package implements the Quartz algorithm described in:
Y. William Yu, Deniz Yorukoglu, Jian Peng, and Bonnie Berger.
"Quality Score Compression Improves Genotyping"
Manuscript submitted for publication.

http://quartz.csail.mit.edu/

-----------------------------
Package contents:

Following programs take input in FASTQ format:
misra_gries_dict:  builds a dictionary of common k-mers from a corpus by
                   approximate counting
    dict_txt2bin:  converts text dictionary to binary format
  sort_dict_file:  sorts a binary dictionary and outputs it
                   also swaps the high and low order bits, sorts it, and outputs
				   the swapped dictionary as well
          quartz:  uses dictionary to smooth quality values for high confidence
                   calls as measured by k-mer Hamming distance to a default Q.

-----------------------------
Requirements:
	64 GiB RAM

Compilation:
	make all

Preprocessing:
	./misra_gries_dict MINCOUNT -1 dict.db *.fastq
	./decrement_misra_gries.py dict.db > dict_flt.db
	./dict_txt2bin dict_flt.db dict.bin
	./sort_dict_file dict.bin dict.bin.sorted dict.bin.sorted.swapped

Quickstart:
	./quartz dict.bin.sorted "Q" 8 *.fastq

We also provide a testsuite/ directory with example FASTQ files to
demonstrate operation. To run the commented example scripts:
	cd testsuite/
	./generate_dictionary.sh
	./run_quartz.sh

Note: dictionary generation can be very expensive. The authors have already
generated a high quality human genome dictionary, available for download at
	http://yunwilliamyu.net/quartz/dec200.bin.sorted.gz
	MD5: b8f8409d7bd7fd2beea1ee9d5b68b6d1

	http://yunwilliamyu.net/quartz/dec200.bin.sorted.swapped.gz
	MD5: e7f0dc501ee05dba12cdd2f5ae8540ad
Naturally, dict.bin.sorted should be replaced by dec200.bin.sorted when using
this dictionary. Note that when generating from a large number of FASTQ files,
it is sometimes necessary to generate the dictionary in chunks by using the 2nd
command-line argument of misra_gries_dict, and then combining the resulting files.

Quartz performs best on recent NGS Illumina FASTQ reads. As an example, the
following files from 1000 Genomes Project, for NA12878, were used for the
manuscript:
	ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/NA12878/sequence_read/SRR622461_1.filt.fastq.gz
	ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/NA12878/sequence_read/SRR622461_2.filt.fastq.gz
	ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/data/NA12878/sequence_read/SRR622461.filt.fastq.gz
	
-----------------------------
The main program is quartz:

quartz:
Discards likely non-SNP quality scores for known reads.
Usage: ./quartz dictionary_file [QUAL] [num_threads] input_file(s)
	Input is assumed to be a standard FASTQ file.

	Output will be modified FASTQ files named [input_file].filtered_[QUAL],
	identical to the original except with changed quality scores.

	The quality of nearly all bases corresponding to a 32-mer listed in the 1st
	column of [dictionary_file] will be set to [QUAL]. Correspondance shall be
	defined as a Hamming distance of less than or equal to 1 in the 32-mer.
	Note that a 32-mer can correspond to multiple 32-mers in the database.
	The exceptions to the setting of quality values above will be all bases
	that are different to any of the corresponding 32-mers in the database,
	which will retain their original quality value data.

	Example w/ 8-mers: suppose "AAAAAAAA" and "TAAAAAAT" are in the database,
	and we have the 8-mer "TAAAAAAA" in the input_file with original qualities
	"ABCDEFGH". Then the new quality values will be "A~~~~~~H", where QUAL="~"

	Note that for reads longer than 32 bases, quartz will only check 32-mers
	for matches at most once every 16 bases, and once more at the end.

	Quality values of 2 (i.e. "#") are ignored and not modified.

	NB: For "~" in particular, be careful about shell expansion, so you will
	wnat to use single quotes '~'.

	NB: quartz requires that dictionary_file.swapped be the swapped version
	of the dictionary in the same directory as dictionary_file.

	NB: [num_threads] specifies the number of threads that should be used.
	On our machines, disk I/O becomes a bottleneck at around 8 threads. Using
	more threads than your disk I/O can support is not recommended and may have
	a detrimental effect on speed.

Preprocessor is centered around generating the dictionary:

misra_gries_dict:
Usage: ./misra_gries_dict MINCOUNT MASK output_file input_file(s)
	Approximates the number of times frequent 32-mers appear, outputs two column
	format, with the first column specifying the 32-mer and the second column
	the number of times it appears in the corpus + decrement. MINCOUNT specifies
	the minimum number of appearances to be output in the database.

	MASK \in [0,255] specifies the middle 4 base pairs. Set to -1 to get all.

	Note also that the first line specifies the decrement. To get actual
	approximate Misra-Gries counts, you'll need to subtract the decrement from
	the second column. A simple Python script is provided for this purpose,
	decrement_misra_gries.py, which takes as first argument the output of this
	program.

dict_txt2bin:
Converts text dictionary to binary
Usage: ./dict_txt2bin text_dict.db binary_dict.db
	Converts the 32-mer in the first column of text_dict.db to a 64-bit integer.
	binary_dict.db has the total number of 32-mers in the dictionary as the
	first 8 bytes, and each of the rest of the entries following in 8 bytes each.

sort_dict_file:
Usage: ./sort_dict_file input_dict.bin output_dict.bin output_dict_swapped.bin
	Sorts a dictionary file, outputs it, swaps high and low bits, sorts it, 
	and outputs it.

	Note that for use with quartz, the name of the swapped dictionary *must*
	be the name of the first sorted dictionary + ".swapped".

