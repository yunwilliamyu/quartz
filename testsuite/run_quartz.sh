#!/bin/bash

if [ ! -f ../quartz ];
then
	echo "quartz not found. Please run 'make all' in the root directory and rerun this script."
	exit
fi

echo In this example, we compress the corpus we used to generate the dictionary:
echo "(Quartz here is run in low memory mode; to run in high-memory mode, replace 0 with 1 in the following line)"
echo "../quartz dict.bin.sorted 'S' 8 0 *.fastq"
../quartz dict.bin.sorted 'S' 8 0 *.fastq
echo

echo Compression is now done, but let\'s see how well we did using BZIP2:
echo "We'll splice out the quality scores and pipe it through BZIP2"
echo and then compute the bits per quality value needed.
for file in *.filtered_S
do
	cat $file  | awk '{if (NR%4==0) print $0}' > $file.qual
	orig_size=`wc -c < $file.qual`
	orig_lines=`wc -l < $file.qual`
	orig_size=`echo "$orig_size - $orig_lines" | bc`
	bzip2 -f $file.qual
	new_size=`wc -c < $file.qual.bz2`
	echo -e $file:'\t' `echo "scale=4; 1/( $orig_size / ( $new_size * 8)) " | bc` bits / quality score
done
echo

echo "All done!"
echo
echo "If you want to clean up all the files generated in this example, just run"
echo "./cleanup.sh"
