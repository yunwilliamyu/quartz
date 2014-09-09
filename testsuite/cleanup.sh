#!/bin/bash
for file in example?.fastq
do
	echo Cleaning up files generated from $file
	rm -rf	$file.\
			$file.filtered_S \
			$file.filtered_S.qual \
			$file.filtered_S.qual.bz2
done
echo "Removing generated dictionary files"
rm -rf dict.db dict_flt.db dict.bin dict.bin.sorted dict.bin.sorted.swapped
