#!/bin/bash

if [ ! -f ../misra_gries_dict ];
then
	echo "misra_gries_dict not found. Please run 'make all' in the root directory and rerun this script."
	exit
fi

echo Generate our dictionary with read multiplicity r=2:
echo "../misra_gries_dict 2 -1 dict.db *.fastq"
../misra_gries_dict 2 -1 dict.db *.fastq
echo

echo Stripping away the decrement counter line:
echo "../decrement_misra_gries.py dict.db > dict_flt.db"
../decrement_misra_gries.py dict.db > dict_flt.db
echo

echo Converting to binary format:
echo "../dict_txt2bin dict_flt.db dict.bin"
../dict_txt2bin dict_flt.db dict.bin
echo

echo Sorting and swapping dictionary:
echo "../sort_dict_file dict.bin dict.bin.sorted dict.bin.sorted.swapped"
../sort_dict_file dict.bin dict.bin.sorted dict.bin.sorted.swapped
echo
