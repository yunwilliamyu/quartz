/* Copyright (c) 2013 Y. William Yu. All rights reserved. */
/* sort_dict_file.cpp
 * Sorts a (binary) dictionary file
*/

#include "global.h"
#include <string>
#include <sstream>
#include <parallel/algorithm>

#include "jumpgate.h"

int main(int argc, char *argv[])
{
	if (argc < 4) {
		std::cerr << "Usage: " << argv[0] << " input_dict.bin output_dict.bin output_dict_swapped.bin" << std::endl;
		std::cerr << "\tSorts a dictionary file, outputs it, swaps high and low bits, sorts it, and outputs it" << std::endl;
		exit(-1);
	}
	std::string f_dict_name = argv[1];
	std::string f_out_name = argv[2];
	std::string f_out_swap_name = argv[3];


	std::cerr << "Loading dictionary" << std::endl;
	std::string fn = f_dict_name;
	// Read number of entries in dictionary
	std::ifstream f_dict_initial_read (fn, std::ios::in | std::ios::binary);
	unsigned long num_dict_entries = 0;
	f_dict_initial_read.read( (char*) &num_dict_entries, sizeof(num_dict_entries));

	// Read dictionary itself (note, dictionary is 0-indexed here)
	std::vector<uint64_t> dict;
	dict.resize(num_dict_entries);
	f_dict_initial_read.read( (char*) (&dict[0]), num_dict_entries*sizeof(uint64_t));

	std::cerr << "Sorting dictionary" << std::endl;
	__gnu_parallel::sort(dict.begin(),dict.end());
	std::ofstream f_out (f_out_name, std::ios::out | std::ios::binary | std::ios::trunc);
	f_out.write( (char*) &num_dict_entries, sizeof(num_dict_entries));
	std::cerr << "Writing dictionary to disk" << std::endl;
	f_out.write( (char*) (&dict[0]), num_dict_entries*sizeof(readseq));
	f_out.close();

	std::cerr << "Creating swapped dictionary ... ";
#pragma omp parallel for
	for (unsigned int i=0; i<num_dict_entries; ++i) {
		dict[i] = swap_halves(dict[i]);
	}
	std::cerr << " sorting ... " << std::endl;
	__gnu_parallel::sort(dict.begin(),dict.end());
	std::ofstream f_out_swap (f_out_swap_name, std::ios::out | std::ios::binary | std::ios::trunc);
	f_out_swap.write( (char*) &num_dict_entries, sizeof(num_dict_entries));
	std::cerr << "Writing dictionary to disk" << std::endl;
	f_out_swap.write( (char*) (&dict[0]), num_dict_entries*sizeof(readseq));
	f_out_swap.close();
	/*
	std::cerr << "Counting low bit buckets" << std::endl;
	std::vector<uint32_t> low_bit_counts;
	low_bit_counts.resize((1UL<<32)-1);
	for (uint64_t i = 0; i < dict.size(); ++i) {
		low_bit_counts[low_bits(dict[i])]++;
	}
	std::vector<uint32_t> jumpgate = low_bit_counts;
	for (uint64_t i = 1; i < jumpgate.size(); ++i) {
		jumpgate[i]+=jumpgate[i-1];
	}
	std::vector<uint32_t> dict_high;
	dict_high.resize(num_dict_entries);
	for (uint64_t i = 1; i < 
*/
	

//	std::vector<uint64_t> dict_swapped;


	
	
}


