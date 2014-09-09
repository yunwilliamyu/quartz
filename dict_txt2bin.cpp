/* Copyright (c) 2013 Y. William Yu. All rights reserved. */
/* dict_txt2bin.cpp
 * Converts text based dictionary to binary 64-bit numbers (no counts)
 */

#include <stdio.h>
#include <fstream>
#include <string>
#include <global.h>

int main (int argc, char *argv[])
{
	if (argc < 3) {
		std::cerr << "Usage: " << argv[0] << " text_dict.db binary_dict.db" << std::endl;
	}

	std::ifstream text_dict (argv[1]);
	std::ofstream binary_dict (argv[2], std::ios::binary);
	std::string line;
	std::vector<readseq> curr_32mer;
	char line_string[4096];
	char *char_rep;
	readseq reverse;
	readseq final;

	unsigned long count = 0;
	char_rep = (char*)&count;
	binary_dict.write(char_rep,8);
	while(std::getline(text_dict, line))
	{
		++count;
		strcpy(line_string,line.c_str());
		curr_32mer = encode_read_vector(line_string);
		reverse = rev_compl(curr_32mer[0]);
		final = curr_32mer[0] < reverse ? curr_32mer[0] : reverse;
		char_rep = (char*)&final;
		binary_dict.write(char_rep,8);
	}
	binary_dict.seekp(0);
	char_rep = (char*)&count;
	binary_dict.write(char_rep,8);

}
