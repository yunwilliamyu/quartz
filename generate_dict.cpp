/* Copyright (c) 2013 Y. William Yu. All rights reserved. */
/* generate_dict.cpp
   Counts the number of times 32-mers appear.
   Outputs all 32-mers appearing with multiplicity greater than some threshold.
 */

#include <stdio.h>
#include <fstream>
#include <string>
#include <list>
#include <cstdio>
#include <unordered_map>
#include <utility>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <global.h>

int main( int argc, char *argv[])
{
	int64checker();
	if (argc < 4) {
		std::cerr<< "Usage: " << argv[0] << " MINCOUNT output_file input_file(s)" << std::endl;
		std::cerr<< "\tCounts the number of times 32-mers appear, and outputs it in a two column" << std::endl;
		std::cerr<< "\tformat, with the first column specifying the 32-mer and the second column" << std::endl;
		std::cerr<< "\tthe number of times it appears in the corpus. MINCOUNT specifies" << std::endl;
		std::cerr<< "\tthe minimum number of appearances to be output in the database."  << std::endl;
		exit(0);
       	}

	unsigned int *count;
	unsigned long countsize = 1024;
	count = (unsigned int*) malloc(countsize*sizeof(unsigned int));
	unsigned long charcount = 0;

	unsigned int mincount = atoi(argv[1]);

	readseq a;
	readseq b;
	std::deque<readseq> mer_list;
	std::unordered_map<readseq,unsigned long> db;

	std::FILE * pFile;
	int c;
	char str[1024];
	
	int i = 1;
	int j = 0;
	for (int argwalker = 3; argwalker < argc; ++argwalker) {
		pFile=fopen(argv[argwalker],"r");
		std::cout << "Processing: " << argv[argwalker] << std::endl;
		if (pFile==NULL) perror ("Error opening file");
		else
		{
			while ( (c=getc(pFile))!=EOF )
			{
				switch (c) {
					case '\t':
						if (++i == 11) {
							str[j]='\0';
							mer_list = encode_read(str);
							while (!mer_list.empty()) {
								a = mer_list.back();
								mer_list.pop_back();
								if (db.count(a)) {
									count[db.at(a)]+=1;
								}
								else {
									b = rev_compl( a );
									if (db.count(b)) {
										count[db.at(b)]+=1;
									}
									else {
										count[charcount]=1;
										db.emplace (a,charcount++);
										if (charcount == countsize) {
											count = (unsigned int*) realloc(count, countsize*2*sizeof(unsigned int));
											if (count!=NULL) {
												countsize*=2;
												std::cerr << "Allocating " << countsize << " bytes." <<  std::endl;
											}
											else {
												std::cerr << "Not enough memory. Could not allocate ";
												std::cerr << countsize*2 << " unsigned ints." << std::endl;
												exit(-2);
											}
										}
									}
								}
							}
						}
						break;
					case '\n':
						i=1;
						j=0;
						break;
					default:
						if (i==10)
						{
							switch (c) {
								case 'A':
								case 'C':
								case 'G':
								case 'T':
									str[j++]=(char)c;
									break;
								default:
									i=12;
									break;
							}
						}
						break;
				}
			}
		}
		fclose (pFile);
	}

	std::ofstream outfile;
	outfile.open(argv[2]);
	typedef std::unordered_map<readseq,unsigned long>::iterator it_type;
	for(it_type it=db.begin(); it != db.end(); ++it) {
		if (count[it->second] >= mincount) {
			readseq tmp = it->first;
			outfile << decode_read(tmp) << '\t';
			outfile << (unsigned int) count[it->second] << std::endl;
		}
	}
	outfile.close();

	return 0;
}

