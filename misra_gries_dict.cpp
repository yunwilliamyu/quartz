/* Copyright (c) 2013 Y. William Yu. All rights reserved. */
/* misra_gries_dict.cpp
   Does a misra-gries sketch with 2^32 counters to count the most frequent 32-mers appear.
   Outputs all 32-mers appearing with high frequency
 */
// For basically all masks except for the first two COUNTER_QUOTIENT = 64 is sufficient.
// Change COUNTER_QUOTIENT to 32 if you end up getting a high decrement rate
#define COUNTER_QUOTIENT 16
//#define COUNTER_QUOTIENT 64

#include <stdio.h>
#include <fstream>
#include <string>
#include <list>
#include <cstdio>
#include <utility>
#include <boost/functional/hash.hpp>
#include <unordered_set>
#include <unordered_map>
#include <map>
#include <set>
#include <vector>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef uint64_t readseq;
// Masks out to include only things whose middle four bases (8 bits) match mid. Set mid=-1 to avoid.
// Note, we use the old big endian versions of everything in this program, so it doesn't match with the newer little endian encodings
std::vector<readseq> encode_read_low(const char *bases, const int mid)
{
	std::vector<readseq> anslist;
	readseq ans=0;
	readseq rev=0;
	const readseq rev_mult =0x4000000000000000;
	int i=0;
	while ((*bases))
	{
		switch (*bases++) {
			case 'A': case 'a': // store as 00
				ans <<= 2;
				rev = (rev>>2)+3*rev_mult;
				break;
			case 'C': case 'c': // store as 01
				ans = (ans<<2)+1;
				rev = (rev>>2)+2*rev_mult;
				break;
			case 'G': case 'g': // store as 10
				ans = (ans<<2)+2;
				rev = (rev>>2)+1*rev_mult;
				break;
			case 'T': case 't': // store as 11
				ans = (ans<<2)+3;
				rev >>= 2;
				break;
			case '\n':
				return anslist;
			default: // if improper input, return empty
				std::vector<readseq> temp;
				return temp;
		}
		if (++i>=32) {
//			std::cerr << decode_read(ans) << std::endl;
//			std::cerr << decode_read(rev) << std::endl;
			//anslist.push_back(ans );
			readseq true_ans = ans < rev ? ans : rev;
			if ((mid < 0) || (mid > 255))
				anslist.push_back(true_ans);
			else if (((true_ans >> 28) & 0xFF) == (unsigned)mid)
				anslist.push_back(true_ans);
		}
	}
	return anslist;
}

// We use the old big endian version of everything in this program.
std::string decode_read(readseq e)
{
	char a[33];
	a[32]=0;
	for (int i=31; i>=0; --i)
	{
		switch (e & 3) {
			case 0:
				a[i]='A';
				break;
			case 1:
				a[i]='C';
				break;
			case 2:
				a[i]='G';
				break;
			case 3:
				a[i]='T';
				break;
		}
		e >>= 2;
	}
	std::string ans = a;
	return ans;
}

// This function *only* works on unsigned 64-bit integers!!!
readseq rev_compl(const readseq orig)
{
	static uint16_t table[65536] = {};
	static bool initialised = false;
	if (!initialised) {
		initialised = true;
		uint16_t x;
		uint16_t y;
		for (unsigned int i = 0; i<65536; i++) {
			x = 0;
			y = i;
			for (unsigned int j=0; j<8; j++) {
				x<<=2;
				switch ( y&3) {
					case 0:
						x+=3;
						break;
					case 1:
						x+=2;
						break;
					case 2:
						x+=1;
						break;
					case 3:
						break;
				}
				y>>=2;
			}
			table[i]=x;
			//printf (" 0x%X ,", x);
			//if (i%16 == 15) printf("\n");
		}
	}

	readseq ans;
	uint16_t *c_orig = (uint16_t*)(&orig);
	uint16_t *c_ans = (uint16_t*)(&ans);
	for (unsigned int i=0; i<4; ++i) {
		c_ans[i] = table[c_orig[3-i]];
	}
	return ans;
}

/*
//https://gist.github.com/badboy/6267743
uint32_t hash (uint64_t key)
{
  key = (~key) + (key << 18); // key = (key << 18) - key - 1;
  key = key ^ (key >> 31);
  key = key * 21; // key = (key + (key << 2)) + (key << 4);
  key = key ^ (key >> 11);
  key = key + (key << 6);
  key = key ^ (key >> 22);
  return (uint32_t) key;
}
*/

int main( int argc, char *argv[])
{
	//int64checker();
	if (argc < 4) {
		std::cerr<< "Usage: " << argv[0] << " MINCOUNT MASK output_file input_file(s)" << std::endl;
		std::cerr<< "\tApproximates the number of times frequent 32-mers appear, outputs two column" << std::endl;
		std::cerr<< "\tformat, with the first column specifying the 32-mer and the second column" << std::endl;
		std::cerr<< "\tthe number of times it appears in the corpus + decrement. MINCOUNT specifies" << std::endl;
		std::cerr<< "\tthe minimum number of appearances to be output in the database."  << std::endl << std::endl;
		std::cerr<< "\tMASK \\in [0,255] specifies the middle 4 base pairs. Set to -1 to get all."  << std::endl;
		exit(0);
       	}

	//omp_set_num_threads(16);
	unsigned int mincount = atoi(argv[1]);
	int mid = atoi(argv[2]);

	uint64_t max_el = 0x100000000; // total number of counters shouldn't exceed 2^32
	if ((mid >= -1) && (mid <= 255))
		max_el = max_el / COUNTER_QUOTIENT; // Decrease total number of counters if using mask

	uint32_t d = 0; // decrement counter
	std::vector<readseq> mer_list;
	/*
	std::map<uint64_t,uint16_t> db; // Keeps track of which B_list each k-mer is in
	std::set<uint64_t> B[65536]; */
	
	std::unordered_map<uint64_t,uint16_t, boost::hash<uint64_t> > db; // Keeps track of which B_list each k-mer is in
	db.rehash(max_el);
	db.reserve(max_el);
	std::unordered_set<uint64_t, boost::hash<uint64_t> > B[65536];
	
	/*
	google::sparse_hash_map<uint64_t,uint16_t, boost::hash<uint64_t> > db; // Keeps track of which B_list each k-mer is in
	db.set_deleted_key(0xFFFFFFFFFFFFFFFF);
//	db.resize(max_el);
	google::sparse_hash_set<uint64_t, boost::hash<uint64_t> > B[65536]; */
	/*for(int i=0;i<65536; i++)
		B[i].resize(65536);
	B[d+1].resize(max_el); */


	std::FILE * pFile;
	char line[4096];
	
	for (int argwalker = 4; argwalker < argc; ++argwalker) {
		if (argv[argwalker][0] != '-' ) {
			pFile=fopen(argv[argwalker],"r");
			setvbuf (pFile, NULL, _IOFBF, BUFSIZ);
			std::cout << "Processing: " << argv[argwalker] << std::endl;
		} else {
			std::cout << "Processing: standard input, mask " << mid << std::endl;
			pFile = stdin;
		}
		if (pFile==NULL) perror ("Error opening file");
		else
		{
			long int readcounter = 0;
			printf("Read number: %016d",0);
			while ( fgets (line, 4096, pFile) != NULL )
			{
				if (true) {
					if ((readcounter++%65536)==0) {
						printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%016ld", readcounter);
						fflush(stdout);
					}
					mer_list = encode_read_low(line, mid);
					//std::vector<readseq>::iterator it;
					unsigned int it;
					for (it = 0; it<mer_list.size(); ++it) {
						readseq i = mer_list[it];
						// If i already in hash, move the counter up one
						if (db.count(i)==1) {
							uint16_t x = db[i];
							if (x < 65535) {
								if (B[x].size()>1)
									B[x].erase(i);
								else
									B[x] = std::unordered_set<uint64_t, boost::hash<uint64_t> >();
								db[i] = x + 1;
								B[x+1].insert(i);
							}
						} else if (db.size()<max_el) { // If we have uninitialized counters
					//		#pragma omp critical
							{
							B[1].insert(i);
							db[i]=1;
							}
						} else {
							while (B[d].size()==0) { // If no more value "d" counters
								B[d] = std::unordered_set<uint64_t,  boost::hash<uint64_t> >();
								++d;
							}
							readseq t = *B[d].begin();
							db.erase(t);
							B[d].erase(t);
							B[d+1].insert(i);
							db[i]=d+1;
						}
					}
				}
			}
		}
		fclose (pFile);
		printf("\n");
	}

	std::ofstream outfile;
	outfile.open(argv[3]);
	outfile << "Dcremented:\t" << d << std::endl;
	//typedef google::sparse_hash_map<uint64_t,uint16_t, boost::hash<uint64_t> >::iterator it_type;
	//typedef std::map<uint64_t,uint16_t>::iterator it_type;
	typedef std::unordered_map<uint64_t,uint16_t, boost::hash<uint64_t> >::iterator it_type;
	for(it_type it=db.begin(); it != db.end(); ++it) {
		if (it->second - d >= mincount) {
			readseq tmp = it->first;
			outfile << decode_read(tmp) << '\t';
			outfile << (unsigned int) it->second << std::endl;
		}
	}
	outfile.close();

	return 0;
}

