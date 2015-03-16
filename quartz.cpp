/* Copyright (c) 2013-2014 Y. William Yu. Released under CC0 1.0 Universal. */
/* jumpgate_sparsify.cpp
 * Sparsifies the quality vector for a read based on k-mer distance from dictionary
*/

//#define VERBOSE
#define BUFFER_SIZE 65536

#include "global.h"

#include <fstream>

#include "jumpgate.h"
#include <omp.h>

// Populates alter_this with the correct places based on the read specified in string
void compute_alter_this(bool *alter_this, char *str, read_entry_database &red) {
	readseq a;
	int lasty;
	std::vector<readseq> mer_list;
	mer_list = encode_read_vector(str);
	lasty=-100;
	for (unsigned int y=0; y<mer_list.size(); ++y) {
		if (((alter_this[y])||(y-lasty<16))&&(!(y==(mer_list.size()-1)))) {
			// Silently do nothing
		} else {
			a = mer_list[y];
			lasty=y;
			std::vector<int> snp;
			readseq ar = rev_compl(a);
			ar = rev_compl(a);
			if (a < ar) {
				snp = red.check_hamming_neighbors(a);
			} else {
				snp = red.check_hamming_neighbors(ar);
				for (unsigned int i = 0; i < snp.size(); ++i) {
					if (snp[i]>=0)
						snp[i] = 31-snp[i];
				}
			}
			if (snp.size()==0) {
				// Silently do nothing
			} else {
				std::vector<int> locations;
				locations.resize(32);
				for (auto it = snp.begin(); it!= snp.end(); ++it) {
					if (*it >= 0)
						locations[*it]=1;
				}
				for (int x=0; x<32; ++x) {
					alter_this[x+y] = alter_this[x+y] || !(locations[x]);
				}
			}
		}
	}
}

// Modifies read for FASTQ files
void fastq_walker (char * read, read_entry_database &red_whole, char qual) {
	char *c = &read[0];
	bool pass_through_read = false;
	int j = 0;
	int col = 1; // col here refers actually to a line of the 4 line FASTQ read entries
	bool alter_this[4096] = { 0 };
	char *str = '\0';
	
	while ( col <= 4 )
	{
		if (pass_through_read == true) {
			break;
		} else if (col == 2) {
			switch (*c) {
				case 'A': case 'a':
				case 'C': case 'c':
				case 'G': case 'g':
				case 'T': case 't':
				//case 'N':
					if (str == '\0') str = c;
					break;
				case '\t': case '\0': case '\n':
					//str[j]='\0';
					compute_alter_this(alter_this,str,red_whole);
					j=0;
					++col;
					break;
				default:
					// If something other than A,T,C,G (or N) appears, just skip line.
					pass_through_read = true;
					break;
			}
		} else if (col == 4) {
			switch (*c) {
				case '\n':
					pass_through_read = true;
					break;
	//			case '#':
	//				break;
				default:
					//if (((alter_this[j++])||(*c>qual))&&(*c!='#')) {
					if (((alter_this[j++]))&&(*c!='#')) {
					//if ((alter_this[j++])||(*c>qual)) {
						//*c = 126;
						*c = qual;
					}
					break;
			}

		} else {
			if (*c == '\n')
				++col;
		}
		c++;
	} 
}

// Modifies line for SAM files
// If qual='~', then we'll do an averaging instead
void linewalker (char * line, read_entry_database &red_whole, char qual) {
	// If qual is less than '~', then we'll jump ahead to the second pass
	bool second_pass;
	if (qual=='~')
		second_pass = false;
	else
		second_pass = true;

	
	char *c = &line[0];
	bool pass_through_line = false;
	int j = 0;
	int col = 1;
	bool alter_this[4096] = { 0 };
	char *read = '\0';
	char *qual_str = '\0';
	int qual_sum = 0;
	int qual_count = 0;

	
	// So we let headers pass through unmodified
	if (*c == '@') {
		pass_through_line = true;
	}
	while ( *c != '\n' )
	{
		if (pass_through_line == true) {
			break;
		} else if (col == 10) {
			switch (*c) {
				case 'A': case 'a':
				case 'C': case 'c':
				case 'G': case 'g':
				case 'T': case 't':
				//case 'N':
					if (read == '\0') read = c;
					break;
				case '\t': case '\0':
					//read[j]='\0';
					compute_alter_this(alter_this,read,red_whole);
					j=0;
					++col;
					break;
				case '\n':
					// This should NOT happen
					std::cerr << "Malformed read. Not enough columns." << std::endl;
					// Deliberate fall-through here.
				default:
					// If something other than A,T,C,G (or N) appears, just skip line.
					pass_through_line = true;
					break;
			}
		} else if (col == 11) {
			switch (*c) {
				case '\t':
					if (second_pass) {
						pass_through_line = true;
					} else {
						second_pass = true;
						c = qual_str-1;
						qual = qual_sum / qual_count + 1;
					}
					break;
				//case '#':
				//	break;
				default:
					if (qual_str == '\0') qual_str = c;
					if (second_pass) {
						//if (((alter_this[j++])||(*c>qual))&&(*c!='#')) {
						if (((alter_this[j++]))&&(*c!='#')) {
							//*c = 126;
							*c = qual;
						}
					} else {
						qual_sum += *c;
						qual_count++;
					}
					break;
			}

		} else {
			if (*c == '\t')
				++col;
		}
		c++;
	} 
}

// Finds all newlines in a string
inline std::vector<int> loc_newlines (char *str, int size) {
	std::vector<int> ans;
	for (int i=0; i<size; ++i) {
		if (str[i]=='\n') {
			ans.push_back(i);
		}
	}
	return ans;
}

// Get every fourth newline for SAM files
inline std::vector<int> keep_every_fourth (std::vector<int> &old) {
	std::vector<int> ans;
	for (unsigned int i=0; i<old.size(); ++i) {
		if (i%4==3)
			ans.push_back(old[i]);
	}
	return ans;
}

// Detects if file is a SAM file. This function is *very* brittle, and simply counts columns (tabs) in the line
bool is_sam_file(char *line) {
	char *c = line;
	int col = 1;
	while ((*c != '\0')&&(*c != '\n')) {
		if (*c++ == '\t') ++col;
	}
	return (col > 10);
}

int main( int argc, char *argv[])
{
	int64checker();
	// Make sure static variables are all properly initialized first
	if (argc < 5) {
		std::cerr << "Discards non-SNP quality values for known reads." << std::endl;
		std::cerr << "Usage: " << argv[0] << " [dictionary_file] [qual] [num_threads] [mem_option] input_file(s)" << std::endl;
		std::cerr << "\tInput is assumed to be either a FASTQ file." << std::endl
				  << std::endl
				  << "\tOutput will be modified FASTQ files named [input_file].filtered_[qual]," <<std::endl
				  << "\tidentical to the original except with sparsified quality scores." << std::endl << std::endl
				  << "\tThe quality of nearly all bases corresponding to a 32-mer listed in the" << std::endl
				  << "\tbinary [dictionary_file] will be set to [qual]. Correspondance shall be" << std::endl
				  << "\tdefined as a Hamming distance of less than or equal to 1 in the 32-mer." << std::endl
				  << "\tNote that a 32-mer can correspond to multiple 32-mers in the dictionary." << std::endl
				  << "\tThe exceptions to the setting of quality values above will be all bases" << std::endl
				  << "\tthat are different to any of the corresponding 32-mers in the dictionary," << std::endl
				  << "\twhich will retain their original quality value data." << std::endl << std::endl
				  << "\tThis implements the method described in [Yu, Yorukoglu, Peng, Berger, 2014]" << std::endl
				  << std::endl
				  << "\t[num_threads] specifies the number of threads that should be used. Be aware" << std::endl
				  << "\tthat disk I/O is usually the bottleneck for high thread number. In our" << std::endl
				  << "\ttests we found [num_thread]=8 to fully max out disk I/O. YMMV" <<std::endl
				  << std::endl
				  << "\t[mem_option] specifies the memory mode that should be used. If set to" << std::endl
				  << "\t'0', will use low-memory mode (<64GB); any other option will use high-" << std::endl
				  << "\tmemory-mode, which requires ~70GB of RAM" <<std::endl
				  << std::endl
				  << "This package is made available solely for academic use." << std::endl
				  << "Copyright (c) 2015 Yun William Yu. All rights reserved" << std::endl
				  << "Quartz version 0.2.1" << std::endl;
		exit(-1);
	}
	rev_compl(0);
	subst_find(0,0);

    bool lowmem = false;
    if (atoi(argv[4])==0)
        lowmem = true;

	std::cerr << time_now() << std::endl;
	std::string dict_filename = argv[1];
	read_entry_database red_whole(dict_filename, lowmem);

	
	std::FILE * inFile;
	std::FILE * outFile;

	char qual = argv[2][0];

	int num_threads = atoi(argv[3]);

	std::vector<readseq> mer_list;
	for (int argwalker = 5; argwalker < argc; ++argwalker) {
		std::string infile_name = argv[argwalker];

		if (infile_name == "-") {
			std::cerr << "Processing: Standard Input";
			inFile = stdin;
			infile_name = "stdin-";
			infile_name.append(time_now("%F_%H-%M-%S"));
		} else {
			std::cerr << "Processing: " << infile_name;
			inFile=fopen(infile_name.c_str(),"r");
			if (inFile==NULL) {
				perror ("Error opening input file.");
				exit(-3);
			}
		}
		std::string outfile_name = infile_name;
		outfile_name.append(".filtered_");
		outfile_name.append(argv[2]);
		// File reading buffers
		char * buffer;
		char * tmp_buffer;
		buffer = (char*) malloc(BUFFER_SIZE*2);
		tmp_buffer = (char*) malloc(BUFFER_SIZE*2);
		int num_read_bytes = 0;
		int buffer_filled = 0; // Number of active filled characters in buffer
		// Do a first read in order to determine if SAM file or FASTQ file
		/*
		num_read_bytes = fread(buffer + buffer_filled, 1, BUFFER_SIZE, inFile);
		buffer_filled = buffer_filled + num_read_bytes;
		std::vector<int> read_breaks = loc_newlines(buffer,buffer_filled);
		bool is_sam = is_sam_file(buffer+read_breaks[read_breaks.size()-2]); */
		// Autodetection isn't working right now. Assume SAM
		int thread_num;
		bool is_sam = false; // Or assume FASTQ
		//bool is_sam = true; // Or assume SAM
		if (is_sam) {
			thread_num = num_threads;
			std::cerr << " (SAM, using "<< thread_num <<" threads) ";
		} else {
			thread_num = num_threads;
			std::cerr << " (FASTQ, using "<< thread_num <<" threads) ";
		}
		omp_set_num_threads(thread_num);

		std::cerr << time_now() << " --> " ;
		setvbuf (inFile, NULL, _IOFBF, BUFSIZ);
		outFile=fopen(outfile_name.c_str(),"w");
		if (outFile==NULL) {
			perror ("Error opening output file.");
			exit(-4);
		}
		setvbuf (outFile, NULL, _IOFBF, BUFSIZ );

		// Progress counter
		long unsigned int counter = 0;
		unsigned char counter_mod = 0; // Because the loop isn't big enough
		printf("%016lu",counter);
		while (!feof(inFile)) { // This is deliberate, we only stop *after* we've tried to read past the end
			num_read_bytes = fread(buffer + buffer_filled, 1, BUFFER_SIZE, inFile);
			buffer_filled = buffer_filled + num_read_bytes;
			std::vector<int> read_breaks = loc_newlines(buffer,buffer_filled);
			if (!is_sam)
				read_breaks = keep_every_fourth(read_breaks);
#pragma omp parallel for
			for (unsigned int i=0; i<read_breaks.size(); ++i) {
				int start;
				if (i==read_breaks.size()-1)
					start = 0;
				else
					start = read_breaks[i]+1;
				if (!is_sam) {
					fastq_walker (buffer+start, red_whole, qual);
				} else {
					linewalker (buffer+start, red_whole, qual);
				}
			}
			int processed_bytes = read_breaks.back()+1;
			fwrite(buffer, 1, processed_bytes, outFile);
			memcpy(tmp_buffer, buffer + processed_bytes, buffer_filled - processed_bytes);
			memcpy(buffer,tmp_buffer, buffer_filled - processed_bytes);
			buffer_filled = buffer_filled - processed_bytes;
		//	if(counter>5000000) { std::cerr << " " << time_now() << std::endl; exit(0); }
			counter += read_breaks.size();
			if (!counter_mod++) {
				printf("\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b%016lu",counter);
				fflush(stdout);
			}
		}
		// Timer
		std::cerr << " " << time_now() << std::endl;

		fclose(inFile);
		fclose(outFile);
	}
	return 0;
};

