/* Copyright (c) 2013-2014 Y. William Yu. Released under CC0 1.0 Universal. */
#ifndef _JUMPGATE_H_
#define _JUMPGATE_H_

#include "global.h"
#include <fstream>
#define BUFFER_SIZE 65536

// Creating jumpgate vector (for quick traversal of sorted dictionary)
std::vector<uint32_t> create_jumpgate(uint64_t *dict);

inline uint64_t swap_halves(uint64_t i) {
	return (i<<32)+(i>>32);
}
inline bool weird_lessthan ( uint64_t i, uint64_t j ) {
	return swap_halves(i) < swap_halves(j);
}
inline uint32_t low_bits(uint64_t i) {
	return (i&((1UL << 32)-1));
}
inline uint32_t high_bits(uint64_t i) {
	return (i>>32);
}

class read_entry_database {
	private:
		std::vector<uint32_t> dict_low; // Indexed by high bits
		std::vector<uint32_t> dict_high; // Indexed by low bits
		std::vector<uint32_t> jumpgate_low; // Indexed by high bits
		std::vector<uint32_t> jumpgate_high; // Indexed by low bits
		void dictionary_load(std::string dict_fn, std::vector<uint32_t> & dict_vec, std::vector<uint32_t> & jumpgate);
	public:
		// Note that dict_fn must contain a *sorted* dictionary of 64-bit integers, with first value being the number of remaining intgers in the file (i.e. size of dictionary) Also, dict_fn.swapped should be the same dictionary, but with low and higher bits swapped and appropriately resorted
		read_entry_database (std::string dict_fn, bool lm);
		int count_low (uint64_t x) const;
		int count_high (uint64_t x) const;
		std::vector<int> check_hamming_neighbors(const readseq work);
        // If lowmem = true, then count_high is just a call to count_low and saves memory
        bool lowmem;
};


// Returns a vector V with non-negative substitutions denoting positions of all
// new SNPs (i.e. some substitution at that point would make work match something
// in the dictionary). If there is a substitution at point i, then V[i]=i.
// Additionally, if there is an exact match, we use V[32]=-100). Non-specified points
// will be denoted with a -1. Thus, this will generally be a length 33 vector.
//
// Additionally, if there are no exact matches or Hamming neighbors in the dictionary,
// then we'll return an empty vector.
std::vector<int> check_red_vec(read_entry_database &red, const readseq work);

// Function which returns a set listing all new positions which are SNPs
// (i.e. the database has an entry exactly one substitution away from tested
// bitpattern shows currently set SNPs
//
// In the event that there are no hits in the database, return an empty vector
// Exact matches are denoted by the presence of a -1 in the resulting answer set.
std::set<int> check_red(read_entry_database &red, const readseq work);


#endif // #idndef _JUMPGATE_H_
