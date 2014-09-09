/* Copyright (c) 2013 Y. William Yu. All rights reserved. */
/* library.cpp */

#include "global.h"

std::string time_now( const char* format )
{
	    std::time_t t = std::time(0) ;
	    char cstr[128] ;
	    std::strftime( cstr, sizeof(cstr), format, std::localtime(&t) ) ;
	    return cstr ;
}

void int64checker()
{
	if (CHAR_BIT!=8) {
		std::cerr << "Error: CHAR_BIT !=8." << std::endl;
		std::cerr << "This implementation assumes 8 bits per byte." << std::endl;
		std::cerr << "Your system is not compatibile. Please run on a standard x86-64 machine." << std::endl;		
		exit(-8);
	}
	if ((sizeof(unsigned long) * CHAR_BIT )!=64) {
		std::cerr << "Error: sizeof(unsigned long) != 64." << std::endl;
		std::cerr << "This implementation assumes 64 bit machines." << std::endl;
		std::cerr << "Your system is not compatibile. Please run on a standard x86-64 machine." << std::endl;		
		exit(-64);
	}
}

std::vector<readseq> encode_read_vector(const char *bases)
{
	std::vector<readseq> anslist;
	anslist.reserve(100);
	readseq ans=0;
	int i=0;
	while ((*bases))
	{
		switch (*bases++) {
			case 'A': case 'a': // store as 00
				ans >>= 2;
				break;
			case 'C': case 'c': // store as 01
				ans = (ans>>2)+(1ul<<62);
				break;
			case 'G': case 'g': // store as 10
				ans = (ans>>2)+(2ul<<62);
				break;
			case 'T': case 't': // store as 11
				ans = (ans>>2)+(3ul<<62);
				break;
			case '\n': case '\t': case 'N':
				return anslist;
			default: // if improper input, return empty
				std::vector<readseq> temp;
				return temp;
		}
		if (++i>=32)
			anslist.push_back(ans);
	}
	return anslist;
}

std::deque<readseq> encode_read(const char *bases)
{
	std::deque<readseq> anslist;
	readseq ans=0;
	int i=0;
	while ((*bases))
	{
		switch (*bases++) {
			case 'A': case 'a': // store as 00
				ans >>= 2;
				break;
			case 'C': case 'c': // store as 01
				ans = (ans>>2)+(1ul<<62);
				break;
			case 'G': case 'g': // store as 10
				ans = (ans>>2)+(2ul<<62);
				break;
			case 'T': case 't': // store as 11
				ans = (ans>>2)+(3ul<<62);
				break;
			case '\n': case '\t':
				return anslist;
			default: // if improper input, return empty
				std::deque<readseq> temp;
				return temp;
		}
		if (++i>=32)
			anslist.push_back(ans);
	}
	return anslist;
}

std::string decode_read(readseq e)
{
	char a[33];
	a[32]=0;
	for (int i=31; i>=0; --i)
	{
		switch ((e & (3ul<<62))>>62) {
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
		e <<= 2;
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

// Return negative number of differences if more than one diff
// Return position of difference if exactly one difference (0-indexed)
// Return -1 if two sequences are identical
int subst_find(const readseq a, const readseq b)
{
	static unsigned char table[65536] = {};
	static bool initialised = false;
	if (!initialised) {
		initialised = true;
		table[0]=0;
		for (unsigned int i = 0; i<65536; i++) {
			table[i]= (i & 1) + table[i / 2];
		}
	}

	// Don't need all this logic if the two are identical
	if (a == b)
		return -1;

	//printf("%lx\n%lx\n", a, b);
	readseq c;
	readseq d;
	c = a ^ b;	// XOR to find the differences
	d = c >> 1;		// SHIFT a copy by 1
	//printf("c=a^b:\t%lx\n d=c>>1: %lx\n", c, d);
	c |= d;		// OR the two together
	//printf("c |= d: %lx\n", c);
	//printf("mask: %lx\n", 0x5555555555555555);
	c &= 0x5555555555555555; // Mask out every other bit
	//printf("masked c: %lx\n", c); fflush(stdout);
	int count;
	uint16_t *c_16 = (uint16_t*)(&c);
	count = table[c_16[0]] + table[c_16[1]] + table[c_16[2]] + table[c_16[3]];
	//printf("count: %i\n", count);
	if (count>1)
		return (-count);
	else if (count<1)
		return -1;
	else {
		// Find location of set bit if exactly 1 difference
		for (unsigned int i = 0; i<32; ++i) {
			if (c&1)
				return (i);
			else 
				c >>= 2;
		}
		// Control should NEVER get here.
		std::cerr << "Internal sanity check failed in function subst_find" << std::endl;
		exit(-10);
	}

}

