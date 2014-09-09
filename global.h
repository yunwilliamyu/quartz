/* Copyright (c) 2013 Y. William Yu. All rights reserved. */
#ifndef _GLOBAL_H_
#define _GLOBAL_H_

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <limits.h>
//#include <boost/dynamic_bitset.hpp>
//#include <bitset>
#include <algorithm>
#include <unordered_set>
#include <set>
#include <stdint.h>
#include <deque>
#include <ctime>

#include <fstream>

typedef uint64_t readseq;

/*
 * A = 00
 * C = 01
 * G = 10
 * T = 11
*/ 
std::deque<readseq> encode_read(const char *bases);
std::vector<readseq> encode_read_vector(const char *bases);
std::string decode_read(const readseq e);
readseq rev_compl(const readseq orig);
int subst_find(const readseq a, const readseq b);
//int subst_count(const std::bitset<MAXBITS> &a, const std::bitset<MAXBITS> &b);

void int64checker();

std::string time_now( const char* format = "%F %X"  );


#endif // #idndef _GLOBAL_H_
