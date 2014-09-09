#!/usr/bin/env python

import sys
from string import split

f = open(sys.argv[1], 'r')
decrement_string = f.readline().split()
dec_count = int(decrement_string[1])
for line in f:
	x = line.split()
	print x[0] + "\t" + str(int(x[1])-dec_count)


