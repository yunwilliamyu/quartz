# FAQ
## How do I generate a high quality dictionary for compression?
In spirit, all that you need is a list of common 32-mers in a large corpus of reads drawn from the species of your choice. How you get that list doesn't matter so much really. A low quality dictionary you can get just by taking all the 32-mers in the reference genome for your species. Using a large corpus of reads and filtering for heavy hitters effectively gives you access to common variations too, which is useful for maintaining/improving accuracy.

The authors provide two tools for generating this list from a corpus of FASTQ files. generate_dict.cpp does this naively by counting all 32-mers and then filtering by some user-determined threshold for frequency. This works for small datasets, but will require extremely large amounts of RAM.

misra_gries_dict.cpp implements a version of the Misra-Gries approximate counting scheme and also allows counting of only a subset of 32-mers in each pass, which significantly decreases memory usage. The counting is no longer exact, but should still suffice for heavy hitters.

As a first approximation, given a corpus of reads from many individuals of the species with total depth-of-coverage D, choosing all 32-mers that appear around D/5 times seems to give relatively reasonable results.

dict_txt2bin.cpp converts a newline-separated list of 32-mers into the binary format Quartz uses.
