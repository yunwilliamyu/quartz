IDIR=.
CC=g++
LIB=/home/ywy/bin/lib
#CFLAGS=-I${IDIR} -I/usr/local/includes -std=c++11 -L${LIB} -Wall -O3 -DNDEBUG
CFLAGS=-I${IDIR} -I/home/ywy/bin/include -std=c++11 -pipe -L${LIB} -Wall -O3 -DNDEBUG -fopenmp \
#CFLAGS=-I${IDIR} -I/home/ywy/bin/include -std=c++11 -pipe -march=native -L${LIB} -Wall -O3 -DNDEBUG -fopenmp \
# -pg \
	   #-ftree-vectorizer-verbose=2

OBJS = $(patsubst %.cpp,%.o,$(wildcard *.cpp))

all: ${OBJS} misra_gries_dict dict_txt2bin sort_dict_file quartz 
	echo "All made."

quartz: quartz.o library.o jumpgate.o jumpgate.h
	${CC} -D_GLIBCXX_PARALLEL -o $@ ${CFLAGS} $@.o library.o jumpgate.o

sort_dict_file: sort_dict_file.o library.o jumpgate.o
	${CC} -D_GLIBCXX_PARALLEL -o $@ ${CFLAGS} $@.o library.o jumpgate.o

misra_gries_dict: misra_gries_dict.o
	${CC} -o $@ ${CFLAGS} $@.o

dict_txt2bin: dict_txt2bin.o library.o
	${CC} -o $@ ${CFLAGS} $@.o library.o

%.o: %.cpp 
	${CC} ${CFLAGS} -c -o $@ $<

clean:
	rm -f quartz misra_gries_dict dict_txt2bin ${OBJS} sort_dict_file 
	@echo "All cleaned up!"
