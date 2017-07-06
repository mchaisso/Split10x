HTSLIB=$(HOME)/software/htslib
ZLIB=$(HOME)/software/lib
all: split10x
#splitread

split10x: Split10xByPhase.cpp
	g++ -O3 -g -std=c++11 -D__STDC_LIMIT_MACROS -static  $^  -I $(HTSLIB)/htslib -L$(HTSLIB) -lhts  -lm -L$(ZLIB) -lz -lpthread  -o $@ 



