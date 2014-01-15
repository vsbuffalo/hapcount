CC=clang++
CFLAGS=-g -O0

BAMTOOLS_ROOT=../bamtools
INCLUDE = -I$(BAMTOOLS_ROOT)/src
LIBS=-L$(BAMTOOLS_ROOT)/lib -lz -lm

all: hapcount

# builds bamtools static lib, and copies into root
$(BAMTOOLS_ROOT)/lib/libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && $(MAKE)

VariantProcessor.o: VariantProcessor.cpp VariantProcessor.h
	$(CC) $(CFLAGS) $(INCLUDE) -c VariantProcessor.cpp

hapcount.o: hapcount.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c hapcount.cpp

Variant.o: Variant.cpp Variant.h
	$(CC) $(CFLAGS) $(INCLUDE) -c Variant.cpp

ReadHaplotype.o: ReadHaplotype.cpp ReadHaplotype.h
	$(CC) $(CFLAGS) $(INCLUDE) -c ReadHaplotype.cpp

hapcount: hapcount.o VariantProcessor.o ReadHaplotype.o Variant.o $(BAMTOOLS_ROOT)/lib/libbamtools.a
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

clean:
	rm -f hapcount hapcount.o VariantProcessor.o ReadHaplotype.o Variant.o
