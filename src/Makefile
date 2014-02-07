CC=clang++
CFLAGS=-g -O0 --std=c++11

BAMTOOLS_ROOT=../bamtools
INCLUDE = -I$(BAMTOOLS_ROOT)/src
LIBS=-L$(BAMTOOLS_ROOT)/lib -lz -lm

OBJ=hapcount.o VariantProcessor.o ReadHaplotype.o Variant.o \
			bamalignment_utilities.o VariantAnnotatedAlignment.o \
			$(BAMTOOLS_ROOT)/lib/libbamtools.a
all: hapcount

# builds bamtools static lib, and copies into root
$(BAMTOOLS_ROOT)/lib/libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && $(MAKE)

bamalignment_utilities.o: bamalignment_utilities.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

VariantAnnotatedAlignment.o: VariantAnnotatedAlignment.cpp VariantAnnotatedAlignment.h Variant.h
	$(CC) $(CFLAGS) $(INCLUDE) -c $<

VariantProcessor.o: VariantProcessor.cpp VariantProcessor.h Variant.h
	$(CC) $(CFLAGS) $(INCLUDE) -c VariantProcessor.cpp

hapcount.o: hapcount.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c hapcount.cpp

Variant.o: Variant.cpp Variant.h
	$(CC) $(CFLAGS) $(INCLUDE) -c Variant.cpp

ReadHaplotype.o: ReadHaplotype.cpp ReadHaplotype.h
	$(CC) $(CFLAGS) $(INCLUDE) -c ReadHaplotype.cpp

hapcount: $(OBJ)
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

clean:
	rm -f hapcount hapcount.o VariantProcessor.o ReadHaplotype.o Variant.o
