CC=clang++
CFLAGS=-g -O0

BAMTOOLS_ROOT=../bamtools
INCLUDE = -I$(BAMTOOLS_ROOT)/src
LIBS=-L$(BAMTOOLS_ROOT)/lib -lz -lm

all: hapcount

# builds bamtools static lib, and copies into root
$(BAMTOOLS_ROOT)/lib/libbamtools.a:
	cd $(BAMTOOLS_ROOT) && mkdir -p build && cd build && cmake .. && $(MAKE)

AlignmentProcessor.o: AlignmentProcessor.cpp AlignmentProcessor.h
	$(CC) $(CFLAGS) $(INCLUDE) -c AlignmentProcessor.cpp

hapcount.o: hapcount.cpp 
	$(CC) $(CFLAGS) $(INCLUDE) -c hapcount.cpp

hapcount: hapcount.o AlignmentProcessor.o $(BAMTOOLS_ROOT)/lib/libbamtools.a
	$(CC) $(CFLAGS) $(INCLUDE) -o $@ $^ $(LIBS)

clean:
	rm -f hapcount.o AlignmentProcessor.o
