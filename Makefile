################################################################################
# Makefile for dagSW software
# Gabor T. Marth 2015
# University of Utah
################################################################################

################################################################################
# To run (examples)
################################################################################

# optimized
# make CFLAGS="-O9 -static"

# optimized + static
# make CFLAGS="-O9 -static"

# optimized + static + 32bit
# make CFLAGS="-O9 -static -m32"

# to profile
# make CFLAGS="-O9 -static -pg"

################################################################################
# Variables
################################################################################

# Compiler
CC=g++

# Compiler flags
CFLAGS=

################################################################################
# Program builds
################################################################################
all: sw gsw

sw: sw.o Class-Alignment.o Class-Node.o
	$(CC) $(CFLAGS) sw.o Class-Alignment.o -lboost_regex -lpthread -o sw -g

gsw: gsw.o Class-GraphAlignment.o Class-Node.o
	$(CC) $(CFLAGS) gsw.o Class-Node.o Class-GraphAlignment.o -lboost_regex -lpthread -o gsw -g

################################################################################
# Compilation
################################################################################

# Main program compilations
sw.o: sw.cpp
	$(CC) $(CFLAGS) -c sw.cpp -o sw.o

gsw.o: gsw.cpp
	$(CC) $(CFLAGS) -c gsw.cpp -o gsw.o -g

Class-Alignment.o: Class-Alignment.cpp
	$(CC) $(CFLAGS) -c Class-Alignment.cpp -o Class-Alignment.o

Class-GraphAlignment.o: Class-GraphAlignment.cpp
	$(CC) $(CFLAGS) -c Class-GraphAlignment.cpp -o Class-GraphAlignment.o

Class-Node.o: Class-Node.cpp
	$(CC) $(CFLAGS) -c Class-Node.cpp -o Class-Node.o

################################################################################
# Archival
################################################################################
archive:	sw.cpp Class-Alignment.cpp Class-Alignment.h \
	gsw.cpp Class-GraphAlignment.cpp Class-GraphAlignment.h \
	 Class-Node.cpp  Class-Node.h Makefile
	tar cvf GSW.tar sw.cpp Class-Alignment.cpp Class-Alignment.h \
	gsw.cpp Class-GraphAlignment.cpp Class-GraphAlignment.h \
	Class-Node.cpp  Class-Node.h Makefile; \
	gzip GSW.tar

################################################################################
# Installation
################################################################################

################################################################################
# House-keeping
################################################################################
clean:
	rm -rf *.o sw gsw

