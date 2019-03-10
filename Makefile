CXX = g++
# FLAGS  = -march=native -O3 -funroll-loops -w
FLAGS  = -march=native -O3 -funroll-loops -Wno-unused-function -w -std=c++11 -pthread# -fopenmp
FLAGS1  = -march=native -O3 -funroll-loops -Wno-unused-function -w -std=c++11
# FLAGS  = -g -march=native -O3 -funroll-loops -Wno-unused-function -w -std=c++11 -pthread
# FLAGS  = -march=native -O3 -Wno-unused-function -w -std=c++11 -pthread
HEADER= src/BloomFilter.hpp src/ntHashIteratorSimple.hpp src/ntHashIterator.hpp \
		src/nthash.hpp src/khash.h src/ksort.h src/kvec.h

all: bfmem formatconvert compareMEM

bfmem: bfmem.o StopWatch.o
	$(CXX) $(FLAGS) src/bfmem.o src/StopWatch.o -o bfmem

bfmem.o: src/bfmem.cpp $(HEADER)
	$(CXX) -c $(FLAGS) src/bfmem.cpp  -o src/bfmem.o

StopWatch.o: src/StopWatch.cpp src/StopWatch.h
	$(CXX) -c $(FLAGS) src/StopWatch.cpp  -o src/StopWatch.o

formatconvert: formatconvert.o
	$(CXX) $(FLAGS1) src/formatconvert.o -o formatconvert

formatconvert.o: 
	$(CXX) -c $(FLAGS1) src/formatconvert.cpp -o src/formatconvert.o

compareMEM: compareMEM.o
	$(CXX) $(FLAGS1) src/compareMEM.o -o compareMEM

compareMEM.o: 
	$(CXX) -c $(FLAGS1) src/compareMEM.cpp -o src/compareMEM.o

clean:
	rm -f src/*.o a.out bfmem

all: bfmem 
