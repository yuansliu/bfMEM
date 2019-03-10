# include <iostream>
# include <cstdlib>
# include <fstream>
# include <algorithm>
# include <numeric>
# include <string>
# include <tuple>
# include <vector>
# include <queue>
# include <map>
# include <cassert>
# include <ctime>
# include <thread>
# include <mutex>
# include <unistd.h>
# include "BloomFilter.hpp"
# include "ntHashIterator.hpp"
# include "ntHashIteratorSimple.hpp"
# include "bfmem.h"
# include "StopWatch.h"
using namespace std;

typedef std::pair<std::string, size_t> SequenceItem;
typedef std::vector<SequenceItem> SequenceVector;

typedef std::tuple<std::string, char*, size_t> SequenceItem2;
typedef std::vector<SequenceItem2> SequenceVector2;

typedef std::tuple<size_t, char*, SequenceVector, std::vector<size_t> > GenomeData; //<size of buffer, memory pointer, starting pointer, sequence list>
//GenomeData is a tuple<size of buffer, memory pointer, information about such as `>chr1'>; the last is pointer to file

enum revcomp { no, yes, both };

vector<pair<string, string> > resultfile; //first is the filename, second is 
vector<int> splitfilenum;
map<string, pair<size_t, size_t> > smallfile;
mutex resultfilemtx, smallfilemtx, mapmtx;

const size_t MATCHES_BUFFER_SIZE = 1ULL << 24;
const size_t BUFFER_SIZE = 1ULL << 18;

int L, kmer, n_threads;

const int HASH_BITS = 20;
const int HASH_SIZE = 1 << HASH_BITS;
const int mask = HASH_SIZE - 1;

std::string matchesFN, smallMF;

std::ofstream f_matches;
std::ofstream fsmo;
mutex fmathcesmtx;
std::string R_FN;
std::string Q_FN;
string folder;

revcomp isRC;
char complement[256];

static const char alphanum[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
inline string generateString(const string &chr, int length = 5) {
	string res = chr + "_";
	for (int i = 0; i < length; ++i) {
		res += alphanum[rand() % 62];
	}
	return res;
}

inline void initPars() {
	srand(time(0));
	L = 100;
	// kmer = 60; //kmer is 4*i = kmer
	isRC = no;
	n_threads = 10;
	folder = generateString("bfmem", 10); //creat a temp folder for current input

	char *cmd = new char[200];
	sprintf(cmd, "mkdir -p %s", folder.c_str());
	system(cmd);
	delete[] cmd;

	folder += "/";
	smallMF = folder + "smallMatches";

	memset(complement, 125, 256);
	complement['A'] = 'T';
	complement['a'] = 'T';
	complement['C'] = 'G';
	complement['c'] = 'G';
	complement['G'] = 'C';
	complement['g'] = 'C';
	complement['T'] = 'A';
	complement['t'] = 'A';
	complement['N'] = 'N';
	complement['n'] = 'N';
}

std::ostream *v1logger;

inline void displayHelp(const char* prog) {
	printf("bfMEM v1.1, by Yuansheng Liu, March 2019.\n");
	printf("Usage: %s -r <ref_genome> -q <query_genome> -o <output> [option parameters]\n", prog);
	printf("\t options:\n \t\t -l <length> -k <length of k-mer> -s <strands> -t <threads>\n\n");
	// printf("-----------\n");
	printf("\t\t -r is the reference genome, a multi-FASTA file\n");
	printf("\t\t -q is the query genome, a multi-FASTA file\n");
	printf("\t\t -o is the output file\n");
	printf("\t\t -l is the minimum length of matches; default value is 100\n");
	printf("\t\t -k is the length of k-mer\n");
	printf("\t\t -t is the number of threads\n");
	printf("\t\t -s is the strands; default is forward;\n");
	printf("\t\t\t 'r' is reverse-complement;\n");
	printf("\t\t\t 'b' is both forward and reverse-complement;\n");
	printf("\t\t -h print help message\n");

	printf("Example:\n\t\t");
	printf("./bfmem -r H.all.fa -q M.all.fa -o hm-100.txt\n\n");
}

inline void displayParams() {
	printf("Reference file is: %s\n", R_FN.c_str());
	printf("Query file is: %s\n", Q_FN.c_str());
	printf("Minimum length of matches L = %d\n", L);
	if (isRC == no) {
		printf("Strands is forward\n");
	} else 
	if (isRC == yes) {
		printf("Strands is reversecomp complement\n");
	} else 
	if (isRC == both) {
		printf("Strands is both forward and reversecomp complement\n");
	}
}

inline void getPars(int argc, char* argv[]) {
	v1logger = &std::cout;
	bool is1 = false, is2 = false, is3 = false, iskmer = false; //four
	int oc;
	while ((oc = getopt(argc, argv, "r:q:o:l:s:k:t:hf")) >= 0) {
		switch (oc) {
			case 'l':
				L = atoi(optarg);
				// is0 = true;
				break;
			case 'r':
				R_FN = optarg;
				is1 = true;
				break;
			case 'q':
				Q_FN = optarg;
				is2 = true;
				break;
			case 'o':
				matchesFN = optarg;
				is3 = true;
				break;
			case 'k':
				kmer = atoi(optarg);
				iskmer = true;
				break;
			case 's':
				if (strcmp(optarg, "r") == 0) {
					isRC = yes;
				} else 
				if (strcmp(optarg, "b") == 0) {
					isRC = both;
				} else {
					std::cerr << "Error parameters.\n Please run 'bfmem -h'\n";
					exit(1);
				}
				break;
			case 't':
				n_threads = atoi(optarg);
				break;
			case 'h':
				displayHelp(argv[0]);
				exit(0);
			case '?':
				std::cerr << "Error parameters.\n Please run 'bfmem -h'\n";
				exit(1);
				break;
		}
	}

	if (!is1 || !is2 || !is3) {
		fprintf(stderr, "Required parameters are not provided!!\n\n");
		exit(1);
	}
	
	std::ifstream f;

	f.open(R_FN);
	if (f.fail()) {
		fprintf(stderr, "Reference file '%s' does not exist.\n", R_FN.c_str());
		exit(1);
	}
	f.close();

	f.open(Q_FN);
	if (f.fail()) {
		fprintf(stderr, "Query file '%s' does not exist.\n", Q_FN.c_str());
		exit(1);
	}
	f.close();

	f_matches.open(matchesFN);
	if (f_matches.fail()) {
		fprintf(stderr, "The output file '%s' can not be created.\n", matchesFN.c_str());
		exit(1);
	}
	f_matches.close();

	if (!iskmer) { //k is not set by user; set default value;
		// kmer = 60; //kmer is 4*i = kmer
		if (L >= 300) {
			kmer = 200;
		} else 
		if (L >= 200) {
			kmer = 160;
		} else 
		if (L >= 150) {
			kmer = 100;
		} else 
		if (L >= 100) {
			kmer = 60;
		} else 
		if (L >= 80) {
			kmer = 52;
		} else 
		if (L >= 50) {
			kmer = 40;
		} else 
		if (L >= 40) {
			kmer = 32;
		} else 
		if (L >= 34) {
			kmer = 28;
		} else {
			printf("L is too small, k must be provided by user.\n");
			exit(1);
		}
	}
	assert(kmer <= L && kmer%4==0);
}

/* this function is from https://github.com/wbieniec/copmem/blob/master/CopMEM.src/CopMEM.cpp
 * authors: Szymon Grabowski and Wojciech Bieniecki
 * changes: no
 */
inline bool arrcmp(size_t* p1, size_t* p2) {
	if (p1[1] < p2[1]) return true; 
	if (p1[1] > p2[1]) return false;
	return p1[0] < p2[0];
}
 
inline void dumpMEMTight(std::string &matchesBuffer, size_t* match) {
	matchesBuffer.append(std::to_string(match[0]) + " " +
						 std::to_string(match[1]) + " " + 
						 std::to_string(match[2]) + "\n");
}

/* this function is from https://github.com/wbieniec/copmem/blob/master/CopMEM.src/CopMEM.cpp
 * authors: Szymon Grabowski and Wojciech Bieniecki
 * changes: combine four append() to one
 */
inline void dumpMEMTight(std::string &matchesBuffer, SequenceItem& item1, size_t* match, size_t counter) {
	size_t baseindex1 = match[0] - item1.second;
	size_t baseindex2 = match[1] - counter;
	matchesBuffer.append(" " + item1.first + "\t" + 
						 std::to_string(baseindex1) + "\t" +
						 std::to_string(baseindex2) + "\t" +
						 std::to_string(match[2]) + "\n");
}

inline void dumpMEM(std::string &matchesBuffer, SequenceItem& item1, size_t* match, size_t counter) {
	size_t baseindex1 = match[0] - item1.second;
	size_t baseindex2 = match[1] - counter;
	// matchesBuffer.append(to_string(baseindex2) + " " + 
	matchesBuffer.append(
				item1.first + "\t" + 
				to_string(baseindex1) + "\t" + 
				to_string(baseindex2) + "\t" + 
				to_string(match[2]) + "\n");
}

/* this function is from https://github.com/wbieniec/copmem/blob/master/CopMEM.src/CopMEM.cpp
 * authors: Szymon Grabowski and Wojciech Bieniecki
 * changes: no
 */
inline bool lowerBoundComp(const SequenceItem &lhs, const SequenceItem &rhs) {
	return lhs.second < rhs.second;
}

/* this function is from https://github.com/wbieniec/copmem/blob/master/CopMEM.src/CopMEM.cpp
 * authors: Szymon Grabowski and Wojciech Bieniecki
 * changes: no
 */
inline SequenceItem findSeqDesc(size_t index, SequenceVector& seq) {
	SequenceItem dummySequenceItem = { "", index };
	SequenceItem item = seq[0];
	auto lower = std::lower_bound(seq.begin(), seq.end(), dummySequenceItem, lowerBoundComp);
	size_t diff = lower - seq.begin();
	return seq[diff - 1];
}

/* this function is from https://github.com/wbieniec/copmem/blob/master/CopMEM.src/CopMEM.cpp
 * authors: Szymon Grabowski and Wojciech Bieniecki
 * changes: delete some codes
 */
inline void postProcessTight(std::string &matchesBuffer, std::vector<size_t*> &matches, GenomeData *t1, std::string seqName, const size_t &counter) {
	matchesBuffer.append("> ");
	matchesBuffer.append(seqName);
	matchesBuffer.append("\n");

	if (matches.size() == 0) {
		// displayMatchInfo(seqName, 0);
		return;
	}
	char* gen1 = std::get<1>(*t1);
	SequenceVector& seq1 = std::get<2>(*t1);
	if (matches.size() > 1)
		std::sort(matches.begin(), matches.end(), arrcmp);

	SequenceItem seq1item;

	size_t* prev_match = matches[0];
	seq1item = findSeqDesc(prev_match[0], seq1);

	// auto foundPos = std::find(gen1 + prev_match[0], gen1 + prev_match[0] + prev_match[2], 'N');
	// if (foundPos == gen1 + prev_match[0] + prev_match[2])
	dumpMEMTight(matchesBuffer, seq1item,  prev_match, counter);

	std::uint64_t count = 1ULL;
	for (auto match: matches) {
		if (memcmp(prev_match, match, 3 * sizeof(size_t)) != 0) {
			// auto foundPos = std::find(gen1 + match[0], gen1 + match[0] + match[2], 'N');
			// if (foundPos != gen1 + match[0] + match[2])
				// continue;
			if (prev_match[0] != match[0])
				seq1item = findSeqDesc(match[0], seq1);
			dumpMEMTight(matchesBuffer, seq1item,  match, counter);
			++count;
		}
		prev_match = match;
	}
	// displayMatchInfo(seqName, count);
	for (auto match: matches)
		delete[] match;
	matches.clear();
}

/* this function change from previous one inline void postProcessTight(std::string &matchesBuffer, ...) */
inline void postProcessTight(std::string &matchesBuffer, std::vector<size_t*> &matches, GenomeData *t1, const size_t &counter) {
	if (matches.size() == 0) {
		// displayMatchInfo(seqName, 0);
		return;
	}
	char* gen1 = std::get<1>(*t1);
	SequenceVector& seq1 = std::get<2>(*t1);
	if (matches.size() > 1)
		std::sort(matches.begin(), matches.end(), arrcmp);

	SequenceItem seq1item;

	size_t* prev_match = matches[0];
	seq1item = findSeqDesc(prev_match[0], seq1);

	// auto foundPos = std::find(gen1 + prev_match[0], gen1 + prev_match[0] + prev_match[2], 'N');
	// if (foundPos == gen1 + prev_match[0] + prev_match[2])
	dumpMEMTight(matchesBuffer, seq1item,  prev_match, counter);

	std::uint64_t count = 1ULL;
	for (auto match: matches) {
		if (memcmp(prev_match, match, 3 * sizeof(size_t)) != 0) {
			// auto foundPos = std::find(gen1 + match[0], gen1 + match[0] + match[2], 'N');
			// if (foundPos != gen1 + match[0] + match[2])
				// continue;
			if (prev_match[0] != match[0])
				seq1item = findSeqDesc(match[0], seq1);
			dumpMEMTight(matchesBuffer, seq1item,  match, counter);
			++count;
		}
		prev_match = match;
	}
	// displayMatchInfo(seqName, count);
	for (auto match: matches)
		delete[] match;
	matches.clear();
}

/* this function change from previous one inline void postProcessTight(std::string &matchesBuffer, ...) */
inline void postProcessTight(std::string &matchesBuffer, std::vector<size_t*> &matches) {
	if (matches.size() == 0) {
		// displayMatchInfo(seqName, 0);
		return;
	}
	if (matches.size() > 1)
		std::sort(matches.begin(), matches.end(), arrcmp);

	size_t* prev_match = matches[0];
	dumpMEMTight(matchesBuffer, prev_match);

	std::uint64_t count = 1ULL;
	for (auto match: matches) {
		if (memcmp(prev_match, match, 3 * sizeof(size_t)) != 0) {
			dumpMEMTight(matchesBuffer, match);
		}
		prev_match = match;
	}
	// displayMatchInfo(seqName, count);
	for (auto match: matches)
		delete[] match;
	matches.clear();
}

inline void sortMatches(std::vector<size_t*> &matches) {
	if (matches.size() <= 1) return;
	std::sort(matches.begin(), matches.end(), arrcmp);
	size_t pre_idx = 0, matchessize = matches.size();
	for (size_t id = 1; id < matchessize; ++id) {
		if (memcmp(matches[pre_idx], matches[id], 3 * sizeof(size_t)) != 0) {
			matches[++pre_idx] = matches[id];
		} else {
			delete[] matches[id];
		}
	}
	matches.resize(pre_idx + 1);
}

/* this function is from https://github.com/wbieniec/copmem/blob/master/CopMEM.src/CopMEM.cpp
 * authors: Szymon Grabowski and Wojciech Bieniecki
 * changes: no
 */
inline void reverseComplement(char* start, const std::size_t N) {
	char* left = start + 1; // sequence starts from paddingChar
	char* right = start + N - 1;
	while (right > left) {
		char tmp = complement[*left];
		*left = complement[*right];
		*right = tmp;
		++left;
		--right;
	}
	if (left == right)
		*left = complement[*left];
}

/* the function change from inline void reverseComplement(char* start, const std::size_t N) */
inline void reverseComplement(char* start, char* end) {
	char* left = start + 1; // sequence starts from paddingChar
	char* right = end - 1;
	while (right > left) {
		char tmp = complement[*left];
		*left = complement[*right];
		*right = tmp;
		++left;
		--right;
	}
	if (left == right)
		*left = complement[*left];
}

/* this function is from https://github.com/wbieniec/copmem/blob/master/CopMEM.src/CopMEM.cpp
 * authors: Szymon Grabowski and Wojciech Bieniecki
 * changes: no
 */
inline void replaceBadSymbol(char* gen, char* dst, char symbol, char paddingChar) {
	char* movingPtr = gen;
	while (1) {
		char* tempPtr = std::find(movingPtr, dst, symbol);
		while (*tempPtr == symbol) {
			*tempPtr = paddingChar;
			++tempPtr;
		}
		if (tempPtr == dst)
			break;
		movingPtr = tempPtr + 1;
	}
}

/* this function is from https://github.com/wbieniec/copmem/blob/master/CopMEM.src/CopMEM.cpp
 * authors: Szymon Grabowski and Wojciech Bieniecki
 * changes: add a file point (the last element of GenomeData) used to reads chr
 */
inline GenomeData readMultiFasta(std::string fn, const char paddingChar, const bool removeNs, const char* seqType) {
	CStopWatch stopWatch;
	stopWatch.start();
	const char beginChar = '>';
	const char terminatorChar = 0;
	const char eolChar1 = 10;
	const char eolChar2 = 13;
	// const char spaceChar = ' ';
	const int paddingSize = L + 1;
	std::vector<size_t> ptr;

	*v1logger << "Reading " << seqType << " genome ...\n";  // seqType is "Reference" or "Query"

	//create a buffer for the whole file + padding at left and right
	std::ifstream f(fn, std::ios::ate | std::ios::binary);
	if (f.fail()) {
		std::cerr << "\nFile '" << fn << "' does not exist. Quit.";
		exit(1);
	}
	size_t N = f.tellg();

	f.seekg(0, std::ios::beg);
	char* buf1 = new char[N + 2 * paddingSize];
	if (buf1 == nullptr) {
		std::cerr << "\nFile '" << fn << "' is too large. Quit.";
		exit(1);
	}
	memset(buf1, paddingChar, paddingSize);
	f.read(buf1 + paddingSize, N);
	// *v2logger << "\treadMultiFasta: Reading file from disk " << stopWatch.stop() << "\n";
	stopWatch.resume();
	
	buf1[paddingSize + N] = terminatorChar; // null-terminate the string
	memset(buf1 + paddingSize + N + 1, paddingChar, paddingSize - 1);
	memset(buf1 + paddingSize + N, terminatorChar, 10);  // >= sizeof(uint_64) is enough
	f.close();

	char* gen = buf1 + paddingSize;
	SequenceVector seq;

	char* dst = gen;
	char* src = gen;

	char *tempLine = new char[512];  // FASTA lines shouldn't be that long (even 128 should be ok)

	while (1) {
		if (*src == beginChar) {
			ptr.push_back(src - gen);
			size_t idx = 0;
			while (*src != eolChar1 && *src != eolChar2 && *src != ' ' && *src != terminatorChar) {
				tempLine[idx] = *src++;
				++idx;
			}
			tempLine[idx] = 0;
			seq.push_back({ tempLine + 1, (dst - gen) });  // + 1, as we omit the starting '>'
														   //search for EOL
			while (*src != eolChar1 && *src != eolChar2)
				src++;

			*dst++ = paddingChar;
		}
		else {
			while (*src == eolChar1 || *src == eolChar2) {
				++src;
			}
			if (*src == beginChar)
				continue;
			if (*src == terminatorChar)
				break;

			uint64_t temp2;
			memcpy(&temp2, src, 8);
			while (((~temp2) & 0x4040404040404040) == 0) {
				temp2 = temp2 & 0xDFDFDFDFDFDFDFDF;
				memcpy(dst, &temp2, 8);
				dst += 8;
				src += 8;
				memcpy(&temp2, src, 8);
			}
			while (((~(*src)) & (char)0x40) == 0) {
				*dst++ = (*src++) & (char)0xDF;
			}
		}
	}
	ptr.push_back(src - gen);

	memset(dst, paddingChar, N + paddingSize - (dst - gen));

	if (removeNs == true) {
		replaceBadSymbol(gen, dst, 'N', paddingChar);
		replaceBadSymbol(gen, dst, 'n', paddingChar);
	}
	// *v2logger << "\treadMultiFasta: Analysing file " << stopWatch.stop() << "\n";
	// *v2logger << "\treadMultiFasta: Genome data size " << (dst - gen) << std::endl;
	// *v2logger << "\treadMultiFasta: Sequences " << seq.size() << std::endl;

	delete[] tempLine;
	return { (dst - gen), gen, seq, ptr };
}

/* this function change from the function readBlock(...) of https://github.com/wbieniec/copmem/blob/master/CopMEM.src/CopMEM.cpp
 * authors of readBlock: Szymon Grabowski and Wojciech Bieniecki
 */
inline SequenceItem2 readChr(std::ifstream &f, const size_t &beg, const size_t &bytesNum, const char &paddingChar, bool removeNs) {
	const char terminatorChar = 0;
	const char eolChar1 = 10;
	const char eolChar2 = 13;
	const char beginChar = '>';
	const int paddingSize = L + 1;
	char *buf1 = new char[bytesNum + paddingSize];
	if (buf1 == NULL) {
		std::cerr << "ERROR: buf1 is NULL\n";
		exit(1);
	}
	f.seekg(beg, std::ios::beg);
	f.read(buf1, bytesNum);
	memset(buf1 + bytesNum, 0, paddingSize);
	// fprintf(stderr, "%s\n", buf1);
	char *tempLine = new char[1024];
	string seqName;
	char *gen = buf1;
	char *src = gen;
	char *dst = gen;
	char *lastChar = gen + bytesNum;
	if (*src != beginChar) {
		std::cerr << "Invalid fasta File\n";
		exit(1);
	}
	size_t seqSize;

	size_t idx = 0;
	while (*src != eolChar1 && *src != eolChar2 && *src != ' ' && src != lastChar) { //reading header line
		tempLine[idx] = *src++;
		++idx;
	}
	tempLine[idx] = 0;
	seqName = tempLine + 1;

	while (*src != eolChar1 && *src != eolChar2) //search for EOL after header
			src++;
	char *startPtr = dst;
	*dst++ = paddingChar;  //padding char at the begin to make indexing matches from 1

	while (1) { //scanning the sequence
		while (*src == eolChar1 || *src == eolChar2) {
			if (src == lastChar)
				break;
			++src;
		}
		if (src == lastChar)
			break;

		uint64_t temp2; //scan 8 bytes at once
		memcpy(&temp2, src, 8);
		while (((~temp2) & 0x4040404040404040) == 0) {
			temp2 = temp2 & 0xDFDFDFDFDFDFDFDF;
			memcpy(dst, &temp2, 8);
			dst += 8;
			src += 8;
			memcpy(&temp2, src, 8);
		}
		while (((~(*src)) & (char)0x40) == 0) {
			*dst++ = (*src++) & (char)0xDF;
		}
	}
	// *dst = paddingChar;
	seqSize = dst - startPtr;
	memset(dst, paddingChar, bytesNum + paddingSize - seqSize);

	if (removeNs) {
		replaceBadSymbol(startPtr, dst, 'N', paddingChar);
		replaceBadSymbol(startPtr, dst, 'n', paddingChar);
	}
	delete[] tempLine;
	return {seqName, startPtr, seqSize};
}

/* this function is from https://github.com/wbieniec/copmem/blob/master/CopMEM.src/CopMEM.cpp
 * authors: Szymon Grabowski and Wojciech Bieniecki
 * changes: no
 */
inline void deleteReading(GenomeData & r) {
	delete[] (std::get<1>(r) - (L + 1));
}

inline void deleteChr(SequenceItem2 & r) {
	delete[] (std::get<1>(r));
}

// #define DEBUG

#ifdef DEBUG
int main(int argc, char* argv[]) {
	char seq[] = "AAAAAAAAAAAAAAAACATTAGATGAAATGCCCTCCAGTGGGGAGAGTGGAATTGTAGAGTTCTCCTCCCAGCAATCCTCTGAGAGTGGAATTGTAGAGTTCTCCTCCCAGCAATCCTCCAAAAAAAAAAAAAAAACATTAGATGAAATGCCCTCCAGTGGGGAGAGTGGAATTGTAGAGTTCTCCTCCCAGCAATCCTCTGAGAGTGGAATTGTAGAGTTCTCCTCCCAGCAATCCTC";
	// char seq[] = "AAAAAAAAAAAAAAAACATTAGATGAAATGCCCTCCAGTGGGGAGAGTGGAATTGTAGAGTTCTCCTCCCAGCAATCCTCTAAAAAAAAAAAAAAAACATTAGATGAAATGCCCTCCAGTGGGGAGAGTGGAATTGTAGAGTTCTCCTCCCAGCAATCCTCT";
	GenomeData  rGenome, qGenome;
	fprintf(stderr, "DEBUG\n");
	// initPars();
	L = 200;
	v1logger = &std::cout;
	// v2logger = &null_stream;

	fprintf(stderr, "%s\n", argv[1]);
	R_FN = argv[1];

	rGenome = readMultiFasta(R_FN, 123, true, "Reference");
	
	fprintf(stderr, "length of rGenome: %u\n", std::get<0>(rGenome));

	SequenceVector rseq = std::get<2>(rGenome);
	vector<size_t> ptr = get<3>(rGenome);

	std::ifstream f(R_FN, std::ios::binary | std::ifstream::in);

	for (int i = 0; i < rseq.size(); ++i) {
		// if (rseq[i].first == "IWGSC_CSS_6BS_scaff_2983071") {
		cout << rseq[i].first << endl;
		if (false && rseq[i].first == "TGAC_WGS_durum_v1_contig_3935499") {
			fprintf(stderr, "%s, %lu\n", (rseq[i].first).c_str(), rseq[i].second);
			fprintf(stderr, "%s, %lu\n", (rseq[i+1].first).c_str(), rseq[i+1].second);
			fprintf(stderr, "%lu, %lu\n", ptr[i], ptr[i+1]);
			SequenceItem2 chr = readChr(f, ptr[i], ptr[i+1] - ptr[i], 125, true);
			char *qstart = get<1>(chr);
			size_t qseqlen = get<2>(chr);
			fprintf(stderr, "%s\n", qstart);
			fprintf(stderr, "%lu\n", qseqlen);
		}
	}
	return 0;

	/*if (argc < 3) {
		displayHelp(argv[0]);
		return 1;
	} 
	initPars();
	getPars(argc, argv);
	qGenome = readMultiFasta(Q_FN, 125, true, "Query");
	char *qstart = std::get<1>(qGenome);
	fprintf(stderr, "%s\n", qstart);
	std::vector<size_t> ptr = std::get<3>(qGenome);
	for (int i = 0; i < ptr.size(); ++i) {
		fprintf(stderr, "%lu\n", ptr[i]);
	}

	std::ifstream f(Q_FN, std::ios::binary);
	// f.seekg(ptr[1], std::ios::beg);
	// char *genTemp = new char[templen+L];
	// f.read(genTemp, templen);
	// genTemp[templen] = '\0';
	// fprintf(stderr, "%s\n", genTemp);

	for (int i = 0; i < ptr.size() - 1; ++i) {
		size_t templen = ptr[i+1] - ptr[i];
		fprintf(stderr, "templen: %lu\n", templen);
		SequenceItem2 chr = readChr(f, ptr[i], templen, 125, true);
		fprintf(stderr, "---\n%s\n%s\n%lu\n---\n", std::get<0>(chr).c_str(), get<1>(chr), get<2>(chr));
		deleteChr(chr);
	}
	deleteReading(qGenome);*/
	// SequenceVector rseq = std::get<2>(qGenome);
	
	// BloomFilter(size_t expectedElemNum, double fpr, unsigned hashNum, unsigned kmerSize) :
	/*size_t seq_N = std::get<0>(qGenome);
	char *qstart = std::get<1>(qGenome);
	fprintf(stderr, "%s\n", qstart);
	const unsigned k = 81;
	const unsigned numHashes = 5;
	ntHashIterator itr(seq, numHashes, k, strlen(seq));
	// ntHashIterator itr(qstart, numHashes, k, seq_N);
	while (itr != itr.end()) {
		fprintf(stderr, "%lu %llu\n", itr.pos(), itr.getHash());
		++itr;
	}*/
}
#else

mutex *mtx, *rmtx, file_mtx;

struct SEGMENT {
	size_t b, e;
	size_t ptrb, ptre;
	SEGMENT(): b(0), e(0), ptrb(0), ptre(0) {}
	SEGMENT(const size_t &_b, const size_t &_e, const size_t &_ptrb, const size_t &_ptre): b(_b), e(_e), ptrb(_ptrb), ptre(_ptre)  {}
};

bool cmp(const SEGMENT &a, const SEGMENT &b) {
	return a.e-a.b > b.e-b.b;
}

// mutex nthashcntmtx;
// size_t nthashcnt = 0;

void iterate0(char *qstart, const int &qnumHashes, const int &L, const int &kmer, BloomFilter *qbf, BloomFilter *rqbf, vector<SEGMENT> arr) {
	for (int i = 0; i < arr.size(); ++i) {
		if (isRC != yes) {
			ntHashIteratorSimple qitr(qstart, arr[i].b, qnumHashes, L, kmer, arr[i].e);
			while (qitr != qitr.end()) {
				qbf->insert(*qitr);
				++qitr;
				
				// for test
				// nthashcntmtx.lock();
				// ++nthashcnt;
				// nthashcntmtx.unlock();
			}
		}

		if (isRC != no) { //isRC is yes
			// reverseComplement(qstart + arr[i].b, arr[i].e - arr[i].b + 1);
			reverseComplement(qstart + arr[i].b, qstart + arr[i].e);
			ntHashIteratorSimple qitr0(qstart, arr[i].b, qnumHashes, L, kmer, arr[i].e);
			while (qitr0 != qitr0.end()) {
				rqbf->insert(*qitr0);
				++qitr0;
			}
			// reverseComplement(qstart + arr[i].b - 1, arr[i].e - arr[i].b + 1); //
		}
	}
}

int idx = 0;

// mutex indexedMutex;
// size_t indexedNum = 0;

void iterate1(const char *rstart, vector<size_t> *th_pos, const int &qnumHashes, const int &kmer, BloomFilter *qbf, BloomFilter *rqbf, BloomFilter *rbf, BloomFilter *rrbf, mm_idx_t *htidx, mm_idx_t *rhtidx) {
	bool con;
	while (1) {
		int tidx = __sync_fetch_and_add(&idx, 1);
		if (tidx >= (*th_pos).size() - 1) break;

		ntHashIterator ritr(rstart, (*th_pos)[tidx], qnumHashes, kmer, (*th_pos)[tidx + 1] + kmer - 1); 
		// here, must + kmer - 1; as the class ntHashIterator can not process
		// the tail. Fix it is better. //

		while (ritr != ritr.end()) {
			if (isRC != yes) {
				con = qbf->contains(*ritr);
				if (con) {
					rbf->insert(*ritr);

					mm128_t mm128val(ritr.getHash(), ritr.pos());
					int val = mm128val.x & mask;
					mm128_v *p = &htidx->B[val].a;
					mtx[val].lock();
					kv_push(mm128_t, *p, mm128val);
					mtx[val].unlock();
					
					// indexedMutex.lock();
					// ++ indexedNum;
					// indexedMutex.unlock();
				} 
			}

			if (isRC != no) {
				con = rqbf->contains(*ritr);
				if (con) {
					rrbf->insert(*ritr);

					mm128_t mm128val(ritr.getHash(), ritr.pos());
					int val = mm128val.x & mask;
					mm128_v *p = &rhtidx->B[val].a;
					rmtx[val].lock();
					kv_push(mm128_t, *p, mm128val);
					rmtx[val].unlock();
					
					// indexedMutex.lock();
					// ++ indexedNum;
					// indexedMutex.unlock();
				} 
			}
			++ritr;
		}
	}
}

// mutex memMtx;
// size_t currSize;
// size_t maxChrSizeByte;

void iterate2(GenomeData *rGenome, const int &rnumHashes, const int &L, const int &kmer, BloomFilter *rbf, BloomFilter *rrbf, vector<SEGMENT> arr, mm_idx_t *htidx, mm_idx_t *rhtidx) {
	std::ifstream f(Q_FN, std::ios::binary | std::ifstream::in);
	char* rstart = std::get<1>(*rGenome);
	bool con;
	size_t arrsize = arr.size();
	bool *flag = new bool[arrsize];
	for (int i = 0; i < arrsize; ++i) {
		flag[i] = true;
	}

	size_t bptr, eptr;

	bool allflag = true;
	int cnt = 0, n;
	while (allflag) {
		allflag = false;
		for (int i = 0; i < arrsize; ++i) {
			if (flag[i]) {
				allflag = true;
				// memMtx.lock();
				// if (currSize + arr[i].ptre - arr[i].ptrb <= maxChrSizeByte) {
				// 	currSize += arr[i].ptre - arr[i].ptrb;
				// 	// fprintf(stderr, "i: %d; currSize: %lu\n", i, currSize);
				// }
				// memMtx.unlock();
				flag[i] = false;

				if (!flag[i]) {
					std::vector<size_t*> matches;
					std::string matchesBuffer;
					matchesBuffer.reserve(BUFFER_SIZE);
	
					SequenceItem2 chr = readChr(f, arr[i].ptrb, arr[i].ptre - arr[i].ptrb, 125, true);
					char *qstart = get<1>(chr);
					size_t qseqlen = get<2>(chr);
						// if (std::get<0>(chr) == "chr10") {
						// fprintf(stderr, "qseqlen: %lu\n", qseqlen);
						// int totalSum = 0, simpleSum = 0, hashSum = 0, kmerSum = 0;
					if (isRC != yes) { //&& get<0>(chr) == "IWGSC_CSS_6BS_scaff_2983071") {
						//for debug
						// CStopWatch stopwatch;
						// stopwatch.start();
						//for small chr write to final file directly
						std::ofstream f_seq_matches;
						string seqname = std::get<0>(chr);
						// bool debug = false;
						// if (seqname == "TGAC_WGS_durum_v1_contig_2580515") {
						// if (seqname == "IWGSC_CSS_6BS_scaff_2983071") {
						// 	debug = true;
						// }
						string seqfilename = folder + generateString(seqname);
						bool isopen = false; //first not open the file

						ntHashIteratorSimple fqitr(qstart, 0, rnumHashes, L, kmer, qseqlen);

						while (fqitr != fqitr.end()) {
							con = rbf->contains(*fqitr);
							if (con) {
								const uint64_t *r;
								r = mm_idx_get(htidx, fqitr.getHash(), &n);
								for (int i = 0; i < n; ++i) {
									char* op1 = rstart + (size_t)r[i];
									char* op2 = qstart + fqitr.pos();
									if (memcmp(op1, op2, kmer) == 0){
											// ++kmerSum;
										char* p1 = op1 + kmer;
										char* p2 = op2 + kmer;
											// ++totalSum;
										while (memcmp(p1, p2, 4) == 0) {
											p1 += 4, p2 += 4;
										}

										while (*p1 == *p2) {
											++ p1, ++ p2;
										}
										char *right = p1;

										p1 = op1; p2 = op2;
										while (*p1 == *p2) {
											-- p1, -- p2;
										}

										if (right - p1 >= L + 1) {
											size_t* tempMatch = new size_t[3];
											tempMatch[0] = p1 + 1 - rstart;
											tempMatch[1] = p2 + 1 - qstart;
											tempMatch[2] = right - p1 - 1;

											//a simple remove duplicate, compare with previous one
											// matches.push_back(tempMatch);
											
											int matchessize = matches.size();
											// if (matchessize == 0 || (matchessize > 0 && memcmp(tempMatch, matches[matchessize - 1], 3 * sizeof(size_t)) != 0)) {
												matches.push_back(tempMatch);
											// } else {
											// 	delete[] tempMatch;
											// }
											
											// if (debug) fprintf(stderr, "%lu\n", matches.size());
											if (matches.size() > 500) {
												sortMatches(matches);
												if (matches.size() > 300) {
													postProcessTight(matchesBuffer, matches);

													if (!isopen) {
														resultfilemtx.lock();
														resultfile.push_back({seqname, seqfilename});
														resultfilemtx.unlock();

														f_seq_matches.open(seqfilename);
														isopen = true;
													}

													f_seq_matches << matchesBuffer;
													matchesBuffer.clear();
												}
											}
										}
									}
								}
							}
							++fqitr;
						}
						// if (debug) fprintf(stderr, "finale matches.size(): %lu\n", matches.size());

						if (isopen) {
							if (matches.size() > 0) {
								postProcessTight(matchesBuffer, matches);
								f_seq_matches << matchesBuffer;
								matchesBuffer.clear();
							}
							f_seq_matches.close();
						} else {
							if (matches.size() > 0) {
								postProcessTight(matchesBuffer, matches, rGenome, 0);
								
								smallfilemtx.lock();
								bptr = fsmo.tellp();
								fsmo << matchesBuffer;
								eptr = fsmo.tellp();
								smallfilemtx.unlock();

								matchesBuffer.clear();

								mapmtx.lock();
								smallfile[seqname] = {bptr, eptr - bptr};
								mapmtx.unlock();
							}
						}
						
						// if (isopen) {
							// if (debug) fprintf(stderr, "isopen\n");
						/*if (matches.size() > 0) {
							if (!isopen) {
								resultfilemtx.lock();
								resultfile.push_back({seqname, seqfilename});
								resultfilemtx.unlock();

								f_seq_matches.open(seqfilename);
								isopen = true;
							}

							postProcessTight(matchesBuffer, matches);
							f_seq_matches << matchesBuffer;
							matchesBuffer.clear();
						}
						if (isopen) {
							f_seq_matches.close();
						}*/
						// } else { //do not open a file
						// 	// if (debug) fprintf(stderr, "not open\n");
						// 	postProcessTight(matchesBuffer, matches, rGenome, seqname, 0);
						// 	fmathcesmtx.lock();
						// 	f_matches << matchesBuffer;
						// 	fmathcesmtx.unlock();
						// }

						// *v1logger << "Time of " << std::get<0>(chr).c_str() << " = " << stopwatch.stop() << std::endl;
					}
	
					if (isRC != no) {///
						reverseComplement(qstart, qseqlen);
						ntHashIteratorSimple fqitr0(qstart, 0, rnumHashes, L, kmer, qseqlen);

						std::ofstream f_seq_matches;
						string seqname = std::get<0>(chr) + " Reverse";
						string seqfilename = folder + generateString(std::get<0>(chr) + "_Reverse");

						bool isopen = false;

						while (fqitr0 != fqitr0.end()) {
							con = rrbf->contains(*fqitr0);
							if (con) {
								const uint64_t *r;
								r = mm_idx_get(rhtidx, fqitr0.getHash(), &n);
								for (int i = 0; i < n; ++i) {
									// if ?
									char* op1 = rstart + (size_t)r[i];
									char* op2 = qstart + fqitr0.pos();
									if (memcmp(op1, op2, kmer) == 0){
										char* p1 = op1 + kmer;
										char* p2 = op2 + kmer;

										while (memcmp(p1, p2, 4) == 0) {
											p1 += 4, p2 += 4;
										}

										while (*p1 == *p2) {
											++ p1, ++ p2;
										}
										char *right = p1;

										p1 = op1; p2 = op2;
										while (*p1 == *p2) {
											-- p1, -- p2;
										}
										if (right - p1 >= L + 1) {
											size_t* tempMatch = new size_t[3];
											tempMatch[0] = p1 + 1 - rstart;
											tempMatch[1] = p2 + 1 - qstart;
											tempMatch[2] = right - p1 - 1;
											
											// matches.push_back(tempMatch);

											int matchessize = matches.size();
											// if (matchessize == 0 || (matchessize > 0 && memcmp(tempMatch, matches[matchessize - 1], 3 * sizeof(size_t)) != 0)) {
											matches.push_back(tempMatch);
											// } else {
											// 	delete[] tempMatch;
											// }
											
											// if (debug) fprintf(stderr, "%lu\n", matches.size());
											if (matches.size() > 500) {
												sortMatches(matches);
												if (matches.size() > 300) {
													postProcessTight(matchesBuffer, matches);

													if (!isopen) {
														resultfilemtx.lock();
														resultfile.push_back({seqname, seqfilename});
														resultfilemtx.unlock();

														f_seq_matches.open(seqfilename);
														isopen = true;
													}

													f_seq_matches << matchesBuffer;
													matchesBuffer.clear();
												}
											}
										}
									}
								}
							}
							++fqitr0;
						}
						// if (isopen) {
						// 	if (matches.size() > 0) {
						// 		postProcessTight(matchesBuffer, matches);
						// 		f_seq_matches << matchesBuffer;
						// 	}
						// 	f_seq_matches.close();
						// } else {
						// 	postProcessTight(matchesBuffer, matches, rGenome, seqname, 0);
						// 	fmathcesmtx.lock();
						// 	f_matches << matchesBuffer;
						// 	fmathcesmtx.unlock();
						// }
						// matchesBuffer.clear();

						if (isopen) {
							if (matches.size() > 0) {
								postProcessTight(matchesBuffer, matches);
								f_seq_matches << matchesBuffer;
								matchesBuffer.clear();
							}
							f_seq_matches.close();
						} else {
							if (matches.size() > 0) {
								postProcessTight(matchesBuffer, matches, rGenome, 0);
								
								smallfilemtx.lock();
								bptr = fsmo.tellp();
								fsmo << matchesBuffer;
								eptr = fsmo.tellp();
								smallfilemtx.unlock();

								matchesBuffer.clear();

								mapmtx.lock();
								smallfile[seqname] = {bptr, eptr - bptr};
								mapmtx.unlock();
							}
						}
					}

					deleteChr(chr);
				}
			}
		}
	}
}

#define sort_key_128x(a) ((a).x)
KRADIX_SORT_INIT(128x, mm128_t, sort_key_128x, 8) 
KSORT_INIT_GENERIC(uint32_t)

mutex fidmtx;
int fid = 0;

inline void split(string filename, SequenceVector &seq1, int &filenum);

void splitTempFileThread(SequenceVector *seq1) {
	int tfid;
	while (true) {
		tfid = __sync_fetch_and_add(&fid, 1);
		if (tfid >= resultfile.size()) break;
		split(resultfile[tfid].second, *seq1, splitfilenum[tfid]);
	}
}

// int(*hashFunc[130])(const size_t &a, const size_t &b, const size_t &c);
int(*hashFunc[130])(const size_t &b);

inline void initHashFunc() {
	// hashFunc[1] = hash1;
	hashFunc[2] = hash2;
	hashFunc[4] = hash4;
	hashFunc[8] = hash8;
	hashFunc[16] = hash16;
	hashFunc[32] = hash32;
	hashFunc[64] = hash64;
	hashFunc[128] = hash128;
}

vector<string> filesvec;
mutex filesvecmtx;
size_t filesvecidx = 0;

inline void split(string filename, SequenceVector &seq1, int &filenum) {
	std::ifstream f(filename, std::ios::ate | std::ios::binary);
	size_t N = f.tellg(); //bytes
	f.close();
	filenum = 0; // no matches
	if (N == 0) return;

	N /= 1000000; //kb --> MB
	N /= 300;//300M
	N += 1;
	
	filenum = 1;
	// if (N >= 64) { // or setting 128
	// 	filenum = 64;
	// }
	if (N >= 128) { // if error for opening much files; change 128 to 64 or "ulimit -SHn 65535"
		filenum = 128;
	}
	 else {
		while (filenum < N) {
			filenum <<= 1;
		}
	}

	if (filenum == 1) {	// sort and findSeqDesc
		size_t a, b, c;
		std::vector<size_t*> matches;
		FILE *fp = fopen(filename.c_str(), "r");
		while (fscanf(fp, "%lu%lu%lu", &a, &b, &c) != EOF) {
			size_t *tempMatch = new size_t[3];
			tempMatch[0] = a, tempMatch[1] = b, tempMatch[2] = c;
			matches.push_back(tempMatch);
		}
		fclose(fp);

		if (matches.size() > 1) {
			std::sort(matches.begin(), matches.end(), arrcmp);
		}

		// findSeqDesc
		if (matches.size() > 0) {
			string matchesBuffer;
			matchesBuffer.reserve(MATCHES_BUFFER_SIZE);

			std::ofstream f_chr_matches;
			f_chr_matches.open(filename);

			size_t *prev_match = matches[0];

			SequenceItem seq1item = findSeqDesc(prev_match[0], seq1);
			dumpMEMTight(matchesBuffer, seq1item, prev_match, 0);
			// fprintf(splitfp, "%lu %lu %lu\n", prev_match[0], prev_match[1], prev_match[2]);
			for (auto match: matches) {
				if (memcmp(prev_match, match, 3 * sizeof(size_t)) != 0) {
					// fprintf(splitfp, "%lu %lu %lu\n", match[0], match[1], match[2]);
					seq1item = findSeqDesc(match[0], seq1);
					dumpMEMTight(matchesBuffer, seq1item, match, 0);
					if (matchesBuffer.size() > (size_t(MATCHES_BUFFER_SIZE*0.95))) {
						f_chr_matches << matchesBuffer; 
						matchesBuffer.clear();
					}
				}
				prev_match = match;
			}
			if (matchesBuffer.size() > 0) {
				f_chr_matches << matchesBuffer; 
				matchesBuffer.clear();
			}
			f_chr_matches.close();
		}
		for (auto match: matches) {
			delete[] match;
		}
		matches.clear();
	} else { //filenum > 1; //only split
		string *splitfilename = new string[filenum];
		FILE **splitfp = new FILE*[filenum];
		string subfolder = filename + "_split";

		char *cmd = new char[1024];
		sprintf(cmd, "mkdir -p %s", subfolder.c_str());
		system(cmd);
		delete[] cmd;

		for (int i = 0; i < filenum; ++i) {
			splitfilename[i] = subfolder + "/part" + to_string(i);
			splitfp[i] = fopen(splitfilename[i].c_str(), "w");
			if (splitfp[i] == NULL) {
				fprintf(stderr, "open file fail!\n");
				exit(1);
			}
		}

		filesvecmtx.lock();
		for (int i = 0; i < filenum; ++i) {
			filesvec.push_back(splitfilename[i]);
		}
		filesvecmtx.unlock();
		// fprintf(stderr, "open parts file success!\n");

		size_t a, b, c;
		FILE *fp = fopen(filename.c_str(), "r");
		// if (fp == NULL) {
		// 	fprintf(stderr, "fp is NULL\n");
		// 	exit(1);
		// }

		while (fscanf(fp, "%lu%lu%lu", &a, &b, &c) != EOF) {
			// int val = hashFunc[filenum](a, b, c);
			int val = hashFunc[filenum](b);
			fprintf(splitfp[val], "%lu %lu %lu\n", a, b, c);
		}
		for (int i = 0; i < filenum; ++i) {
			fclose(splitfp[i]);
		}
		fclose(fp);

		delete[] splitfilename;
		delete[] splitfp;
	}	

}

const size_t MATCHES_BUFFER_SIZE1 = 1 << 28;

inline void sortAndFindSeqInfo(string &filename, SequenceVector &seq1) {
	size_t a, b, c;
	std::vector<size_t*> matches;
	FILE *fp = fopen(filename.c_str(), "r");
	while (fscanf(fp, "%lu%lu%lu", &a, &b, &c) != EOF) {
		size_t *tempMatch = new size_t[3];
		tempMatch[0] = a, tempMatch[1] = b, tempMatch[2] = c;
		matches.push_back(tempMatch);
	}
	fclose(fp);

	if (matches.size() > 1) {
		std::sort(matches.begin(), matches.end(), arrcmp);
	}

	// findSeqDesc
	if (matches.size() > 0) {
		string matchesBuffer;
		matchesBuffer.reserve(MATCHES_BUFFER_SIZE);

		std::ofstream f_chr_matches;
		f_chr_matches.open(filename);

		size_t *prev_match = matches[0];

		SequenceItem seq1item = findSeqDesc(prev_match[0], seq1);
		dumpMEMTight(matchesBuffer, seq1item, prev_match, 0);
		// dumpMEM(matchesBuffer, seq1item, prev_match, 0);
		// fprintf(splitfp, "%lu %lu %lu\n", prev_match[0], prev_match[1], prev_match[2]);
		for (auto match: matches) {
			if (memcmp(prev_match, match, 3 * sizeof(size_t)) != 0) {
				// fprintf(splitfp, "%lu %lu %lu\n", match[0], match[1], match[2]);
				seq1item = findSeqDesc(match[0], seq1);
				dumpMEMTight(matchesBuffer, seq1item, match, 0);
				// dumpMEM(matchesBuffer, seq1item, match, 0);
				if (matchesBuffer.size() > (size_t(MATCHES_BUFFER_SIZE*0.95))) {
					f_chr_matches << matchesBuffer; 
					matchesBuffer.clear();
				}
			}
			prev_match = match;
		}
		if (matchesBuffer.size() > 0) {
			f_chr_matches << matchesBuffer; 
			matchesBuffer.clear();
		}
		f_chr_matches.close();
	}
	for (auto match: matches) {
		delete[] match;
	}
	matches.clear();
}

void sortAndFindSeqInfoThread(SequenceVector *seq1) {
	size_t tfilesvecidx;
	while (true) {
		tfilesvecidx = __sync_fetch_and_add(&filesvecidx, 1);
		if (tfilesvecidx >= filesvec.size()) break;
		sortAndFindSeqInfo(filesvec[tfilesvecidx], *seq1);
	}
}

void mergeSortedFiles(SequenceVector &qseq) {
	size_t bytesNum = 1UL<<29;
	char* buf1 = new char[bytesNum + 2];
	map<string, int> chrname2k;
	for (int k = 0; k < resultfile.size(); ++k) {
		chrname2k[resultfile[k].first] = k;
	}

	ifstream smfi;
	smfi.open(smallMF);
	pair<size_t, size_t> ptr;

	for (int ii = 0; ii < qseq.size(); ++ii) {
		if (isRC != yes) {

			string mbuf;
			mbuf.append("> ");
			mbuf.append(qseq[ii].first);
			mbuf.append("\n");
			f_matches << mbuf;
			mbuf.clear();

			if (chrname2k.find(qseq[ii].first) != chrname2k.end()) {

				int k = chrname2k[qseq[ii].first];

				if (splitfilenum[k] > 1) {
					// string folder = resultfile[k].second + "_split";
					// vector<string> files;
					// for (int i = 0; i < splitfilenum[k]; ++i) {
					// 	files.push_back(folder + "/part" + to_string(i));
					// }
					// mergeSortedFiles(files);

					// ----
					string folder = resultfile[k].second + "_split";
					for (int i = 0; i < splitfilenum[k]; ++i) {
						string splitfilename = folder + "/part" + to_string(i);
						std::ifstream f(splitfilename, std::ios::ate | std::ios::binary);
						if (f.fail()) {
							fprintf(stderr, "error\n");
							exit(0);
						}
						size_t N = f.tellg();
						// std::cerr << "N: " << N << "\n";
						size_t t = 0, tmp;
						f.seekg(0, std::ios::beg);	
						while (t < N) {
							tmp = N - t;
							if (tmp > bytesNum) {
								tmp = bytesNum;
							}
							t += tmp;
							f.read(buf1, tmp);
							buf1[tmp] = '\0';

							f_matches << buf1;
						}
						f.close();
					}
				} else 
				if (splitfilenum[k] == 1) {
					std::ifstream f(resultfile[k].second, std::ios::ate | std::ios::binary);
					if (f.fail()) {
						fprintf(stderr, "error\n");
						exit(0);
					}
					size_t N = f.tellg();
					// std::cerr << "N: " << N << "\n";
					size_t t = 0, tmp;
					f.seekg(0, std::ios::beg);	
					while (t < N) {
						tmp = N - t;
						if (tmp > bytesNum) {
							tmp = bytesNum;
						}
						t += tmp;
						f.read(buf1, tmp);
						buf1[tmp] = '\0';

						f_matches << buf1;
					}
					f.close();
				}
			} else 
			if (smallfile.find(qseq[ii].first) != smallfile.end()) {
				ptr = smallfile[qseq[ii].first];
				// fprintf(stderr, "%lu %lu\n", ptr.first, ptr.second);

				smfi.seekg(ptr.first, std::ios::beg);
				smfi.read(buf1, ptr.second);
				buf1[ptr.second] = '\0';

				f_matches << buf1;
			}
			
		}

		if (isRC != no) {

			string fn = qseq[ii].first + " Reverse";
			string mbuf;
			mbuf.append("> ");
			mbuf.append(fn);
			mbuf.append("\n");
			f_matches << mbuf;
			mbuf.clear();

			if (chrname2k.find(fn) != chrname2k.end()) {
				int k = chrname2k[fn];
				if (splitfilenum[k] > 1) {
					// string folder = resultfile[k].second + "_split";
					// vector<string> files;
					// for (int i = 0; i < splitfilenum[k]; ++i) {
					// 	files.push_back(folder + "/part" + to_string(i));
					// }
					// mergeSortedFiles(files);

					// 
					string folder = resultfile[k].second + "_split";
					for (int i = 0; i < splitfilenum[k]; ++i) {
						string splitfilename = folder + "/part" + to_string(i);
						std::ifstream f(splitfilename, std::ios::ate | std::ios::binary);
						if (f.fail()) {
							fprintf(stderr, "error\n");
							exit(0);
						}
						size_t N = f.tellg();
						// std::cerr << "N: " << N << "\n";
						size_t t = 0, tmp;
						f.seekg(0, std::ios::beg);	
						while (t < N) {
							tmp = N - t;
							if (tmp > bytesNum) {
								tmp = bytesNum;
							}
							t += tmp;
							f.read(buf1, tmp);
							buf1[tmp] = '\0';

							f_matches << buf1;
						}
						f.close();
					}
				} else 
				if (splitfilenum[k] == 1) {
					std::ifstream f(resultfile[k].second, std::ios::ate | std::ios::binary);
					if (f.fail()) {
						fprintf(stderr, "error\n");
						exit(0);
					}
					size_t N = f.tellg();
					// std::cerr << "N: " << N << "\n";
					size_t t = 0, tmp;
					f.seekg(0, std::ios::beg);	
					while (t < N) {
						tmp = N - t;
						if (tmp > bytesNum) {
							tmp = bytesNum;
						}
						t += tmp;
						f.read(buf1, tmp);
						buf1[tmp] = '\0';

						f_matches << buf1;
					}
					f.close();
				}
			} else 
			if (smallfile.find(fn) != smallfile.end()) {
				ptr = smallfile[fn];
				// fprintf(stderr, "%lu %lu\n", ptr.first, ptr.second);

				smfi.seekg(ptr.first, std::ios::beg);
				smfi.read(buf1, ptr.second);
				buf1[ptr.second] = '\0';

				f_matches << buf1;
			}
		}
	}
	smfi.close();
	delete[] buf1;
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		displayHelp(argv[0]);
		return 1;
	}
	initPars();
	initHashFunc();
	getPars(argc, argv);
	displayParams();

	CStopWatch stopwatch;
	stopwatch.start();
	GenomeData rGenome, qGenome;
	
	qGenome = readMultiFasta(Q_FN, 125, true, "Query");
	SequenceVector seq2 = std::get<2>(qGenome);
	// fprintf(stderr, "length of qGenome: %lu\n", std::get<0>(qGenome));
	*v1logger << "Time of read Query MultiFasta = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	// BloomFilter(size_t expectedElemNum, double fpr, unsigned hashNum, unsigned kmerSize) :
	size_t qseqlen = std::get<0>(qGenome);
	size_t qexpectedElemNum = (qseqlen - kmer + 1)/(L - kmer + 1); //
	double prob = 0.01; //fpr
	BloomFilter *qbf, *rqbf;
	if (isRC != yes) {
		qbf = new BloomFilter(qexpectedElemNum, prob, 0, kmer);
	}
	if (isRC != no) {
		rqbf = new BloomFilter(qexpectedElemNum, prob, 0, kmer);
	}

	unsigned int qnumHashes;
	if (isRC != yes) {
		qnumHashes = qbf->getHashNum();
	} else
	if (isRC != no) {
		qnumHashes = rqbf->getHashNum();
	}
	// 2^32/8/1000/1000
	// fprintf(stderr, "qexpectedElemNum: %u\n", qexpectedElemNum);
	// fprintf(stderr, "Q hash number: %d\n", qnumHashes);
	// fprintf(stderr, "FilterSize: %lu\n", qbf->getFilterSize());
	// unsigned int qnumHashes = 5;
	char *qstart = std::get<1>(qGenome);
	#ifdef DEBUG
	fprintf(stderr, "%s\n", qstart);
	#endif

	vector<SEGMENT> seg;
	SequenceVector qseq = std::get<2>(qGenome);
	std::vector<size_t> ptr = std::get<3>(qGenome);
	size_t qseqsize = qseq.size();
	for (int i = 0; i < qseqsize; ++i) {
		// fprintf(stderr, "%lu\n", qseq[i].second);
		seg.push_back(SEGMENT(qseq[i].second, (i<qseqsize-1?qseq[i+1].second:qseqlen), ptr[i], ptr[i+1]));
	}
	sort(seg.begin(), seg.end(), cmp);

	//assume n_threads < qseqsize
	vector<SEGMENT> *thr_seg = new vector<SEGMENT>[n_threads];
	size_t *sum_seg = new size_t[n_threads];

	for (int i = 0; i < n_threads && i < qseqsize; ++i) {
		thr_seg[i].push_back(seg[i]);
		sum_seg[i] = seg[i].e - seg[i].b;
	}
	for (int i = n_threads; i < qseqsize; ++i) {
		int idx = 0;
		size_t min = sum_seg[0];
		for (int j = 1; j < n_threads; ++j) {
			if (sum_seg[j] < min) {
				idx = j;
				min = sum_seg[j];
			}
		}
		thr_seg[idx].push_back(seg[i]);
		sum_seg[idx] += seg[i].e - seg[i].b;
	}

	vector<thread> threadVec0;
	for (int i = 0; i < n_threads; ++i) {
		threadVec0.push_back(std::thread(iterate0, qstart, qnumHashes, L, kmer, qbf, rqbf, thr_seg[i]));
	}
	std::for_each(threadVec0.begin(), threadVec0.end(), [](std::thread & thr) {
		thr.join();
	});
	
	deleteReading(qGenome); 

	// fprintf(stderr, "nthashcnt: %lu\n", nthashcnt);

	*v1logger << "Time of ntHashIteratorSimple = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	BloomFilter *rbf, *rrbf;
	if (isRC != yes) {
		rbf = new BloomFilter(qexpectedElemNum, prob, qnumHashes, kmer);
	}
	if (isRC != no) {
		rrbf = new BloomFilter(qexpectedElemNum, prob, qnumHashes, kmer);
	}

	unsigned int rnumHashes;
	if (isRC != yes) {
		rnumHashes = rbf->getHashNum();
	} else
	if (isRC != no) {
		rnumHashes = rrbf->getHashNum();
	}
	// fprintf(stderr, "R hash number: %d\n", rnumHashes);

	rGenome = readMultiFasta(R_FN, 123, true, "Reference");
	
	// fprintf(stderr, "length of rGenome: %u\n", std::get<0>(rGenome));

	*v1logger << "Time of read Reference MultiFasta = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	char* rstart = std::get<1>(rGenome);

	int arrsize = HASH_SIZE;
	
	mm_idx_t *htidx, *rhtidx; //hash table //?the type (size_t) can be changed by the data size like copMEM.
	if (isRC != yes) { //isRC = no or both
		htidx = mm_idx_init(HASH_BITS);
 		mtx = new mutex[arrsize];
	} 
	if (isRC != no) {
		rhtidx = mm_idx_init(HASH_BITS);
		rmtx = new mutex[arrsize];
	}

	// for multithread
	size_t rseqlen = std::get<0>(rGenome);
	SequenceVector rseq = std::get<2>(rGenome);

	vector<size_t> th_pos; //at most n_threads 
	int seg_num = n_threads;
	th_pos.push_back(0);
	size_t avg = rseqlen/seg_num;
	// rseqlen - avg*n_threads
	for (int i = 1; i < seg_num; ++i) {
		th_pos.push_back(i*avg);
	}
	th_pos.push_back(rseqlen);

	// currSize = 0;

	vector<thread> threadVec;
	for (int i = 0; i < n_threads; ++i) {
		threadVec.push_back(std::thread(iterate1, rstart, &th_pos, qnumHashes, kmer, qbf, rqbf, rbf, rrbf, htidx, rhtidx));
	}
	std::for_each(threadVec.begin(), threadVec.end(), [](std::thread & thr) {
		thr.join();
	});

	// fprintf(stderr, "indexedNum: %lu\n", indexedNum); //for output number of indexed k-mers
	if (isRC != yes) {
		delete qbf;
	}
	if (isRC != no) {
		delete rqbf;
	}

	*v1logger << "Time of finding regions on Reference = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	int j, start_a, start_p, n, n_keys;
	idxhash_t *h;

	if (isRC != yes) {
		for (int i = 0; i < HASH_SIZE; ++i) {
			mm_idx_bucket_t *b = &htidx->B[i];
			if (b->a.n > 0) {
				radix_sort_128x(b->a.a, b->a.a + b->a.n);
				// count and preallocate
				// b->n is the number of minimizers appearing >1 times
				for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
					if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
						++n_keys;
						if (n > 1) b->n += n;
						n = 1;
					} else ++n;
				}
				h = kh_init(idx); 
				kh_resize(idx, h, n_keys);
				b->p = (uint64_t*)calloc(b->n, 8);//uint64_t *p; // position array for kmers appearing >1 times

				// create the hash table
				for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
					if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
						khint_t itr;
						int absent;
						mm128_t *p = &b->a.a[j-1];
						itr = kh_put(idx, h, p->x>>HASH_BITS<<1, &absent);
						assert(absent && j - start_a == n);
						if (n == 1) {
							kh_key(h, itr) |= 1;
							kh_val(h, itr) = p->y;
						} else {
							int k;
							for (k = 0; k < n; ++k)
								b->p[start_p + k] = b->a.a[start_a + k].y;
							kh_val(h, itr) = (uint64_t)start_p<<32 | n;
							start_p += n;
						}
						start_a = j, n = 1;
					} else ++n;
				}
				b->h = h;
				assert(b->n == start_p);

				// deallocate and clear b->a
				free(b->a.a);
				b->a.n = b->a.m = 0, b->a.a = 0;
			}
		}
	}

	if (isRC != no) {
		for (int i = 0; i < HASH_SIZE; ++i) {
			mm_idx_bucket_t *b = &rhtidx->B[i];
			if (b->a.n > 0) {
				radix_sort_128x(b->a.a, b->a.a + b->a.n);
				// count and preallocate
				// b->n is the number of minimizers appearing >1 times
				for (j = 1, n = 1, n_keys = 0, b->n = 0; j <= b->a.n; ++j) {
					if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
						++n_keys;
						if (n > 1) b->n += n;
						n = 1;
					} else ++n;
				}
				h = kh_init(idx); 
				kh_resize(idx, h, n_keys);
				b->p = (uint64_t*)calloc(b->n, 8);//uint64_t *p; // position array for minimizers appearing >1 times

				// create the hash table
				for (j = 1, n = 1, start_a = start_p = 0; j <= b->a.n; ++j) {
					if (j == b->a.n || b->a.a[j].x != b->a.a[j-1].x) {
						khint_t itr;
						int absent;
						mm128_t *p = &b->a.a[j-1];
						itr = kh_put(idx, h, p->x>>HASH_BITS<<1, &absent);
						assert(absent && j - start_a == n);
						if (n == 1) {
							kh_key(h, itr) |= 1;
							kh_val(h, itr) = p->y;
						} else {
							int k;
							for (k = 0; k < n; ++k)
								b->p[start_p + k] = b->a.a[start_a + k].y;
							kh_val(h, itr) = (uint64_t)start_p<<32 | n;
							start_p += n;
						}
						start_a = j, n = 1;
					} else ++n;
				}
				b->h = h;
				assert(b->n == start_p);

				// deallocate and clear b->a
				free(b->a.a);
				b->a.n = b->a.m = 0, b->a.a = 0;
			}
		}
	}

	*v1logger << "Time of indexing = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	f_matches.open(matchesFN);
	fsmo.open(smallMF);

	vector<thread> threadVec2;
	for (int i = 0; i < n_threads; ++i) {
		threadVec2.push_back(std::thread(iterate2, &rGenome, rnumHashes, L, kmer, rbf, rrbf, thr_seg[i], htidx, rhtidx));//, &rbf, thr_seg[i], htarr));
	}
	std::for_each(threadVec2.begin(), threadVec2.end(), [](std::thread & thr) {
		thr.join();
	});

	fsmo.close();
	// delete some memory
	for (int i = 0; i < n_threads; ++i) {
		thr_seg[i].clear();
	}
	delete[] thr_seg;

	if (isRC != yes) {
		delete rbf;
		mm_idx_destroy(htidx);
	}
	if (isRC != no) {
		delete rrbf;
		mm_idx_destroy(rhtidx);
	}
	deleteReading(rGenome);

	*v1logger << "Time of matches = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
// exit(0);
	// fprintf(stderr, "rseq.size(): %lu\n", rseq.size());
	fid = 0;
	for (int i = 0; i < resultfile.size(); ++i) {
		splitfilenum.push_back(0);
	}

	std::vector<thread> threadVec3;
	for (int i = 0; i < n_threads; ++i) {
	// for (int i = 0; i < 1; ++i) {
		threadVec3.push_back(std::thread(splitTempFileThread, &rseq));
	}
	std::for_each(threadVec3.begin(), threadVec3.end(), [](std::thread & thr) {
		thr.join();
	});

	if (filesvec.size() > 0) {
		*v1logger << "Number of splitted files = " << filesvec.size() << std::endl;
		// rseq
		*v1logger << "Time of splitting files = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
	}

	filesvecidx = 0;
	std::vector<thread> threadVec4;
	for (int i = 0; i < n_threads; ++i) {
		threadVec4.push_back(std::thread(sortAndFindSeqInfoThread, &rseq));//, &rbf, thr_seg[i], htarr));
	}
	std::for_each(threadVec4.begin(), threadVec4.end(), [](std::thread & thr) {
		thr.join();
	});

	filesvec.clear();
	// rseq
	*v1logger << "Time of sorting & finding seq info = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	mergeSortedFiles(seq2);
	// mergeSortedFiles_old(seq2);

	f_matches.close();
	*v1logger << "Time of merging sorted files = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	char *cmd = new char[128];
	sprintf(cmd, "rm -rf %s", folder.c_str());
	system(cmd);
	delete[] cmd;

	*v1logger << "Time of deleting temporary folder = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
	// printf("matchNum: %d\n", matchNum);

	return 0;
}
#endif
