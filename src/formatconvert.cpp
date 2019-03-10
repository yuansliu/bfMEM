# include <iostream>
# include <cstdlib>
# include <fstream>
# include <algorithm>
# include <numeric>
# include <string>
# include <vector>
# include <queue>
# include <cassert>
# include <cstring>
# include <unistd.h>
using namespace std;

struct VecEle {
	size_t a, b, c;
	char *s;

	VecEle() {} 
	// VecEle(char *_s, char *_s1, size_t _a, size_t _b, size_t _c): a(_a), b(_b), c(_c) {
	// 	s = strdup(_s);
	// 	s1 = strdup(_s1);
	// }
	VecEle(char *_s, size_t _a, size_t _b): a(_a), b(_b) {
		s = strdup(_s);
	}
};

ofstream f_matches;
const size_t MATCHES_BUFFER_SIZE = 1ULL << 24;

inline bool arrcmp(size_t* p1, size_t* p2) {
	if (p1[1] < p2[1]) return true; 
	if (p1[1] > p2[1]) return false;
	return p1[0] < p2[0];
}

bool cmp(const VecEle &a, const VecEle &b) {
	return a.b < b.b;
	// if (a.b < b.b) return true;
	// if (a.b >= b.b) 
		// return false;
	// return a.a > b.a;
}

int main(int argc, char const *argv[]) {
	// fprintf(stderr, "%s\n", argv[1]);
	FILE *fin = fopen(argv[1], "r");
	if (NULL == fin) {
		fprintf(stderr, "File Open Error!\n");
		exit(1);
	}

	f_matches.open(argv[2]);
	if (f_matches.fail()) {
		fprintf(stderr, "The output file '%s' can not be created.\n", argv[2]);
		exit(1);
	}

	string mbuf;
	mbuf.reserve(MATCHES_BUFFER_SIZE);

	char *s = new char[1024];
	char *s1 = new char[512];
	vector<VecEle> v;
	v.reserve(1<<25);
	size_t a, b, c;

	size_t num = 0;
	while (NULL != fgets(s, 1024, fin)) {
		if (s[0] == '>') {
			if (num < 5) {
				printf("%s", s); //for debug
				++num;
			} else 
			if ( num <= 10 || num == 100 || num == 1000 || 
				 num == 10000 || num == 100000 || num == 1000000 ||
				 num == 10000000 || num == 100000000LL)
			{	
				++ num;
				printf(".\n");
			}

			if (v.size() > 0) {
				stable_sort(v.begin(), v.end(), cmp);
				for (size_t i = 0; i < v.size(); ++i) {
					mbuf.append(v[i].s);
					delete[] v[i].s;
					if (mbuf.size() > (size_t(MATCHES_BUFFER_SIZE*0.95))) {
						f_matches << mbuf; 
						mbuf.clear();
					}
				}
				v.clear();
				v.reserve(1<<25);
			}
			mbuf.append(s);
			if (mbuf.size() > (size_t(MATCHES_BUFFER_SIZE*0.95))) {
				f_matches << mbuf; 
				mbuf.clear();
			}
		} else {
			sscanf(s, "%s%lu%lu", s1, &a, &b);
			v.push_back(VecEle(s, a, b));
		}
	}
	fclose(fin);
	if (v.size() > 0) {
		stable_sort(v.begin(), v.end(), cmp);
		for (size_t i = 0; i < v.size(); ++i) {
			mbuf.append(v[i].s);
			delete[] v[i].s;
			if (mbuf.size() > (size_t(MATCHES_BUFFER_SIZE*0.95))) {
				f_matches << mbuf; 
				mbuf.clear();
			}
		}
	}
	if (mbuf.size() > 0) {
		f_matches << mbuf; 
	}
	f_matches.close();
	delete[] s;
	delete[] s1;
	printf(" over.\n");
	return 0;
}