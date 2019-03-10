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

int main(int argc, char const *argv[]) {
	// fprintf(stderr, "%s\n", argv[1]);
	FILE *fa = fopen(argv[1], "r");
	if (NULL == fa) {
		fprintf(stderr, "The file '%s' can not be created.\n", argv[1]);
		exit(1);
	}

	FILE *fb = fopen(argv[2], "r");
	if (NULL == fb) {
		fprintf(stderr, "The file '%s' can not be created.\n", argv[2]);
		exit(1);
	}

	char *sa = new char[1024];
	char *sb = new char[1024];
	char *sa1 = new char[512];
	char *sb1 = new char[512];
	size_t a1, b1, c1;
	size_t a2, b2, c2;

	bool flag = true;
	size_t linenum = 0;
	while (NULL != fgets(sa, 1024, fa)) {
		++ linenum;
		if (NULL == fgets(sb, 1024, fb)) {
			fprintf(stderr, "The file '%s' contains more lines than the file '%s'\n", argv[1], argv[2]);
			flag = false;
			break;
		}
		if (sa[0] == '>') {
			if (strcmp(sa, sb) != 0) {
				flag = false;
				break;
			}
		} else {
			sscanf(sa, "%s%lu%lu%lu", sa1, &a1, &b1, &c1);
			sscanf(sb, "%s%lu%lu%lu", sb1, &a2, &b2, &c2);
			if (strcmp(sa1, sb1) != 0 || a1 != a2 || b1 != b2 || c1 != c2) {
				fprintf(stderr, "line: %lu\n%s%s", linenum, sa, sb);
				flag = false;
				break;
			}
		}
	}
	if (flag && NULL != fgets(sb, 1024, fb)) { // sb not to end
		if (strlen(sb) > 4) {
			fprintf(stderr, "The file '%s' contains more lines than the file '%s'\n", argv[2], argv[1]);
			flag = false;
		}
	}
	fclose(fa);
	fclose(fb);
	delete[] sa;
	delete[] sb;
	delete[] sa1;
	delete[] sb1;

	if (flag) {
		fprintf(stderr, "MEMs in the two file are the same.\n");
	} else {
		fprintf(stderr, "Two files are different.\n");
	}
	return 0;
}