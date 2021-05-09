#include <cstdio>
#include <string>

using namespace std;

const char valid_nuc[] = "ACGT";

inline bool in_valid_nuc(char c) {
	for (int i = 0; i < 4; ++i)
		if (valid_nuc[i] == c) return true;
	return false;
}

inline char toUpper(char c) {
	if (c >= 'a' && c <= 'z') c -= 'a' - 'A';
	return c;
}

void hamming(size_t length, const char *seq1, const char *seq2, const double *poss_freq,
		int *distance, double *max_maf, std::string *diff) {
	*distance = 0; *max_maf = 0; diff->clear();
	for (size_t i = 0; i < length; ++i) {
		char c1 = toUpper(seq1[i]), c2 = toUpper(seq2[i]);
		if (in_valid_nuc(c1) && in_valid_nuc(c2)) {
			if (c1 != c2) {
				(*distance)++;
				if (!diff->empty()) diff->append(",");
				diff->append(to_string(i));
				if (poss_freq[i] > *max_maf)
					*max_maf = poss_freq[i];
			}
		} else {
			if (c1 != c2) {
				if (!diff->empty()) diff->append(",");
				diff->append(to_string(i));
			}
		}
	}
	if (*distance == 0) {
		for (size_t i = 0; i < length; ++i) {
			char c1 = toUpper(seq1[i]), c2 = toUpper(seq2[i]);
			if (c1 != c2 && poss_freq[i] > *max_maf) {
				*max_maf = poss_freq[i];
			}
		}
	}
}

extern "C" {
void compute_hamming_matrix(
		size_t n_seq, size_t seq_len,
		char* seqss, double* poss_freq,
		char* freq_file, const char* out_file) {
	int* dis = new int[n_seq * n_seq];
	double* mmf = new double[n_seq * n_seq];
	string* strs = new string[n_seq * n_seq];
#pragma omp parallel for collapse(2)
	for (size_t i = 0; i < n_seq; ++i) {
		for (size_t j = 0; j < n_seq; ++j) {
			if (j < i) {
				size_t my_idx = (i + 1) * i / 2 + j;
				hamming(seq_len, 
						seqss + i * seq_len, 
						seqss + j * seq_len, 
						poss_freq,
						dis + my_idx,
						mmf + my_idx,
						strs + my_idx);
			}
		}
	}
	FILE* ouf = fopen(out_file, "w");
	for (size_t i = 0; i < n_seq; ++i) {
		for (size_t j = 0; j < n_seq; ++j) {
			if (j < i) {
				size_t my_idx = (i + 1) * i / 2 + j;
				fprintf(ouf, "%d\t%d\t%d\t%f\t%s\n", i, j, dis[my_idx], 
						mmf[my_idx], strs[my_idx]);
			}
		}
	}
	fclose(ouf);
}
};
