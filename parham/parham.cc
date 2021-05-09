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

void idx2ij(size_t idx, size_t n_seq, size_t& i, size_t& j) {
	i = idx / n_seq;
	j = idx % n_seq;
	if (i >= j) {
		i = n_seq - i - 2;
		j = j + i + 1;
	}
}

extern "C" {
void compute_hamming_matrix(
		size_t n_seq, size_t seq_len,
		char* seqss, double* poss_freq,
		const char* out_file) {
	int* dis = new int[n_seq * n_seq];
	double* mmf = new double[n_seq * n_seq];
	string* strs = new string[n_seq * n_seq];
	int n_tasks = (n_seq - 1) * n_seq / 2;
#pragma omp parallel for
	for (size_t idx = 0; idx < n_tasks; ++idx) {
		size_t i, j;
		idx2ij(idx, n_seq, i, j);
		hamming(seq_len, 
				seqss + i * seq_len, 
				seqss + j * seq_len, 
				poss_freq,
				dis + idx,
				mmf + idx,
				strs + idx);
	}
	FILE* ouf = fopen(out_file, "w");
	for (size_t idx = 0; idx < n_tasks; ++idx) {
		size_t i, j;
		idx2ij(idx, n_seq, i, j);
		fprintf(ouf, "%d\t%d\t%d\t%f\t%s\n", i, j, dis[idx], 
				mmf[idx], strs[idx].c_str());
	}
	fclose(ouf);
}
};
