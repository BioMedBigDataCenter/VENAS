#include <cstdio>
#include <string>
#include <vector>
#include <set>
#include <algorithm>

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

int *dis, *dsr;
int* g_posref;
double* mmf;
string* strs;

inline bool cmp_posref(const int& a, const int& b) {
	return g_posref[a] < g_posref[b];
}

size_t g_n_seq;

inline bool cmp_idx(const int& a, const int& b) {
	if (dis[a] != dis[b]) {
		return dis[a] < dis[b];
	} else if (mmf[a] != mmf[b]) {
		return mmf[a] < mmf[b];
	} else {
		// For total numerical equivalence to python's sort.
		size_t ia, ja, ib, jb;
		idx2ij(a, g_n_seq, ia, ja);
		idx2ij(b, g_n_seq, ib, jb);
		if (ia == ib) {
			return ja < jb;
		} else {
			return ia < ib;
		}
	}
}

int getroot(int u) {
	int p, q;
	for (p = dsr[u]; p != dsr[p]; p = dsr[p]);
	for (; u != p; q = dsr[u], dsr[u] = p, u = q);
	return u;
}

extern "C" {
void compute_hamming_matrix(
		size_t n_seq, size_t seq_len,
		char* seqss, double* poss_freq, int* poss_ref,
		const char* out_file, const char* net_file) {
	g_posref = poss_ref;
	int n_tasks = (n_seq - 1) * n_seq / 2;
	dis = new int[n_tasks];
	mmf = new double[n_tasks];
	strs = new string[n_tasks];

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
	if (out_file) {
		FILE* ouf = fopen(out_file, "w");
		for (size_t idx = 0; idx < n_tasks; ++idx) {
			size_t i, j;
			idx2ij(idx, n_seq, i, j);
			fprintf(ouf, "%d\t%d\t%d\t%f\t%s\n", i, j, dis[idx], 
					mmf[idx], strs[idx].c_str());
		}
		fclose(ouf);
	}
	if (net_file) {
        puts("Generate a list of candidate links");
		int* idxs = new int[n_tasks];
		size_t n_clades = 0;
		for (int i = 0; i < n_tasks; ++i) {
			idxs[i] = i;
		}
		g_n_seq = n_seq;
		sort(idxs, idxs + n_tasks, cmp_idx);
        puts("Candidate link ranking");
		vector<int> net_list;
		set<int> added_nodes;
		dsr = new int[n_seq];
		for (int i = 0; i < n_seq; ++i) {
			dsr[i] = i;
		}
		for (size_t iidx = 0; iidx < n_tasks; ++iidx) {
			size_t idx = idxs[iidx];
			size_t i, j;
			idx2ij(idx, n_seq, i, j);
			int in_i = (added_nodes.find(i) != added_nodes.end());
			int in_j = (added_nodes.find(j) != added_nodes.end());
			added_nodes.insert(in_i);
			added_nodes.insert(in_j);
			if (!in_i || !in_j) {
				n_clades += 2 - in_i - in_j;
			}
			int ri = getroot(i);
			int rj = getroot(j);
			if (ri != rj) {
				dsr[ri] = rj;
				net_list.push_back(idx);
			}
			if (added_nodes.size() == n_seq && n_clades < 2) {
				break;
			}
		}
		printf("Net list size %lu, seq n %lu\n", net_list.size(), n_seq);
		FILE* ouf = fopen(net_file, "w");
		for (auto idx : net_list) {
			size_t i, j;
			idx2ij(idx, n_seq, i, j);
			vector<int> diff;
			int difv = -1;
			size_t sl = strs[idx].length();
			for (size_t p = 0; p < sl; ++p) {
				if (isdigit(strs[idx][p])) {
					if (difv == -1) {
						difv = 0;
					}
					difv = difv * 10 + strs[idx][p] - 48;
				} else {
					diff.push_back(difv);
					difv = -1;
				}
			}
			if (difv != -1) {
				diff.push_back(difv);
			}
			sort(diff.begin(), diff.end(), cmp_posref);
			fprintf(ouf, "%lu\t%lu\t%lu\t", i + 1, j + 1, dis[idx]);
			bool is_first = true;
			for (auto posi : diff) {
				if (is_first) {
					is_first = false;
				} else {
					fprintf(ouf, "|");
				}
				fprintf(ouf, "%d:%c:%c", poss_ref[posi],
						seqss[i * seq_len + posi],
						seqss[j * seq_len + posi]);
			}
			fprintf(ouf, "\n");
		}
		fclose(ouf);
		delete[] idxs;
		delete[] dsr;
	}
	delete[] dis;
	delete[] mmf;
	delete[] strs;
}
};
