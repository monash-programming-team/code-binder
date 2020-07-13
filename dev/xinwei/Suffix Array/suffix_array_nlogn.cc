#include <bits/stdc++.h>

using namespace std;

struct suffix {
    int rank[2];
    int idx;
};

class SuffixArray {
private:
    int n;
    void radixSort(vector<suffix> &S, vector<suffix> &tmp, int rk){
	int d = max(256, n);
	vector<int> c(d + 1, 0);
	for (int i = 0; i < n; i++){
	    c[ S[i].rank[rk] ]++;
	}
	for (int i = 0, sum = 0; i <= d; i++){
	    int t = c[i]; c[i] = sum; sum += t;
	}
	for (int i = 0; i < n; i++) tmp[ c[S[i].rank[rk]]++ ] = S[i];
	S = tmp;
    }

public:
    SuffixArray(int _n, string &str, vector<int> &sarray) : n(_n) {
	sarray.resize(n);
	vector<suffix> S(n), tmp(n);
	for (int i = 0; i < n; i++){
	    S[i].idx = i;
	    S[i].rank[0] = (int) str[i] - CHAR_MIN + 1;
	    S[i].rank[1] = 0;
	}
	radixSort(S, tmp, 0);

	for (int k = 1; k < n; k <<= 1){
	    // rank 2nd pair
	    vector<int> A(n);
	    for (int i = 0; i < n; i++) A[ S[i].idx ] = S[i].rank[0];
	    for (int i = 0; i < n; i++){
		int id = S[i].idx;
		S[i].rank[1] = (id + k < n) ? A[id + k] : 0;
	    }
	    radixSort(S, tmp, 1);
	    radixSort(S, tmp, 0);

	    // re-rank
	    S[0].rank[0] = 1;
	    for (int i = 1; i < n; i++){
		S[i].rank[0] = S[i-1].rank[0];
		if (memcmp(tmp[i].rank, tmp[i-1].rank, 8)) S[i].rank[0]++;
/*
		if (tmp[i].rank[0] != tmp[i-1].rank[0] || tmp[i].rank[1] != tmp[i-1].rank[1])
		    S[i].rank[0]++;
*/

	    }
	    if (S[n-1].rank[0] == n) break;
	}
	for (int i = 0; i < n; i++){
	    sarray[i] = S[i].idx;
	}
    }
};


void build_lcp(int n, string &str, vector<int> &sarray, vector<int> &lcp){
    vector<int> index(n); lcp.resize(n);
    for (int i = 0; i < n; i++) index[sarray[i]] = i;
    int len = 0;
    for (int i = 0; i < n; i++){
	int id = index[i];
	if (id){
	    int j = sarray[id - 1];
	    while (i + len < n && j + len < n && str[i + len] == str[j + len]) len++;
	    lcp[id] = len;
	}  
	len = max(len - 1, 0);
    }
    lcp[0] = 0; // invalid value
}

int main(){
    int tc = 1;
    while(tc--){
	string s; cin >> s;
	int n = (int) s.size();
	vector<int> sarray;
	SuffixArray cc(n, s, sarray);
	for (int i = 0; i < n; i++){
	    printf("%d\n", sarray[i]);
	}
    }

    return 0;
}
