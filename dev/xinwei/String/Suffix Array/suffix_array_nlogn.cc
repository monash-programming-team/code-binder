// 'sarray'
// tested on: South Pacific Divisions 2016 - Intuidiff
//            RMI2016 - Frequency
//
// 'lcp'
// tested on: South Pacific Divisions 2016 - Intuidiff
//            RMI2016 - Frequency
//
// 'lce' tested on: South Pacific Divisions 2016 - Intuidiff
// 'find' tested on: South Pacific Divisions 2016 - Intuidiff
// 'max_match' tested on: South Pacific Divisions 2016 - Intuidiff
#include <bits/stdc++.h>

using namespace std;
#define X first
#define Y second
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<vi> vvi;

#define debug(x) cerr << #x << " = " << (x) << endl;
template<typename T>
ostream& operator<<(ostream& o, vector<T>& v) {
    for (auto& x : v) o << x << ' ';
    return o;
}

struct suffix {
    int rank[2], idx;
};

class SuffixArray {
private:
    int n;
    string str;
    vi sarray, indx;
    vi lcp; // OPTIONAL: only if lcp array is needed

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
    SuffixArray(string &_str) : str(_str) {}
    vi build_sarray(){
        n = (int)str.size();
	sarray.resize(n); indx.resize(n);
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
		if (tmp[i].rank[0] != tmp[i-1].rank[0] ||
                    tmp[i].rank[1] != tmp[i-1].rank[1]) S[i].rank[0]++;
	    }
	    if (S[n-1].rank[0] == n) break;
	}
	for (int i = 0; i < n; i++){
	    sarray[i] = S[i].idx;
            indx[S[i].idx] = i;
	}
        return sarray;
    }
    vi build_lcp(){
        lcp.resize(n);
        int len = 0;
        for (int i = 0; i < n; i++){
            int id = indx[i];
            if (id){
                int j = sarray[id - 1];
                while (i + len < n && j + len < n
                       && str[i + len] == str[j + len]) len++;
                lcp[id] = len;
            }
            len = max(len - 1, 0);
        }
        lcp[0] = 0; // invalid value
        return lcp;
    }
    vi get_indx() { return indx; }
    struct Comp {
        const string &s; int m, j;
        Comp(const string &str, int M, int _j=0) : s(str), m(M), j(_j) {}
        bool operator() (int i, const string &p) const { return s.compare(i,m,p,j,m) < 0; }
        bool operator() (const string &p, int i) const { return s.compare(i,m,p,j,m) > 0; }
    };
    pii find(const string &pat, int j=0){
        auto p = equal_range(sarray.begin(), sarray.end(), pat, Comp(str, pat.size(), j));
        return {p.X - sarray.begin(), p.Y - sarray.begin()};
    }
    // returns max substring in 'str' matching 'pat'
    int max_match(const string &pat, int j=0){
        int m = (int)pat.size() - j;
        pii p = find(pat, j);
        if (p.X != p.Y) return m;
        int res = 0;
        if (p.X) {
            int id = sarray[p.X - 1];
            int i;
            for (i = 0; i < min(m, n - id); i++){
                if (pat[j + i] != str[id + i]) break;
            }
            res = max(res, i);
        }
        if (p.X != n){
            int id = sarray[p.X];
            int i;
            for (i = 0; i < min(m, n - id); i++){
                if (pat[j + i] != str[id + i]) break;
            }
            res = max(res, i);
        }
        return res;
    }
    // OPTIONAL: only if lce(i, j) is needed
    /*
    RMQ<int> rmq;
    void init_lce() { rmq = RMQ<int>(lcp); }
    // Longest Common Extension
    // returns longest common substring beginning from i & j
    int lce(int i, int j) {
        if (i == j) return n - i;
        return rmq.query(i + 1, j).first;
    }
    */
};

int main(){
    std::ios_base::sync_with_stdio(false); cin.tie(0);
    string str = "BANANA";
    SuffixArray sa(str);
    vi sarray = sa.build_sarray();
    string patt = "AB";
    cout << sa.max_match(patt) << endl;

    patt = "#";
    cout << sa.max_match(patt) << endl;
}
