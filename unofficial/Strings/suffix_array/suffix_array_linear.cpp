// A linear time suffix array construction with LCP.
//
// Author: Daniel (modified from Darcy's code binder)
// Date: 04/01/2017
// Reliability: 5
// Tested on: Tuesday Tutorial - Suffix Array, SPOJ - SARRAY, SPOJ - DISUBSTR,
//   SPOJ - SUBST1, SPOJ - MINMOVE, SPOJ - LCS, SPOJ - NHAY
//
// Complexity: O(n) to build suffix array and compute LCP.
//             O(m log(n)) to do pattern matching.
#include <bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef pair<int,int> pii;
typedef vector<int> vi;

//listings:suffix_array_linear
// Suffix array construction with LCP. The suffix array is built into sarray and the
// LCP into lcp. NOTE: sarray does not include the empty suffix. lcp[i] is the longest
// common prefix between the strings at sarray[i-1] and sarray[i], lcp[0] = 0.
// Complexity: O(N) or O(N log(N)) for suffix array. O(N) for LCP.
struct suffix_array {
  int n; string str; vi sarray, lcp;
  void bucket(vi& a, vi& b, vi& r, int n, int K, int off=0) {
    vi c(K+1, 0);
    for (int i=0; i<n; i++) c[r[a[i]+off]]++;
    for (int i=0, sum=0; i<=K; i++) { int t = c[i]; c[i] = sum; sum += t; }
    for (int i=0; i<n; i++) b[c[r[a[i]+off]]++] = a[i];
  }
  // Create the suffix array and LCP array of the string s. (LCP is optional)
  suffix_array(string s) : n(s.size()), str(move(s)) { build_sarray(); build_lcp(); }
  // ------------------ OPTION 1: Linear time suffix array ----------------------------
  #define GetI() (SA12[t] < n0 ? SA12[t] * 3 + 1 : (SA12[t] - n0) * 3 + 2)
  typedef tuple<int,int,int> tiii;
  void sarray_int(vi &s, vi &SA, int n, int K) {
    int n0=(n+2)/3, n1=(n+1)/3, n2=n/3, n02=n0+n2, name=0, c0=-1, c1=-1, c2=-1;
    vi s12(n02 + 3, 0), SA12(n02 + 3, 0), s0(n0), SA0(n0);
    for (int i=0, j=0; i < n+(n0-n1); i++) if (i%3 != 0) s12[j++] = i;
    bucket(s12, SA12, s, n02, K, 2), bucket(SA12, s12, s, n02, K, 1);
    bucket(s12, SA12, s, n02, K, 0);
    for (int i = 0; i < n02; i++) {
      if (s[SA12[i]] != c0 || s[SA12[i]+1] != c1 || s[SA12[i]+2] != c2) 
        name++, c0 = s[SA12[i]], c1 = s[SA12[i]+1], c2 = s[SA12[i]+2];
      if (SA12[i] % 3 == 1) s12[SA12[i]/3] = name;
      else s12[SA12[i]/3 + n0] = name;
    }
    if (name < n02) {
      sarray_int(s12, SA12, n02, name);
      for (int i = 0; i < n02; i++) s12[SA12[i]] = i + 1;
    } else for (int i = 0; i < n02; i++) SA12[s12[i] - 1] = i;
    for (int i=0, j=0; i < n02; i++) if (SA12[i] < n0) s0[j++] = 3*SA12[i];
    bucket(s0, SA0, s, n0, K);
    for (int p=0, t=n0-n1, k=0; k < n; k++) {
      int i = GetI(), j = SA0[p];
      if (SA12[t] < n0 ?
          (pii(s[i], s12[SA12[t] + n0]) < pii(s[j], s12[j/3])) :
          (tiii(s[i],s[i+1],s12[SA12[t]-n0+1]) < tiii(s[j],s[j+1],s12[j/3+n0]))) {
        SA[k] = i; t++;
        if (t == n02) for (k++; p < n0; p++, k++) SA[k] = SA0[p];
      } else {
        SA[k] = j; p++;
        if (p == n0) for (k++; t < n02; t++, k++) SA[k] = GetI();
      }
    }
  }
  void build_sarray() {
    if (n <= 1) { sarray.assign(n, 0); return; }
    vi s(n+3, 0); sarray.assign(n+3, 0);
    for (int i=0; i<n; i++) s[i] = (int)str[i] - CHAR_MIN + 1;
    sarray_int(s, sarray, n, 256), sarray.resize(n);
  }
  // ------------------ OPTION 2: O(N log(N)) time suffix array -----------------------
  void build_sarray() {
    sarray.assign(n, 0); vi r(2*n, 0), sa(2*n), tmp(2*n); if (n <= 1) return;
    for (int i=0; i<n; i++) r[i] = (int)str[i] - CHAR_MIN + 1, sa[i] = i;
    for (int k=1; k<n; k *= 2) {
      bucket(sa,tmp,r,n,max(n,256),k), bucket(tmp,sa,r,n,max(n,256),0);
      tmp[sa[0]] = 1;
      for (int i=1; i<n; i++) {
        tmp[sa[i]] = tmp[sa[i-1]];
        if ((r[sa[i]] != r[sa[i-1]]) || (r[sa[i]+k] != r[sa[i-1]+k])) tmp[sa[i]]++;
      }
      copy(tmp.begin(), tmp.begin()+n, r.begin());
    }
    copy(sa.begin(), sa.begin()+n, sarray.begin());
  }
  // --------------------- OPTIONAL: If you need LCP array ----------------------------
  void build_lcp() {  
    int h = 0; vi rank(n); lcp.assign(n, 0);
    for (int i = 0; i < n; i++) rank[sarray[i]] = i;
    for (int i = 0; i < n; i++) {
      if (rank[i] > 0) {
        int j = sarray[rank[i]-1];
        while (i + h < n && j + h < n && str[i+h] == str[j+h]) h++;
        lcp[rank[i]] = h;
      }
      if (h > 0) h--;
    }
  } 
  // OPTIONAL: Pattern matching -- Find all occurrences of pat[j..] in O(m log(n))
  // Returns an iterator pair of the matching locations in the suffix array
  struct Comp {
    const string& s; int m, j;
    Comp(const string& str,int m, int j) : s(str), m(m), j(j) { }
    bool operator()(int i, const string& p) const { return s.compare(i,m,p,j,m) < 0; }
    bool operator()(const string& p, int i) const { return s.compare(i,m,p,j,m) > 0; }
  };
  auto find(const string& pat, int j=0) { 
    return equal_range(sarray.begin(), sarray.end(), pat, Comp(str,pat.size(),j));
  }
};
//listings:/suffix_array_linear

/*
// OPTIONAL: Partial pattern matching -- Find the longest prefix of pat[j..] that
// occurs in str. Returns max length and an iterator pair of the matching locations.
struct PosComp {
  const string& s; int n, off;
  PosComp(const string& str, int off=0) : s(str), n(s.size()), off(off) { }
  bool operator()(int i, char c) { return i+off<n ? s[i+off]<c : 1; }
  bool operator()(char c, int i) { return i+off<n ? c<s[i+off] : 0; }
};
auto longest_match(const string& pat, int pos) {
  int i = 0; auto cmp = PosComp(str); auto lo = sarray.begin(), hi = sarray.end(); 
  while (pos + i < (int)pat.size()) {
    auto res = equal_range(lo, hi, pat[pos+i], cmp);
    if (res.X != res.Y) lo = res.X, hi = res.Y, i++, cmp.off++; else break;
  }
  return make_pair(i, make_pair(lo, hi));
}
*/

// Tuesday tutorial: Suffix Array
// Verdict: AC
// Print the suffix array
void solve_TT_A() {
  string s; cin >> s;
  suffix_array sa(s);
	for (int x : sa.sarray) cout << x << ' '; cout << '\n';
}
 
// Verdict: AC (100 points)
// Print the suffix array
void solve_SPOJ_sarray() {
  string s; cin >> s;
  suffix_array sa(s);
	for (int x : sa.sarray) cout << x << '\n';
}

// Verdict: AC
// NOTE: Also solves SPOJ - SUBST1 (Same I/O format)
// Count the number of distinct substrings
void solve_SPOJ_DISUBSTR() {
  int T; cin >> T;
  while (T--) {
    string s; cin >> s;
    ll n = (ll)s.size();
    suffix_array sa(s);
    ll ans = n * (n + 1) / 2;
    for (int i=0; i<n; i++) ans -= sa.lcp[i];
    cout << ans << '\n';
  }
}

// Verdict: AC
// Find the minimum number of rotations to achieve the
// lexicographically least rotation
void solve_SPOJ_MINMOVE() {
  string S; cin >> S; int n = (int)S.size();
  suffix_array sa(S+S);
  auto i = distance(sa.sarray.begin(),
    find_if(sa.sarray.begin(), sa.sarray.end(), [=](int x) { return x < n; }));
  int ans = sa.sarray[i];
  while (i < 2*n && sa.lcp[++i] >= n) ans = min(ans, sa.sarray[i]);
  cout << ans << endl;
}

// Verdict: AC
// Find the longest common substring of two strings
// Requires RMQ
template<typename T> struct RMQ {
	vector<vector<pair<T, int>>> sptable;  vi lg;
	RMQ(vector<T> &v)  {
		int n = (int) v.size();
		lg.assign(n + 1, 0);
		for (int i = 2; i <= n; i++) lg[i] = lg[i/2] + 1;
		sptable.assign(lg[n] + 1, vector<pair<T, int>>(n));
		for (int i = 0; i < n; i++) sptable[0][i] = {v[i], i};
		for (int i = 1; i <= lg[n]; i++) for (int j = 0; j + (1 << i) - 1 < n; j++)
			sptable[i][j] = min(sptable[i-1][j], sptable[i-1][j + (1 << (i-1))]);
	}
	// Find the minimum element and its position in the range [L, R]
	pair<T, int> query(int L, int R) {
		int k = lg[R - L + 1];
		return min(sptable[k][L], sptable[k][R - (1 << k) + 1]);
	}
};
void solve_SPOJ_LCS() {
  string s, t; cin >> s >> t;
  int n = (int)s.size(), m = (int)t.size();
  suffix_array sa(s + "$" + t);
  RMQ<int> rmq(sa.lcp);
  vi lhs, rhs; int ans = 0;
  for (int i=0; i<n+m+1; i++)
    if (sa.sarray[i] < n) lhs.push_back(i);
    else if (sa.sarray[i] > n) rhs.push_back(i);
  for (int i : lhs) {
    auto it = lower_bound(rhs.begin(), rhs.end(), i);
    if (it != rhs.end()) ans = max(ans, rmq.query(i+1, *it).X);
    if (it != rhs.begin()) ans = max(ans, rmq.query(*(it-1)+1, i).X);
  }
  cout << ans << endl;
}

// Verdict: AC
// Find all occurrences of a pattern in a string
void solve_SPOJ_NHAY() {
  int n; bool first = true;
  while (cin >> n) {
    if (!first) cout << '\n'; first = false;
    string pat, str;
    cin >> pat >> str;
    suffix_array sa(str);
    auto res = sa.find(pat);
    vi ans;
    for (auto it = res.X; it != res.Y; it++) ans.push_back(*it);
    sort(ans.begin(), ans.end());
    for (int x : ans) cout << x << '\n';
  }
}

int main() {
  ios_base::sync_with_stdio(0); cin.tie(0);
  //solve_TT_A();
  solve_SPOJ_sarray();
  //solve_SPOJ_DISUBSTR();
  //solve_SPOJ_MINMOVE();
  //solve_SPOJ_LCS();
  //solve_SPOJ_NHAY();
}

