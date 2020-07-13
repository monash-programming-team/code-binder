// KMP pattern matching algorithm
//
// Author: Daniel Anderson
// Date: 27-12-2016
// Reliability: 5
// Tested on: SPOJ-NHAY, UVA10298, UVA455, UVA11362, SPOJ-PERIOD
//
// Usage: Compute the prefix array of your pattern then
// run KMP.
//
//   vi pre = prefix(pat);
//   vi res = find_pattern(str, pat, pre);
//
// Returns the indices of all locations i such that
// str[i...i+m) == pat. If you only need one location,
// modify kmp so that it returns after finding one.
//
// Complexity: O(n + m) where n = |str|, m = |pat|
#include<bits/stdc++.h>
using namespace std;

typedef vector<int> vi;

//listings:kmp
// Compute the prefix array for the pattern pat. The prefix array is for
// each index i, the length of the longest proper suffix of pat[0...i] that
// is also a proper prefix of pat[0...i]. Complexity: O(m)
template<typename T> vi prefix(const T& pat) {
  int m = (int)pat.size();  vi pre(m, 0);
  for (int j=0, i=1; i<m; ) {
    if (pat[i] == pat[j]) pre[i++] = ++j;
    else if (j>0) j = pre[j-1];
    else i++;
  }
  return pre;
}

// Knuth-Morris-Pratt pattern matching. Complexity: O(n)
// Find all occurrences of the pattern pat in the string str using
// the prefix array pre computed by prefix(pat).
template<typename T> vi find_pattern(const T& str, const T& pat, const vi& pre) {
  int n = (int)str.size(), m = (int)pat.size();  vi res;
  for (int i=0, j=0; i<n; i++) {
    while (j > 0 && str[i] != pat[j]) j = pre[j-1];
    if (str[i] == pat[j]) j++;
    if (j == m) res.push_back(i - m + 1), j = pre[j-1];
  }
  return res;
}
//listings:/kmp

// Verdict: AC
void solve_SPOJ_NHAY() {
  int n; bool first = true;
  while (cin >> n) {
    if (!first) cout << '\n'; first = false;
    string pat, str;
    cin >> pat >> str;
    auto pre = prefix(pat);
    auto loc = find_pattern(str, pat, pre);
    for (int x : loc) cout << x << '\n';
  }
}

// Verdict: AC
void solve_UVA10298() {
  string S;
  while (cin >> S, S != ".") {
    vi pre = prefix(S);
    auto loc = find_pattern(S + S, S, pre);
    cout << S.size() / loc[1] << '\n';
  }
}

// Verdict: AC
void solve_UVA455() {
  int N; cin >> N;
  while (N--) {
    string S; cin >> S;
    int n = (int)S.size();
    for (int l = 1; l <= n; l++) {
      string pat = S.substr(0, l);
      vi pre = prefix(pat);
      auto loc = find_pattern(S, pat, pre);
      bool good = true;
      for (int i = 0; i < n; i += l)
        if (!binary_search(loc.begin(), loc.end(), i)) good = false;
      if (good) { cout << l << '\n'; if (N) cout << '\n'; break; }
    }
  }
}

// Verdict: AC
void solve_UVA11362() {
  int t; cin >> t;
  while (t--) {
    int n; cin >> n;
    vector<string> pn(n);
    for (auto& x : pn) cin >> x;
    sort(pn.begin(), pn.end());
    bool good = true;
    for (int i = 0; i < n - 1 && good; i++) {
      vi pre = prefix(pn[i]);
      auto loc = find_pattern(pn[i+1], pn[i], pre);
      if (!loc.empty() && loc.front() == 0) good = false;
    }
    if (good) cout << "YES" << '\n';
    else cout << "NO" << '\n';
  }
}

// Verdict: AC
void solve_SPOJ_PERIOD() {
  int T; cin >> T;
  for (int t=1; t<=T; t++) {
    cout << "Test case #" << t << '\n';
    int N; string S; cin >> N >> S;
    vi pre = prefix(S);
    for (int i=1; i<N; i++) {
      if (pre[i] >= (i + 2) / 2 && pre[i] % (i + 1 - pre[i]) == 0) {
        int K = (i + 1) / (i + 1 - pre[i]);
        cout << (i + 1) << ' ' << K << '\n';
      }      
    }
    cout << '\n';
  }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  //solve_SPOJ_NHAY();
  //solve_UVA10298();
  //solve_UVA455();
  //solve_UVA11362();
  solve_SPOJ_PERIOD();
}
