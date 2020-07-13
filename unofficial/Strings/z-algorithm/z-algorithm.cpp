// The Z-algorithm and its application to pattern matching
//
// Reliability: 5
// Tested on: UVA455, UVA10298, UVA11362, SPOJ-NHAY, UVA11475, UVA12467
//
//
// Complexity: O(n)
#include<bits/stdc++.h>
using namespace std;

typedef vector<int> vi;

//listings:z
// The z-array of the sequence s is for each index i, the length
// of the longest substring beginning at i that is also a prefix
// of s. Complexity: O(n)
template<typename T> vi z_array(const T& s) {
  int n = (int)s.size(), L = 0, R = 0;  vi z(n, n - 1);
  for (int i = 1, j; i < n; i++) {
    j = max(min(z[i-L],R-i),0);
    for (; i + j < n && s[i+j] == s[j]; j++);
    z[i] = j;
    if (i + z[i] > R) R = i + z[i], L = i;
  }
  return z;
}

// Find all occurrences of pat in str using the z-algorithm in O(n + m)
template<typename T> vi find_pattern(const T& str, T pat) {
  int n = (int)str.size(), m = (int)pat.size();
  pat.insert(pat.end(), str.begin(), str.end());
  vi z = z_array(pat), res;
  for (int i = 0; i < n; i++) if (z[i + m] >= m) res.push_back(i);
  return res;
}
//listings:/z

// Verdict: AC
void solve_UVA455() {
  int N; cin >> N;
  while (N--) {
    string S; cin >> S;
    int n = (int)S.size();
    for (int l = 1; l <= n; l++) {
      string pat = S.substr(0, l);
      auto loc = find_pattern(S, pat);
      bool good = true;
      for (int i = 0; i < n; i += l)
        if (!binary_search(loc.begin(), loc.end(), i)) good = false;
      if (good) { cout << l << '\n'; if (N) cout << '\n'; break; }
    }
  }
}

// Verdict: AC
void solve_UVA10298() {
  string S;
  while (cin >> S, S != ".") {
    auto loc = find_pattern(S + S, S);
    cout << S.size() / loc[1] << '\n';
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
      auto loc = find_pattern(pn[i+1], pn[i]);
      if (!loc.empty() && loc.front() == 0) good = false;
    }
    if (good) cout << "YES" << '\n';
    else cout << "NO" << '\n';
  }
}

// Verdict: AC
void solve_SPOJ_NHAY() {
  int n; bool first = true;
  while (cin >> n) {
    if (!first) cout << '\n'; first = false;
    string pat, str;
    cin >> pat >> str;
    auto loc = find_pattern(str, pat);
    for (int x : loc) cout << x << '\n';
  }
}

// Verdict: AC
void solve_UVA11475() {
  string s, t, add;
  while (cin >> s) {
    t = s;
    reverse(t.begin(), t.end());
    vi z = z_array(t + s);
    int n = (int)s.size(), best = -1;
    for (int i = 0; i < n; i++) if (z[i + n] == n - i) { best = i; break; }
    if (best == -1) add = s.substr(0, n - 1);
    else add = s.substr(0, best);
    reverse(add.begin(), add.end());
    s += add;
    cout << s << '\n';
  }
}

// Verdict: AC
void solve_UVA12467() {
  int T; cin >> T;
  while (T--) {
    string S, R; cin >> S;
    int n = (int)S.size();
    R = S; reverse(R.begin(), R.end());
    vi z = z_array(S + R);
    int best = 0;
    for (int i = 1; i < n; i++) if (z[n + i] > z[n + best]) best = i;
    string ans = R.substr(best, z[n + best]);
    reverse(ans.begin(), ans.end());
    cout << ans << '\n';
  }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  solve_SPOJ_NHAY();
}
