// Longest Increasing Subsequence
//
// Author : Daniel Anderson
// Date : 03/12/2016
// Reliability: 5
// Tested On: SPOJ - ELIS, Codeforces 269B, UVA481, UVA11456, UVA1062
//
#include<bits/stdc++.h>
using namespace std;

typedef vector<int> vi;

//listings:lis
// Returns the length of the longest (strictly) increasing subsequence of v
// Reverse the input for longest decreasing subsequence. Complexity: O(n log(n))
template<typename T> int lis_len(vector<T>& v) {
  vector<T> s(v.size()); int k=0;
  for (int i=0; i<(int)v.size(); i++) {  // Change to upper_bound for non-decreasing
  	auto it = lower_bound(s.begin(), s.begin()+k, v[i]); *it = v[i];
  	if (it == s.begin()+k) k++;
  }
  return k;
}

// Returns the longest (strictly) increasing subsequence of v
// Reverse the input for longest decreasing subsequence. Complexity: O(n log(n))
template<typename T> vector<T> lis(vector<T>& v) {
  int n = v.size(), len = 0;  vi tail(n), prev(n); T val[n];
  for (int i=0; i < n; i++) { // Change to upper_bound for non-decreasing
    int pos = lower_bound(val, val + len, v[i]) - val;
    len = max(len, pos + 1), prev[i] = (pos > 0 ? tail[pos - 1] : -1);
    tail[pos] = i;  val[pos] = v[i];
  }
  vector<T> res(len);
  for (int i = tail[len - 1]; i >= 0; i = prev[i]) res[--len] = v[i];
  return res;
}
//listings:/lis

// The non-decreasing version.
template<typename T> int lnds_len(const vector<T>& v) {
  multiset<T> s;
  for (const auto i : v) {
    s.insert(i);
    auto it = s.upper_bound(i);
    if (it != s.end()) s.erase(it);
  }
  return s.size();
}

// Verdict: AC
void solve_UVA1062() {
  string in;
	int no = 1;
	while (cin >> in, in != "end") {
		vi a(in.begin(), in.end());
		cout << "Case " << no++ << ": " << lis_len(a) << '\n';
	}
}

// Verdict: AC
void solve_UVA481() {
  int x; vi a;
  while (cin >> x) a.push_back(x);
  cout << lis_len(a) << '\n' << '-' << '\n';
  for (auto x : lis(a)) cout << x << '\n';
}

// Verdict: AC
void solve_UVA11456() {
  int t; cin >> t;
  while (t--) {
    int n; cin >> n;
    vi a(n); int ans = 0;
    for (int i=0; i<n; i++) cin >> a[i];
    for (int mid=0; mid<n; mid++) {
      vi bigger, smaller;
      for (int j=mid+1; j<n; j++) {
        if (a[j] > a[mid]) bigger.push_back(a[j]);
        if (a[j] < a[mid]) smaller.push_back(a[j]);
      }
      reverse(smaller.begin(), smaller.end());
      ans = max(ans, lis_len(bigger) + lis_len(smaller) + 1);
    }
    cout << ans << '\n';
  }
}

// Verdict: AC
void solve_spoj_elis() {
  int N; cin >> N;
  vi A(N);
  for (int i=0; i<N; i++) cin >> A[i];
  cout << lis_len(A) << endl;
}

// Verdict: AC
void solve_CF269B() {
  int n, m; cin >> n >> m;
  vi s(n); double x;
  for (int i=0; i<n; i++) cin >> s[i] >> x;
  int len = lnds_len(s);
  cout << n - len << endl;
}

int main() {
  //solve_UVA1062();
  //solve_UVA11456();
  //solve_UVA481();
  //solve_CF269B();
  solve_spoj_elis();
  
  // Test on non-numeric input
  //vector<string> a = {"Daniel", "Xin Wei", "Peter"};
  //int len = lis_len(a); auto res = lis(a);
  //assert(len == (int)res.size());
  //for (auto& x : res) cout << x << endl;
}
