// Iterative segment tree for RMQ
//
// Author: Daniel Anderson
// Date: 20-01-2016
// Reliability: 5
// Tested on: SPOJ-RMQSQ, SPOJ-RPLN, CF100971FM, SPOJ-FREQUENT, Randomized cases
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:rmq
// Segment tree for dynamic range minimum query. For maximum query, change min
// to max, and use min() as the identity value. Can also use gcd, lcm, sum etc
// with appropriate identity. Complexity: O(n) to build, O(log(n)) to query.
template<typename T> struct SegmentTree {
  int n; vector<pair<T,int>> st; const pair<T,int> I = {numeric_limits<T>::max(),-1};
  SegmentTree(const vector<T>& A) : n(A.size()), st(2*n, I) {
    for (int i=0;i<n;i++) st[n+i] = {A[i],i};
    for (int i=n-1; i; --i) st[i] = min(st[2*i], st[2*i+1]);
  }
  void update(int i, int val) {         // Set A[i] = val
    for (st[i+=n] = {val,i}; i > 1; i/= 2) st[i/2] = min(st[i], st[i^1]);
  }
  pair<T,int> query(int l, int r) {     // Find min A[l..r]
    pair<T,int> res = I;
    for (l += n, r += n; l <= r; l /= 2, r /= 2) {
      if (l&1) res = min(res, st[l++]);
      if (~r&1) res = min(res, st[r--]);
    }
    return res;
  }
};
//listings:/rmq

// Range max query version
template<typename T> struct SegmentTreeMax {
  int n; vector<pair<T,int>> st; const pair<T,int> I = {numeric_limits<T>::min(),-1};
  SegmentTreeMax(const vector<T>& A) : n(A.size()), st(2*n, I) {
    for (int i=0;i<n;i++) st[n+i] = {A[i],i};
    for (int i=n-1; i; --i) st[i] = max(st[2*i], st[2*i+1]);
  }
  void update(int i, int val) {         // Set A[i] = val
    for (st[i+=n] = {val,i}; i > 1; i/= 2) st[i/2] = max(st[i], st[i^1]);
  }
  pair<T,int> query(int l, int r) {     // Find max A[l..r]
    pair<T,int> res = I;
    for (l += n, r += n; l <= r; l /= 2, r /= 2) {
      if (l&1) res = max(res, st[l++]);
      if (~r&1) res = max(res, st[r--]);
    }
    return res;
  }
};

// Verdict: AC
namespace SPOJ_RMQSQ {
  void solve() {
    int N; cin >> N;
    vi A(N); for (auto& x : A) cin >> x;
    SegmentTree<int> st(A);
    int Q; cin >> Q;
    while (Q--) {
      int i, j; cin >> i >> j;
      cout << st.query(i,j).X << '\n';
    }
  }
}

// Verdict: AC
namespace SPOJ_RPLN {
  void solve() {
    int T; cin >> T;
    for (int t=1; t<=T; t++) {
      cout << "Scenario #" << t << ":\n";
      int N, Q; cin >> N >> Q;
      vi A(N); for (auto& x : A) cin >> x;
      SegmentTree<int> st(A);
      while (Q--) {
        int i, j; cin >> i >> j;
        cout << st.query(i-1,j-1).X << '\n';
      }
    }
  }
}

// Verdict: AC
namespace CF_100971FM {
  vvi char_count;
  int diff_chars(int l, int r) {
    vi cnt(26);
    for (int i=0; i<26; i++) cnt[i] = char_count[r+1][i] - char_count[l][i];
    return count_if(cnt.begin(), cnt.end(), [](int x) { return x > 0; });
  }
  void solve() {
    int k, n; string s; cin >> k >> s;
    n = (int)s.size();
    char_count.assign(n+1, vi(26, 0));
    for (int i=1; i<=n; i++) {
      char_count[i] = char_count[i-1];
      char_count[i][s[i-1]-'a']++;
    }
    vi ans(n, -1);  int L = 0, R = 0;
    // Invariant: [L,R] contains all positions j such that [(j+1)..i] is good
    SegmentTree<int> st(vi(n, INT_MAX));
    for (int i=0; i<n; i++) {
      if (diff_chars(0,i) < k) continue;    // ans = -1
      else if (diff_chars(0,i) == k) {      // ans = 1
        ans[i] = 1;
      } else {
        while (diff_chars(L+1,i) > k) L++;
        while (diff_chars(R+1,i) >= k) R++;
        int best_suffix = st.query(L,R-1).X;
        if (best_suffix < INT_MAX) ans[i] = best_suffix + 1;
      }
      if (ans[i] != -1) st.update(i, ans[i]);
    }
    for (int x : ans) cout << x << ' '; cout << endl;
  }
}

// Verdict: AC
namespace SPOJ_FREQUENT {
  void solve() {
  int n, q;
  while (cin >> n && n) {
    cin >> q;
    vi a(n), freq(n), start;
    int f = 1;
    for (int i=0; i<n; i++) {
      cin >> a[i];
      if (i == 0 || a[i] != a[i-1]) f = 1, start.push_back(i);
      freq[i] = f++;
    }
    SegmentTreeMax<int> st(freq);
    while (q--) {
      int i, j; cin >> i >> j; i--, j--;
      if (a[i] == a[j]) cout << j - i + 1 << '\n';
      else {
        auto boundary = lower_bound(start.begin(), start.end(), i);
        if (boundary == start.end()) cout << j - i + 1 << '\n';
        else cout << max(st.query(*boundary, j).X, *boundary - i) << '\n';
      }
    }
  }
}
  
}

int main() {
  ios::sync_with_stdio(0);
  //SPOJ_RMQSQ::solve();
  //SPOJ_RPLN::solve();
  //CF_100971FM::solve();
  SPOJ_FREQUENT::solve();
}
