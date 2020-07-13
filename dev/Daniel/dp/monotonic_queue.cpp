// Monotonic queue
//
// Author: Daniel Anderson
// Date: 01-02-2016
// Reliability: 1
// Tested on: CF100971M
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long int ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:monotonic_queue
// A queue that supports amortized O(1) insertion and query
// for the minimum element. Change <= to >= for max element.
template<typename T> struct MonotonicQueue {
  deque<pair<int,T>> q, mins; int cnt = 0;
  void push(T x) {
    while (!mins.empty() && x <= mins.back().Y) mins.pop_back();
    mins.emplace_back(cnt,x), q.emplace_back(cnt++,x);
  }
  void pop() {
    if (mins.front().X == q.front().X) mins.pop_front();
    q.pop_front();
  }
  T front() { return q.front().Y; }
  T min() { return mins.front().Y; }
  bool empty() { return q.empty(); }
};
//listings:/monotonic_queue

namespace problems {
  // Verdict: 
  namespace CF_100971FM {
    vvi char_count;
    int diff_chars(int l, int r) { // Different characters in [l,r] inclusive
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
      vi ans(n, INT_MAX);  int L = 0, R = 0;  MonotonicQueue<int> q;
      // Invariant: [L,R] contains all positions j such that [(j+1)..i] is good
      for (int i=0; i<n; i++) {
        if (diff_chars(0,i) < k) continue;    // ans = -1
        else if (diff_chars(0,i) == k) {      // ans = 1
          ans[i] = 1;
        } else {
          while (diff_chars(L+1,i) > k) { L++; if (!q.empty()) q.pop(); }
          while (diff_chars(R+1,i) >= k) if (++R > L) q.push(ans[R-1]);
          int best_suffix = (q.empty() ? INT_MAX : q.min());
          if (best_suffix < INT_MAX) ans[i] = best_suffix + 1;
        }
      }
      for (int x : ans)
        if (x == INT_MAX) cout << -1 << ' ';
        else cout << x << ' ';
      cout << endl;
    }
  }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  problems::CF_100971FM::solve(); 
}
