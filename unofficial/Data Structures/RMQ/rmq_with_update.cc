// Range Minimum/Maximum Query w/ Update
//
// Author      : Xin Wei Chow
// Date        : April 20, 2016
// Reliability : 4
// Tested On   : http://www.spoj.com/problems/RMQSQ/   (tests query only)
//               http://www.spoj.com/problems/RPLN     (tests query only)
//               http://codeforces.com/gym/100971/problem/M
//               ANZAC1 - Freddy and Platforms
//
// Computes the min/max value in an array
// Default is min, 1 change required for max
// Supports point updates of the form v[pos] = val
// 
// Complexity: O(n)     - build segment tree, 
//             O(log n) - query on interval
//             O(log n) - point update
//
// Return:
//
//    pair<T, int> query(int L, int R):
//         Returns a pair of (min/max value, index) from [L..R] inclusive

#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:segment_tree
// Segment tree implementing dynamic range minimum query.
// Complexity: O(n) to build, O(log(n)) to query.
template<typename T> struct SegmentTree {
  typedef pair<T, int> pti;
  int n;  vector<T> v;  vector<pti> st;

  pti merge(pti p1, pti p2) { 
    if (p1.second == -1) return p2;
    if (p2.second == -1) return p1;
    return min(p1, p2); // change this to max for max query
  }
  void build(int p, int L, int R){
    if (L == R) st[p] = {v[L], L};
    else {
      int mid = (L + R) / 2;
      build(p * 2, L, mid);
      build(p * 2 + 1, mid + 1, R);
      st[p] = merge(st[p * 2], st[p * 2 + 1]);
    }
  }
  pti query(int p, int L, int R, int i, int j){
    if (i > R || j < L) return {-1, -1};
    if (i <= L && j >= R) return st[p];
    int mid = (L + R) / 2;
    pti p1 = query(p * 2, L, mid, i, j);
    pti p2 = query(p * 2 + 1, mid + 1, R, i, j);
    return merge(p1, p2);
  }
  void update(int p, int L, int R, int pos, T val){
    if (pos > R || pos < L) return;
    if (pos == L && pos == R) st[p].first = val;
    else {
      int mid = (L + R) / 2;
      update(p * 2, L, mid, pos, val);
      update(p * 2 + 1, mid + 1, R, pos, val);
      st[p] = merge(st[p * 2], st[p * 2 + 1]);
    }
  }

  SegmentTree(vector<T> &v) : n(v.size()), v(v), st(4 * n) { build(1, 0, n-1); }
  // Find the minimum element and its position in the range [i, j]
  pti query(int i, int j){ return query(1, 0, n-1, i, j); }
  // Set the value at position pos to val
  void update(int pos, T val){ update(1, 0, n-1, pos, val); }
};
//listings:/segment_tree

int main(){
  vi v = { 1, 3, 2, 7, 10, -5 };
  SegmentTree<int> cc(v);
  cout << cc.query(0, 5).second << endl;
  cout << cc.query(0, 3).second << endl;
  cout << cc.query(1, 2).second << endl;

  return 0;
}
