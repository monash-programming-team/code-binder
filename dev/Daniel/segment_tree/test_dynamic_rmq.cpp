// Random tests for segment tree RMQ
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

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

template<typename T> struct rmq_update {
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
  rmq_update(vector<T> &v) : n(v.size()), v(v), st(4 * n) { build(1, 0, n-1); }
  // Find the minimum element and its position in the range [i, j]
  pti query(int i, int j){ return query(1, 0, n-1, i, j); }
  // Set the value at position pos to val
  void update(int pos, T val){ update(1, 0, n-1, pos, val); }
};

int main() {
  int T = 1000, N = 50000, Q = 50000;
  for (int t=1; t<=T; t++) {
    srand(time(NULL));
    vi A(N); for (auto& x : A) x = rand();
    SegmentTree<int> st1(A);
    rmq_update<int> st2(A);
    for (int q=1; q<=Q; q++) {
      if (rand()%2) {
        int a = rand() % N, b = rand() % N;
        auto r1 = st1.query(a,b), r2 = st2.query(a,b);
        assert(r1.Y == r2.Y);
        if (r1.Y != -1) assert(r1.X == r2.X);
      } else {
        int i = rand() % N, v = rand();
        st1.update(i,v); st2.update(i,v);
      }
    }
    cout << "Test case " << t << " succeeded" << endl;
  }
}
