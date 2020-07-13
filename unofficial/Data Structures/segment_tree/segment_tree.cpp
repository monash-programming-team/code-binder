// Modular Segment Tree with lazy propagation.
//
// Look at the plethora of examples to learn how to use it.
//
// Author: Daniel Anderson
// Date: 17-01-2017
// Reliability: 5
// Tested on: SPOJ-RMQSQ (Range min query only)
//            SPOJ-RPLN (Range min query only)
//            SPOJ-FREQUENT(Range max query only)
//            UVA12532 (Range product and point assignment)
//            SPOJ-SEGSQRSS (Range sum-of-squares, assignment and addition)
//            UVA11402 (Range sum, assignment and inverse of binary values)
//            SPOJ-GSS3 (Max contiguous sub-sum in a range and point assignment)
//            IEEE10-Safety (String hashes with ranged replacement and increment)
//
// Complexity: O(n) to build. O(log(n)) to query.
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:lazy_segment_tree
// Modular segment tree with lazy propagation. You must subclass this to implement
// your own queries and updates. Specifically:
//  -- implement a constructor that sets n(size), I(identity) and calls build(array)
//  -- override merge to define how to merge children into their parent node
//  -- override apply to define how to apply an update to a segment
//  -- override compose to define how to compose successive updates
// See RangeSumAddition example. Complexity: O(n) to build, O(log(n)) to query.
// Parameters: T = Segment data type, U = Update data type
template<typename T, typename U> struct SegmentTree {
  int n; const T I;  vi L, R;
  vector<T> st;  vector<U> lazy;  vector<bool> pending;
  virtual T merge(const T& a, const T& b) = 0;     // a and b are ordered left & right
  virtual void compose(U& a, const U& b) = 0;       // update b happens after update a
  virtual void apply(int i, int j, T& seg, const U& upd) = 0; // seg ranges over [i,j]
  void build(const vector<T>& A) { build(0, 0, n-1, A); }
  void build(int p, int i, int j, const vector<T>& A) {
    L[p] = i, R[p] = j;
    if (i == j) { st[p] = A[i]; return; }
    build(2*p+1, i, (i+j)/2, A);
    build(2*p+2, (i+j)/2+1, j, A);
    st[p] = merge(st[2*p+1], st[2*p+2]);
  }
  void push(int p) {
    if (!pending[p]) return;  pending[p] = 0;
    apply(L[p], R[p], st[p], lazy[p]);
    if (L[p] == R[p]) return;
    if (pending[2*p+1]) compose(lazy[2*p+1], lazy[p]);
    else lazy[2*p+1] = lazy[p], pending[2*p+1] = 1;
    if (pending[2*p+2]) compose(lazy[2*p+2], lazy[p]);
    else lazy[2*p+2] = lazy[p], pending[2*p+2] = 1;
  }
  void update(int p, int i, int j, const U& upd) {          
    push(p);
    if (i <= L[p] && R[p] <= j) { lazy[p] = upd, pending[p] = 1; return; }
    else if (i > R[p] || j < L[p]) return;
    update(2*p+1,i,j,upd), update(2*p+2,i,j,upd);
    push(2*p+1), push(2*p+2);
    st[p] = merge(st[2*p+1], st[2*p+2]);
  }
  T query(int p, int i, int j) {                             
    push(p);
    if (i <= L[p] && R[p] <= j) { return st[p]; }
    else if (i > R[p] || j < L[p]) return I;
    return merge(query(2*p+1,i,j),query(2*p+2,i,j));
  }
  SegmentTree(int n, const T& I) : n(n), I(I), L(4*n,-1), R(4*n,-1),
    st(4*n,I), lazy(4*n), pending(4*n) { }
  void update(int i, int j, const U& upd) { update(0, i, j, upd); }   // update [i,j]
  T query(int i, int j) { return query(0, min(i,j), max(i,j)); }    // query on [i,j]
};
//listings:/lazy_segment_tree

// -------------------------------------------------------------------------------------
// Implementation of Segment Tree supporting range minimum query and ranged assignment.
// -------------------------------------------------------------------------------------
struct RangeMinAssign : public SegmentTree<int,int> {
  RangeMinAssign(const vi& A) : SegmentTree(A.size(),INT_MAX) { build(A); }
  int merge(const int& a, const int& b) { return min(a, b); }
  void compose(int& prev, const int& newv) { prev = newv; }
  void apply(int i, int j, int& seg, const int& upd) { seg = upd; }
};

// -------------------------------------------------------------------------------------
//listings:range_sum_addition
// Example: Implementation of Segment Tree supporting range sum and range addition
struct RangeSumAddition : public SegmentTree<int,int> {
  RangeSumAddition(const vi& A) : SegmentTree(A.size(),0) { build(A); }
  int merge(const int& a, const int& b) { return a + b; }
  void compose(int& prev, const int& newv) { prev += newv; }
  void apply(int i, int j, int& seg, const int& upd) { seg += (j - i + 1) * upd; }
};
//listings:/range_sum_addition

// -------------------------------------------------------------------------------------
// Implementation of Segment Tree supporting range maximum query and ranged assignment.
// -------------------------------------------------------------------------------------
struct RangeMaxAssign : public SegmentTree<int,int> {
  RangeMaxAssign(const vi& A) : SegmentTree(A.size(),INT_MIN) { build(A); }
  int merge(const int& a, const int& b) { return max(a, b); }
  void compose(int& prev, const int& newv) { prev = newv; }
  void apply(int i, int j, int& seg, const int& upd) { seg = upd; }
};

// -------------------------------------------------------------------------------------
// Implementation of Segment Tree supporting range product and ranged assignment.
// -------------------------------------------------------------------------------------
struct RangeProductAssign : public SegmentTree<int,int> {
  RangeProductAssign(const vi& A) : SegmentTree(A.size(),1) { build(A); }
  int merge(const int& a, const int& b) { return a * b; }
  void compose(int& prev, const int& newv) { prev = newv; }
  void apply(int i, int j, int& seg, const int& upd) { seg = upd; }
};

// -------------------------------------------------------------------------------------
// Implementation of Segment Tree supporting range sum-of-squares, ranged assignment
// and ranged addition.
// -------------------------------------------------------------------------------------
typedef pair<ll,ll> pll;
struct RangeSumOfSquares : public SegmentTree<pll, pair<ll,bool>> {
  RangeSumOfSquares(const vi& A) : SegmentTree(A.size(),{0,0}) {
    vector<pair<ll,ll>> B;
    for (auto& x : A) B.emplace_back(x, x*x);
    build(B);
  }
  pll merge(const pll& a, const pll& b) { return {a.X + b.X, a.Y + b.Y}; }
  void compose(pair<ll,bool>& prev, const pair<ll,bool>& newv) {
    if (newv.Y) prev.X += newv.X;   // Ranged addition
    else prev.X = newv.X;           // Ranged assignment
  }
  void apply(int i, int j, pll& seg, const pair<ll,bool>& upd) {
    if (upd.Y) {  // Ranged addition
      seg.Y += 2 * seg.X * upd.X + (j - i + 1) * upd.X * upd.X;
      seg.X += (j - i + 1) * upd.X;
    } else {    // Ranged assignment
      seg.X = (j - i + 1) * upd.X;
      seg.Y = (j - i + 1) * upd.X * upd.X;
    }
  }
  void set_range(int i, int j, ll v) { update(i,j,{v,0}); } // {v,false} means set
  void add_range(int i, int j, ll v) { update(i,j,{v,1}); } // {v,true} means add
  ll sum_query(int i, int j) { return query(i,j).Y; }
};

// -------------------------------------------------------------------------------------
// Implementation of Segment Tree supporting range sum, ranged assignment and
// range inverse (binary NOT) of a binary sequence.
// -------------------------------------------------------------------------------------
enum PirateUpdateType { SET, CLEAR, FLIP, NOTHING };
struct RangeBinaryInverse : public SegmentTree<int,PirateUpdateType> {
  RangeBinaryInverse(const vi& A) : SegmentTree(A.size(),0) { build(A); }
  int merge(const int& a, const int& b) { return a + b; }
  void compose(PirateUpdateType& prev, const PirateUpdateType& newv) {
    if (newv == NOTHING) return;      // Do nothing
    if (prev == NOTHING || newv == SET || newv == CLEAR) prev = newv;
    else if (prev == CLEAR) prev = SET;   // Flip a CLEAR -> becomes a SET
    else if (prev == SET) prev = CLEAR;   // Flip a SET -> becomes a CLEAR
    else prev = NOTHING;                  // Flip a FLIP -> becomes NOTHING
  }
  void apply(int i, int j, int& seg, const PirateUpdateType& upd) {
    if (upd == NOTHING) return;
    else if (upd == SET) seg = j - i + 1;   // Set range to 1
    else if (upd == CLEAR) seg = 0;         // Set range to 0
    else seg = (j - i + 1) - seg;           // Flip range
  }
  void set_range(int i, int j) { update(i,j,SET); }
  void clear_range(int i, int j) { update(i,j,CLEAR); }
  void flip_range(int i, int j) { update(i,j,FLIP); }
  int sum_query(int i, int j) { return query(i,j); }
};

// -------------------------------------------------------------------------------------
// Implementation of Segment Tree supporting finding the maximum (non-empty)
// contiguous sub sum in a range and ranged assignment.
// -------------------------------------------------------------------------------------
struct Sums { int total_sum=0, max_sum=0, max_prefix=0, max_suffix=0; };
constexpr Sums IdSums = {0, -10000, -10000, -10000};
struct RangeMaxSubsumAssign : public SegmentTree<Sums,int> {
  RangeMaxSubsumAssign(const vi& A) : SegmentTree(A.size(),IdSums) {
    vector<Sums> B;
    for (auto& x : A) B.push_back({x,x,x,x});
    build(B);
  }
  Sums merge(const Sums& a, const Sums& b) {
    Sums result;
    result.total_sum = a.total_sum + b.total_sum;
    result.max_sum = max({a.max_sum, b.max_sum, a.max_suffix + b.max_prefix});
    result.max_prefix = max(a.max_prefix, a.total_sum + b.max_prefix);
    result.max_suffix = max(b.max_suffix, b.total_sum + a.max_suffix);
    return result;
  }
  void compose(int& prev, const int& newv) { prev = newv; }
  void apply(int i, int j, Sums& seg, const int& upd) {
    int total = (j - i + 1) * upd;
    if (upd >= 0) seg = {total, total, total, total};
    else seg = {total, upd, upd, upd};
  }
  int max_sum_query(int i, int j) { return query(i,j).max_sum; }
};


// -------------------------------------------------------------------------------------
// Implementation of Segment Tree for string hashing. Used in IEEEExtreme 10.0
// problem "Safety" - Supports ranged replacement and ranged cyclic increment
// -------------------------------------------------------------------------------------
// Segment Node - Records 26 hashes, one for each possible increment
struct Segment {  
  int len, cur; vi hashes;
  Segment(int len, int cur) : len(len), cur(cur), hashes(26) { }
};
const Segment IdentitySegment(0,0);
// Update node, records updates and accumulated increments
struct Update { int type, i, j, k, add; }; 
// "Safety" Tree implementation
struct SafetyTree : SegmentTree<Segment, Update> {
  constexpr static int MOD = 1000000637; constexpr static int BASE = 27;
  int n; vvi hashes; vi bpow;
  // Fast modding arithmetic for hashing --------------------------------------------
  inline int add(int a, int b) { a += b; if (a >= MOD) a-= MOD; if (a < 0) a += MOD; return a; }
  inline int add_alph(int a, int b) { a += b; if (a >= 26) a -= 26; return a; }
  // Hash concatenation
  inline int concat_hashes(int a, int b, int len) { return (a + (ll)bpow[len] * b) % MOD; }
  // Merge: Concatenate two adjacent hashes -----------------------------------------
  Segment merge(const Segment& a, const Segment& b) {
    Segment res(a.len+b.len,0);
    for (int i=0; i<26; i++) res.hashes[i] =
      concat_hashes(a.hashes[add_alph(i,a.cur)], b.hashes[add_alph(i,b.cur)], a.len);
    return res;
  }
  // Update: Compose two updates ----------------------------------------------------
  void compose(Update& old, const Update& newv) {
    if (newv.type == 2) {  // Replacement
      old = newv;
    } else {  // Increment
      old.add = add_alph(old.add, newv.add);
    }
  }
  // Update: Apply an update to a segment --------------------------------------------
  void apply(int i, int j, Segment& seg, const Update& upd) {
    if (upd.type == 2) {    // Replace A[i..j] with original A[k..] + merged increments
      int k1 = upd.k + (i - upd.i), k2 = k1 + (j - i);
      for (int o=0; o<26; o++) {
        int lvl = add_alph(o, upd.add);
        int hash = (((ll)hashes[lvl][k2] -
          (k1 ? hashes[lvl][k1-1] : 0) + MOD) * bpow[n - k1 - 1]) % MOD;
          seg.hashes[o] = hash;
      }
      seg.cur = 0;
    } else {  // Increment the range A[i..j]
      seg.cur = add_alph(seg.cur, upd.add);
    }
  }
  // Build the tree from the given string. Pre-computes hashes -----------------------
  SafetyTree(const string& S) : SegmentTree(S.size(),IdentitySegment), n(S.size()),
                                hashes(26, vi(n)), bpow(n,1) {
    for (int i=1; i<n; i++) bpow[i] = ((ll)bpow[i-1] * BASE) % MOD;
    for (int o=0; o<26; o++) for (int i=0; i<n; i++) {
      hashes[o][i] = ((ll)bpow[i] * add_alph(S[i]-'a',o)) % MOD;
      if (i > 0) hashes[o][i] = add(hashes[o][i], hashes[o][i-1]);
    }
    vector<Segment> segments(n, Segment(1,0));
    for (int i=0; i<n; i++) for (int o=0; o<26; o++)
      segments[i].hashes[o] = ((ll)bpow[n-1] * add_alph(S[i]-'a',o)) % MOD;
    build(segments);
  }
  // Query for just the current hash, rather than all 26. This speed-up gets AC ------
  int fast_query(int p, int i, int j) {
    push(p);
    if (i <= L[p] && R[p] <= j) { return st[p].hashes[st[p].cur]; }
    else if (i > R[p] || j < L[p]) return 0;
    int len = (L[2*p+1]==-1) ? 0 : max(0, min(j, R[2*p+1]) - max(i, L[2*p+1]) + 1);
    return concat_hashes(fast_query(2*p+1,i,j),fast_query(2*p+2,i,j),len);
  }
  // Problem queries -----------------------------------------------------------------
  bool equal(int i, int j, int k) {           // Are the ranges A[i..j] == A[k..] ?
    return fast_query(0,i,j) == fast_query(0,k,k+j-i);
  }
  void increment_range(int i, int j) {        // Increment the range A[i..j]
    update(i, j, {3, i, j, -1, 1});
  }
  void replace_range(int i, int j, int k) {   // Replace the range A[i..j] with A*[k..]
    update(i, j, {2, i, j, k, 0});
  }
};

// Make Segment Trees printable via output streams
template<typename T, typename U>
ostream& operator<<(ostream& o, SegmentTree<T,U>& st) {
  for (int i=0; i<st.n; i++) o << st.query(i,i) << ' ';
  return o;
}

// ----------------------------------------------------------------------------
//                                 TEST PROBLEMS
// ----------------------------------------------------------------------------

// Verdict: AC
void solve_SPOJ_RMQSQ() {
  int N; cin >> N;
  vi A(N); for (auto& x : A) cin >> x;
  RangeMinAssign st(A);
  int Q; cin >> Q;
  while (Q--) {
    int i, j; cin >> i >> j;
    cout << st.query(i,j) << '\n';
  }
}

// Verdict: AC
void solve_SPOJ_RPLN() {
  int T; cin >> T;
  for (int t=1; t<=T; t++) {
    cout << "Scenario #" << t << ":\n";
    int N, Q; cin >> N >> Q;
    vi A(N); for (auto& x : A) cin >> x;
    RangeMinAssign st(A);
    while (Q--) {
      int i, j; cin >> i >> j;
      cout << st.query(i-1,j-1) << '\n';
    }
  }
}

// Verdict: AC
void solve_UVA12532() {
  int N, K;
  while (cin >> N >> K) {
    vi A(N); for (auto& x : A) { cin >> x; if (x) x /= abs(x); }
    RangeProductAssign st(A);
    while (K--) {
      char type; cin >> type;
      if (type == 'C') {
        int I, V; cin >> I >> V;
        if (V) V /= abs(V);
        st.update(I-1,I-1,V);
      } else {
        int I, J; cin >> I >> J;
        int res = st.query(I-1,J-1);
        if (res == 1) cout << '+';
        else if (res == -1) cout << '-';
        else cout << '0';
      }
    }
    cout << '\n';
  }
}

// Verdict: AC
void solve_SPOJ_FREQUENT() {
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
    RangeMaxAssign st(freq);
    while (q--) {
      int i, j; cin >> i >> j; i--, j--;
      if (a[i] == a[j]) cout << j - i + 1 << '\n';
      else {
        auto boundary = lower_bound(start.begin(), start.end(), i);
        if (boundary == start.end()) cout << j - i + 1 << '\n';
        else cout << max(st.query(*boundary, j), *boundary - i) << '\n';
      }
    }
  }
}

// Verdict: AC
void solve_SPOJ_SEGSQRSS() {
  int T; cin >> T;
  for (int t=1; t<=T; t++) {
    cout << "Case " << t << ":\n";
    int N, Q; cin >> N >> Q;
    vi a(N); for (auto& x : a) cin >> x;
    RangeSumOfSquares st(a);
    while (Q--) {
      int type; cin >> type;
      if (type == 0) {  // Ranged assignment
        int i, j, x; cin >> i >> j >> x; i--, j--;
        st.set_range(i,j,x);
      } else if (type == 1) { // Ranged addition
        int i, j, x; cin >> i >> j >> x; i--, j--;
        st.add_range(i,j,x);
      } else {  // Query for sum of squares
        int i, j; cin >> i >> j; i--, j--;
        cout << st.sum_query(i,j) << '\n';
      }
    }
  }
}

// Verdict: AC
void solve_UVA11402() {
  int T; cin >> T;
  for (int t=1; t<=T; t++) {
    cout << "Case " << t << ":\n";
    vi pirates;
    int M; cin >> M;
    while (M--) {
      int T; string P; cin >> T >> P;
      while (T--) for (auto c : P) pirates.push_back(c == '1');
    }
    RangeBinaryInverse st(pirates);
    int Q, q = 1; cin >> Q;
    while (Q--) {
      char type; int a, b; cin >> type >> a >> b;
      if (type == 'F')        st.set_range(a,b);    // Set range to 1
      else if (type == 'E')   st.clear_range(a,b);  // Clear range to 0
      else if (type == 'I')   st.flip_range(a,b);   // Flip range
      else cout << "Q" << q++ << ": " << st.sum_query(a,b) << '\n';
    }
  }
}

// Verdict: AC
void solve_SPOJ_GSS3() {
  int N; cin >> N;
  vi a(N);
  for (auto& x : a) cin >> x;
  RangeMaxSubsumAssign st(a);
  int M; cin >> M;
  while (M--) {
    int type, x, y; cin >> type >> x >> y;
    if (type == 0) st.update(x-1,x-1,y);              // Set A[x] = y
    else cout << st.max_sum_query(x-1,y-1) << '\n';   // Print max
  }
}

// Verdict: AC
void solve_IEEE10_safety() {
  string S; cin >> S;
  int M; cin >> M;
  SafetyTree tree(S);
  while (M--) {
    int type; cin >> type;
    if (type == 1) {            // Is A[i..j] == A[k..] ?
      int i, j, k; cin >> i >> j >> k;
      if (tree.equal(i-1,j-1,k-1)) cout << "Y\n";
      else cout << "N\n";
    } else if (type == 2) {     // Set A[i..j] to original A[k..]
      int i, j, k; cin >> i >> j >> k;
      tree.replace_range(i-1,j-1,k-1);
    } else {                    // Increment A[i..j]
      int i, j; cin >> i >> j;
      tree.increment_range(i-1,j-1);
    }
  }
}

int main() {
  ios_base::sync_with_stdio(0); cin.tie(0);
  //solve_SPOJ_RMQSQ();
  //solve_SPOJ_RPLN();
  //solve_UVA12532();
  //solve_SPOJ_FREQUENT();
  //solve_SPOJ_SEGSQRSS();
  //solve_UVA11402();
  //solve_SPOJ_GSS3();
  solve_IEEE10_safety();
}
