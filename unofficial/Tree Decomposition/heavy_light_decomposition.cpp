// Heavy-Light Decomposition
//
// Performs the decomposition into heavy and light edges
// and chains. You must implement your own queries using lambdas.
// See example problems.
//
// Author: Daniel (Based on Xin Wei's HLD)
// Date: 10-01-2017
// Reliability: 5
// Tested on: RCC2016:Q3D, SPOJ - QTREE, TIMMUS1553, SPOJ - GRASSPLA, CodeChef - TAQTREE
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> pii;

//listings:hld
// Heavy-Light Decomposition. Facilitates ranged queries on trees in O(log^2(n))
// time. decompose_tree(root) returns a vector of values that should be
// initialised in a segment tree for queries. ranged_query accepts a path {u..v}
// and a lambda that takes two values i,j that are indices in the segment that
// should be processed in the query. point_query accepts an edge (u, v) and a
// lambda taking the index i corresponding to that edge in the segment tree.
// T is the type of the edge weights.  Complexity: O(n) to build.
template<typename T> struct HeavyLightDecomposition {
  int n; vi heavy, head, par, pos, level; vector<T> cost;
  vector<vector<pair<int,T>>> adj;
  int dfs(int u, int p, int d) {
    int size = 1, max_child = 0, max_child_id = -1;
    par[u] = p, level[u] = d;
    for (auto& child : adj[u]) if (child.X != p) {
      cost[child.X] = child.Y;
      int child_size = dfs(child.X, u, d + 1);
      if (child_size > max_child) max_child = child_size, max_child_id = child.X;
      size += child_size;
    }
    if (max_child * 2 >= size) heavy[u] = max_child_id;
    return size;
  }
  // Create a tree on n vertices -- add edges using add_edge(u, v, cost)
  HeavyLightDecomposition(int n) :
    n(n), heavy(n), head(n), par(n), pos(n), level(n), cost(n), adj(n) { }
  void add_edge(int u, int v, T cost) {
    adj[u].emplace_back(v, cost), adj[v].emplace_back(u, cost);
  }
  vector<T> decompose_tree(int root = 0) {  // Perform HLD.
    vector<T> val(n);  heavy.assign(n, -1);  dfs(root, -1, 0);  int curPos = 0;
    for (int i=0, cur=0; i<n; cur=++i)
      if (par[i] == -1 || heavy[par[i]] != i) while (cur != -1)
        val[curPos] = cost[cur], pos[cur] = curPos++, head[cur] = i, cur = heavy[cur];
    return val;
  }
  template<typename F> void ranged_query(int u, int v, F query) {
    while (head[u] != head[v]) {
      if (level[head[u]] > level[head[v]]) swap(u, v);
      query(pos[head[v]], pos[v]);  v = par[head[v]];
    }
    if (u != v) query(min(pos[u],pos[v])+1, max(pos[u],pos[v]));
  }
  template<typename F> void point_query(int u, int v, F query) {
    query(level[v] > level[u] ? pos[v] : pos[u]);
  }
};
//listings:/hld

// LCA
template<typename T = int> struct LCA {
	const int LOGN = 20;	// works for n <= 10^6. Change appropriately.
	int n;  vi par, lvl;  vvi anc;  vector<T> len;  vector<vector<pair<int,T>>> adj;
	void dfs(int u, int p, int l, T d) {
		par[u] = p, lvl[u] = l, len[u] = d;
		for (auto v : adj[u]) if (v.X != p) dfs(v.X, u, l+1, d+v.Y);
	}
	// Create a tree with n nodes. Add edges then call build(root).
	LCA(int n) : n(n), par(n), lvl(n), len(n), adj(n) { }
  void add_edge(int u, int v, T w = 1) {
    adj[u].emplace_back(v, w), adj[v].emplace_back(u, w);
  }
  void build(int root = 0) {  // Call this before making queries
    dfs(root,-1,0,0), anc.assign(n, vi(LOGN, -1));
		for (int i = 0; i < n; i++) anc[i][0] = par[i];
		for (int k = 1; k < LOGN; k++) for (int i = 0; i < n; i++)
			if (anc[i][k-1] != -1) anc[i][k] = anc[anc[i][k-1]][k-1];
  }
	int query(int u, int v) {  // LCA with respect to original root
		if (lvl[u] > lvl[v]) swap(u, v);
		for (int k = LOGN - 1; k >= 0; k--) 
			if (lvl[v] - (1 << k) >= lvl[u]) v = anc[v][k];
		if (u == v) return u;
		for (int k = LOGN - 1; k >= 0; k--) {
			if (anc[u][k] == anc[v][k]) continue;
			u = anc[u][k]; v = anc[v][k];
		}
		return par[u];
	}
  int query(int u, int v, int root) {  // OPTIONAL: LCA with respect to any root
    int a = query(u, v), b = query(u, root), c = query(v, root);
    if (a == c && c != b) return b;
    else if (a == b && c != b) return c;
    else return a;
  }
	T dist(int u, int v) { return len[u] + len[v] - 2 * len[query(u,v)]; }
};

// Range-Max Segment Tree
template<typename T> struct SegmentTree {
  typedef pair<T, int> pti;
  int n;  vector<T> v;  vector<pti> st;

  pti merge(pti p1, pti p2) { 
    if (p1.second == -1) return p2;
    if (p2.second == -1) return p1;
    return max(p1, p2); // change this to max for max query
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
  pti query(int i, int j){ return query(1, 0, n-1, i, j); }
  void update(int pos, T val){ update(1, 0, n-1, pos, val); }
};

// Interval-sum Segment Tree
class SumSegmentTree {
private:
    int n;
    vi arr, st, lazy;
	
    int merge(int p1, int p2){
		return p1 + p2;
    }
	
    void build(int p, int L, int R){
		if (L == R) st[p] = arr[L];
		else {
			int mid = (L + R) / 2;
			build(p * 2, L, mid);
			build(p * 2 + 1, mid + 1, R);
			st[p] = merge(st[p * 2], st[p * 2 + 1]);
		}
    }
	
    int query(int p, int L, int R, int i, int j){
		if (i > R || j < L) return 0;
		if (i <= L && j >= R) return st[p];
		int mid = (L + R) / 2;
		int p1 = query(p * 2, L, mid, i, j);
		int p2 = query(p * 2 + 1, mid + 1, R, i, j);
		return merge(p1, p2);
    }
	
    void update(int p, int L, int R, int pos, int val){
		if (pos > R || pos < L) return;
		if (pos == L && pos == R) st[p] = val;
		else {
			int mid = (L + R) / 2;
			update(p * 2, L, mid, pos, val);
			update(p * 2 + 1, mid + 1, R, pos, val);
			st[p] = merge(st[p * 2], st[p * 2 + 1]);
		}
    }

public:
    SumSegmentTree(int n, const vi &arr) : n(n), arr(arr) {
		st.assign(4 * n, 0);
		lazy.assign(4 * n, 0);
		build(1, 0, n-1);
    }
    int query(int i, int j){
		return query(1, 0, n-1, i, j);
    }
    void update(int pos, int val){
		update(1, 0, n-1, pos, val);
    }
};

// Ranged update, point query Fenwick tree
struct FenwickTree {
	int N;  vi A;
	FenwickTree(int n): N(n+1), A(N) {}
	void adjust(int b,int v) { for (;b;b-=b&-b) A[b]+=v; }
	void adjust(int a,int b,int v) { adjust(b,v); adjust(a,-v); }
	int pq(int i) { int r=0; for (i++;i<N;i+=i&-i) r+=A[i]; return r;	}
};

// Ranged query, point update Fenwick tree
struct FenwickTree2 {
	int N;  vi A;
	FenwickTree2(int n): N(n+1), A(N) {}
	int rq(int b) { int sum=0; for (;b;b-=b&-b) sum+=A[b]; return sum; }
	int rq(int a,int b) { return rq(b)-rq(a); }
	void adjust(int i,int v) { for (i++;i<N;i+=i&-i) A[i]+=v; }
};

// ---------------------------------------------------------
//              SPOJ - QTREE. Verdict: AC
// Tested against Xin Wei's code on randomly generated test
// cases since SPOJ won't accept C++14 submissions for this
// problem.
// ---------------------------------------------------------
void solve_SPOJ_QTREE() {
  int t; cin >> t;
  while (t--) {
    int N; cin >> N;
    HeavyLightDecomposition<int> hld(N);
    vector<pii> edges;
    for (int i=0; i<N-1; i++) {
      int a, b, c; cin >> a >> b >> c;
      hld.add_edge(a-1,b-1,c);
      edges.emplace_back(a-1,b-1);
    }
    auto arr = hld.decompose_tree();
    SegmentTree<int> st(arr);
    string type;
    while (cin >> type, type != "DONE") {
      if (type == "QUERY") {
        int u, v; cin >> u >> v; u--; v--;
        int best = INT_MIN;
        hld.ranged_query(u, v, [&](int i, int j) { best = max(best, st.query(i,j).X); });
        cout << best << '\n';
      } else {
        int i, ti; cin >> i >> ti; i--;
        hld.point_query(edges[i].X, edges[i].Y, [&](int x) { st.update(x, ti); });
      }
    }
  }
}

// ---------------------------------------------------------
//            RCC2016 - QUAL3: D. Verdict: AC
// ---------------------------------------------------------
void solve_RCC2016Q3D() {
  int n; cin >> n;
  HeavyLightDecomposition<int> hld(n);
  vi par(n, -1);
  for (int i=1; i<n; i++) {
    int p; cin >> p; p--;
    hld.add_edge(i, p, 1);
    par[i] = p;
  }
  auto arr = hld.decompose_tree(0);
  SumSegmentTree st(n, arr);
  int q; cin >> q;
  while (q--) {
    int type; cin >> type;
    if (type == 1) {
      int a, b; cin >> a >> b; a--, b--;
      int ans = 0;
      hld.ranged_query(a, b, [&](int i, int j) { ans += st.query(i,j); });
      cout << ans << '\n';
    } else {
      int v; cin >> v; v--;
      hld.point_query(v, par[v], [&](int i) { st.update(i, 0); });
    }
  }
}

// ---------------------------------------------------------
//              TIMMUS 1552. Verdict: AC
// ---------------------------------------------------------
void solve_timmus1553() {
  int N; cin >> N;
  HeavyLightDecomposition<ll> hld(N);
  LCA<int> lca(N);
  for (int i=0; i<N-1; i++) {
    int a, b; cin >> a >> b; a--, b--;
    hld.add_edge(a, b, 0);
    lca.add_edge(a, b);
  }
  lca.build();  vector<ll> radiation(N);
  auto arr = hld.decompose_tree();  SegmentTree<ll> st(arr);
  int Q; cin >> Q;
  while (Q--) {
    char C; int U, V; cin >> C >> U >> V;
    if (C == 'I') {
      U--; radiation[U] += V;
      if (hld.par[U] != -1) hld.point_query(U, hld.par[U], [&](int i) { st.update(i, radiation[U]); });
    } else {
      ll ans = LLONG_MIN;
      U--, V--;
      hld.ranged_query(U, V, [&](int i, int j) { ans = max(ans, st.query(i,j).X); });
      ans = max(ans, radiation[lca.query(U, V)]);
      cout << ans << '\n';
    }
  }
}

// ---------------------------------------------------------
//              LIGHTOJ - 1348. Verdict: n/a
//   Can't judge because LightOJ does not support C++11.
// ---------------------------------------------------------
void solve_light1348() {
  int T; cin >> T;
  for (int t=1; t<=T; t++) {
    cout << "Case " << t << ":\n";
    int n; cin >> n;
    vi population(n);
    for (auto& x : population) cin >> x;
    HeavyLightDecomposition<int> hld(n);
    LCA<int> lca(n);
    for (int i=0; i<n-1; i++) {
      int u, v; cin >> u >> v;
      hld.add_edge(u, v, 0);
      lca.add_edge(u, v);
    }
    lca.build();
    SumSegmentTree st(n, hld.decompose_tree());
    for (int i=0; i<n; i++) if (hld.par[i] != -1)
      hld.point_query(i, hld.par[i], [&](int j) { st.update(j, population[i]); });
    int q; cin >> q;
    while (q--) {
      int type; cin >> type;
      if (type == 0) {
        int i, j; cin >> i >> j;
        int ans = 0;
        hld.ranged_query(i, j, [&](int k, int m) { ans += st.query(k, m); });
        ans += population[lca.query(i,j)];
        cout << ans << '\n';
      } else {
        int i, v; cin >> i >> v;
        population[i] = v;
        if (hld.par[i] != -1) hld.point_query(i, hld.par[i], [&](int j) { st.update(j, v); });
      }
    }
  }
}

// ---------------------------------------------------------
//              SPOJ - GRASSPLA. Verdict: AC
// ---------------------------------------------------------
void solve_SPOJ_GRASSPLA() {
  int N, M; cin >> N >> M;
  HeavyLightDecomposition<ll> hld(N);
  for (int i=0; i<N-1; i++) {
    int a, b; cin >> a >> b; a--, b--;
    hld.add_edge(a, b, 0);
  }
  hld.decompose_tree();
  FenwickTree ft(N);
  while (M--) {
    char c; cin >> c; int a, b; cin >> a >> b; a--, b--;
    if (c == 'P') hld.ranged_query(a, b, [&](int i,int j) { ft.adjust(i,j+1,1); });
    else {
      int ans = 0;
      hld.point_query(a, b, [&](int i) { ans += ft.pq(i); });
      cout << ans << '\n';
    }
  }
}


// ---------------------------------------------------------
//              CodeChef - TAQTREE. Verdict: AC
// ---------------------------------------------------------
void solve_CC_TAQTREE() {
  int N; cin >> N;
  HeavyLightDecomposition<int> hld(N);
  vector<pii> edges;
  for (int i=0; i<N-1; i++) {
    int a, b, c; cin >> a >> b >> c; a--, b--;
    hld.add_edge(a,b,c);
    edges.emplace_back(a,b);
  }
  auto arr = hld.decompose_tree();
  FenwickTree2 ft(arr.size());
  for (int i=0; i<(int)arr.size(); i++) ft.adjust(i, arr[i]);
  int Q; cin >> Q;
  while (Q--) {
    int type; cin >> type;
    if (type == 1) {
      int i, c; cin >> i >> c; i--;
      hld.point_query(edges[i].X, edges[i].Y, [&](int j) {
        int cur = ft.rq(j,j+1);
        ft.adjust(j,c-cur);
      });
    } else {
      int u, v; cin >> u >> v; u--, v--;
      int ans = 0;
      hld.ranged_query(u, v, [&](int i, int j) { ans += ft.rq(i, j+1); });
      cout << ans << '\n';
    }
  }
}

int main() {
  ios_base::sync_with_stdio(0); cin.tie(0);
  //solve_SPOJ_QTREE();
  //solve_RCC2016Q3D();
  //solve_timmus1553();
  //solve_light1348();
  //solve_SPOJ_GRASSPLA();
  solve_CC_TAQTREE();
}
