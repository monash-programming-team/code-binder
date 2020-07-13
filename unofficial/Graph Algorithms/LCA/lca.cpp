// Lowest Common Ancestor
//
// Author      : Daniel Anderson (Originally based on Xin Wei's LCA)
// Date        : 15/09/2016
// Reliability : 5
// Tested on   : SPOJ - LCA, CodeChef - TALCA, UVA12238, UVA10938, SPOJ - QTREE2
//
#include <bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:lca
// Lowest common ancestor and tree distances using binary lifting.
// Complexity: O(V log(V)) to build, O(log(V)) to query.
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
//listings:/lca

// Verdict: AC
void solve_UVA12238() {
  int N;
  while (cin >> N && N) {
    LCA<ll> lca(N);
    for (int i=1; i<=N-1; i++) {
      int A; ll L; cin >> A >> L;
      lca.add_edge(i,A,L);
    }
    lca.build();
    int Q; cin >> Q;
    for (int q=0; q<Q; q++) {
      if (q) cout << ' ';
      int S, T; cin >> S >> T;
      cout << lca.dist(S,T);
    }
    cout << '\n';
  }
}

// Verdict: AC
void solve_UVA10938() {
  int n;
  while (cin >> n && n) {
    LCA<int> lca(n);
    for (int i=0; i<n-1; i++) {
      int a, b; cin >> a >> b; a--; b--;
      lca.add_edge(a,b);
    }
    lca.build();
    int l; cin >> l;
    while (l--) {
      int a, b; cin >> a >> b; a--; b--;
      int dist = lca.dist(a,b);
      cout << "The fleas ";
      int k = dist / 2;
      int u = lca.lvl[a] > lca.lvl[b] ? a : b;
      while (k > 0) {
        int jump = 0;
        while (lca.anc[u][jump+1] != -1 && k >= (1 << (jump + 1))) jump++;
        u = lca.anc[u][jump];
        k -= (1 << jump);
      }
      if (dist%2==0) cout << "meet at " << u+1;
      else cout << "jump forever between " << min(u, lca.anc[u][0])+1 << " and "
                << max(u, lca.anc[u][0])+1;
      cout << ".\n";
    }
  }
}

// Verdict: AC
void solve_SPOJ_QTREE2() {
  int t; cin >> t;
  while (t--) {
    int N; cin >> N;
    vvi adj(N);
    vector<vector<pair<int,ll>>> wadj(N);
    LCA<ll> wlca(N); LCA<int> uwlca(N);
    for (int i=1; i<=N-1; i++) {
      int a, b, c; cin >> a >> b >> c; a--; b--;
      wlca.add_edge(a,b,c);
      uwlca.add_edge(a,b);
    }
    wlca.build(), uwlca.build();
    string instruction;
    while (cin >> instruction, instruction != "DONE") {
      if (instruction == "DIST") {
        int S, T; cin >> S >> T; S--; T--;
        cout << wlca.dist(S,T) << '\n';
      } else {
        int a, b, k; cin >> a >> b >> k; a--; b--; k--;
        int l = uwlca.query(a, b), dist = uwlca.dist(a, b), u = a;
        if (uwlca.dist(a, l) < k) { u = b; k = dist - k; }
        while (k > 0) {
          int jump = 0;
          while (uwlca.anc[u][jump+1] != -1 && (1 << (jump + 1)) <= k) jump++;
          u = uwlca.anc[u][jump];
          k -= (1 << jump);
        }
        cout << u + 1 << '\n';
      }
    }
    cout << '\n';
  }
}

// Verdict: AC
void solve_SPOJ_LCA() {
  int T; cin >> T;
  for (int t=1; t<=T; t++) {
    cout << "Case " << t << ":\n";
    int N; cin >> N;
    LCA<int> lca(N);
    for (int i=0; i<N; i++) {
      int M; cin >> M;
      for (int j=0; j<M; j++) {
        int u; cin >> u;
        lca.add_edge(i,u-1);
      }
    }
    lca.build();
    int Q; cin >> Q;
    while (Q--) {
      int u, v; cin >> u >> v;
      cout << lca.query(u-1,v-1)+1 << '\n';
    }
  }
}

// Verdict: 
void solve_CC_TALCA() {
  int N; cin >> N;
	LCA<> lca(N);
	for (int i = 0; i < N - 1; i++) {
		int u, v; cin >> u >> v; u--; v--;
		lca.add_edge(u,v);
	}
  lca.build();
	int Q; cin >> Q;
	for (int q = 0; q < Q; q++) {
		int r, u, v; cin >> r >> u >> v;
		r--; u--; v--;
		cout << lca.query(u,v,r)+1 << '\n';
	}
}

int main() {
	//solve_UVA12238();
  //solve_UVA10938();
  //solve_SPOJ_QTREE2();
  //solve_SPOJ_LCA();
  solve_CC_TALCA();
}
