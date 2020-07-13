// Minimum cost arborescence (directed minimum spanning tree) in a directed
// graph.
//
// Author: Daniel Anderson
// Date: 21-01-2017
// Reliability: 1
// Tested on: UVA11182
//
// Complexity: O(VE)
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

// Strongly-connected components. Required to compute minimum cost arborescence.
struct SCC {
	int n, comp;  vvi g, gt;  vi seq, vis;
	void dfs(int u, const vvi &adj) {
		for (int v : adj[u]) if (vis[v] == -1) { vis[v] = comp; dfs(v, adj); }
		seq.push_back(u);
	}
	SCC(int N) : n(N), g(n), gt(n) { }
	void add_edge(int u, int v) { g[u].push_back(v); gt[v].push_back(u); }
	pair<int, vi> find_SCC() {
		vis.assign(n, -1); comp = 0;
		for (int i = 0; i < n; i++) if (vis[i] == -1) { vis[i] = comp; dfs(i, g); }
		vis.assign(n, -1); comp = 0;
		for (int i = n-1; i >= 0; i--) {
			int u = seq[i];
			if (vis[u] == -1) { vis[u] = comp; dfs(u, gt); comp++; }
		}
		return {comp, vis};
	}
};

//listings:mca
// Computes the minimum cost arborescence (directed minimum spanning tree) from root
// in a directed graph. Returns INF if no arborescence exists. Complexity: O(VE)
template<typename T> struct MinCostArborescence {
	typedef vector<vector<pair<int,T>>> Graph;  Graph adj;
	int n; const T INF = numeric_limits<T>::max() / 2;
	MinCostArborescence(int N) : adj(N), n(N) { }
	void add_edge(int u, int v, T w) { adj[u].emplace_back(v, w); }
	T find(int root) { return find(root, adj); }
	T find(int root, const Graph& G) {
		int nv = (int)G.size();  T res = 0;  vector<T> mins(nv, INF);
		for (int v=0; v<nv; v++) for (auto& e : G[v]) mins[e.X]=min(mins[e.X],e.Y);
		for (int v=0; v<nv; v++) if (v != root) {
			if (mins[v] == INF) return INF; else res += mins[v];
		}
		SCC scc(nv);  // Include Strongly-connected components code
		for (int v=0; v<nv; v++) for (auto& e : G[v]) if (e.X != root)
			if (e.Y - mins[e.X] == 0) scc.add_edge(v, e.X);
		int m; vi comp;  tie(m, comp) = scc.find_SCC();  Graph G2(m);
		if (m == nv) return res;
		for (int v=0; v<nv; v++) for (auto& e : G[v]) if (comp[v] != comp[e.X])
			G2[comp[v]].emplace_back(comp[e.X], e.Y - mins[e.X]);
		return min(INF, res + find(comp[root], G2));
	}
};
//listings:/mca

// ----------------------------------------------------------------------------
//                                TEST PROBLEMS
// ----------------------------------------------------------------------------
namespace problems {
  // Verdict: AC
  namespace UVA_11182 {
    void solve() {
      int N; cin >> N;
      for (int t=1; t<=N; t++) {
        cout << "Case #" << t << ": ";
        int n, m; cin >> n >> m;
        MinCostArborescence<int> G(n);
        while (m--) {
          int u, v, w; cin >> u >> v >> w;
          G.add_edge(u,v,w);
        }
        int res = G.find(0);
        if (res == G.INF) cout << "Possums!\n";
        else cout << res << '\n';
      }
    }
  }
}

int main() {
	problems::UVA_11182::solve();
}
