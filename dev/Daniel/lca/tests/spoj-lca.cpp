// Lowest Common Ancestor test on SPOJ - LCA
//
#include <bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:lca
// Lowest common ancestor using binary lifting.
// Complexity: O(V log(V)) to build, O(log(V)) to query.
struct lca {
	const int LOGN = 20;	// works for n <= 10^6
	vvi& adj;  int n;  vi par, lvl;  vvi anc;
	void dfs(int u, int p, int d) {
		par[u] = p; lvl[u] = d;
		for (int v : adj[u]) if (v != p) dfs(v, u, d+1);
	}
	
	lca(vvi& adj, int root = 0) : adj(adj), n(adj.size()), par(n), lvl(n) {
		dfs(root,-1,0);
		anc.assign(n, vi(LOGN, -1));
		for (int i = 0; i < n; i++) anc[i][0] = par[i];
		for (int k = 1; k < LOGN; k++) for (int i = 0; i < n; i++)
			if (anc[i][k-1] != -1) anc[i][k] = anc[anc[i][k-1]][k-1];
	}
	int query(int u, int v) {
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
	int dist(int u, int v) { return lvl[u] + lvl[v] - 2 * lvl[query(u,v)]; }
};


int main() {
	int T; cin >> T;
	int c = 1;
	while (T--) {
		cout << "Case " << c++ << ":" << '\n';
		int N; cin >> N;
		vvi adj(N);
		for (int u = 0; u < N; u++) {
			int M; cin >> M;
			for (int j = 0; j < M; j++) {
				int v; cin >> v; v--;
				adj[u].push_back(v);
			}
		}
		lca T(adj);
		int Q; cin >> Q;
		for (int q = 0; q < Q; q++) {
			int v, w; cin >> v >> w;
			v--; w--;
			int res = T.query(v, w);
			cout << res + 1 << '\n';
		}
	}
	return 0;
}
