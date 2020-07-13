// Lowest Common Ancestor test for CodeChef TALCA
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
	ios::sync_with_stdio(0); cin.tie(0);
	int N; cin >> N;
	vvi adj(N);
	for (int i = 0; i < N - 1; i++) {
		int u, v; cin >> u >> v; u--; v--;
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
	int Q; cin >> Q;
	map<int, vvi> queries;
	for (int q = 0; q < Q; q++) {
		int r, u, v; cin >> r >> u >> v;
		r--; u--; v--;
		queries[r].push_back({u, v, q});
	}
	vi ans(Q);
	
	for (auto& qr : queries) {
		int r = qr.first;
		lca T(adj, r);
		for (auto& query : qr.second) {
			int u = query[0], v = query[1], q = query[2];
			ans[q] = T.query(u, v) + 1;
		}
	}
	
	for (int q = 0; q < Q; q++) cout << ans[q] << '\n';
	return 0;
}
