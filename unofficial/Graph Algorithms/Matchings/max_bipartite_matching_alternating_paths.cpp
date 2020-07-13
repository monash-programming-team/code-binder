// Maximum bipartite matching via alternating paths.
// Theoretically the slow algorithm but very fast in practice.
//
// Author      : Daniel Anderson
// Date        : 5/10/2016
// Reliability : 5
// Tested on   : UVA10092, UVA11138, UVA10080, ANZAC2016 2E, UVA670
//
// Usage:
//  bipartite_matching G(l, r)    create a graph with l left and r right nodes
//  add_edge(u, v)          add an edge from left node u to right node v
//  match()                 returns the size of the max matching and the matches
//                          the matches indicate each left node's mate, or -1.
//
// Complexity: O(VE)
#include<bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:matching
// Maximum bipartite matching using alternating paths. Complexity: O(VE)
// Returns the number of matches and a vector of each left node's mate.
// match[u] == -1 if node u had no match.
struct bipartite_matching {
	int L, R, p;  vi m, vis;  vvi adj;
	int dfs(int v) {
		vis[v] = p;
		for (auto u : adj[v]) if (m[u] < 0 || (vis[m[u]] != p && dfs(m[u])))
			return m[u] = v, 1;
		return 0;
	}

	bipartite_matching(int l, int r) : L(l), R(r), adj(r) { }
	void add_edge(int u, int v) { adj[v].push_back(u); }
	pair<int, vi> match() {
		int res = 0;  m.assign(L, -1), vis.assign(R, -1);
		for (p = 0; p < R; p++) res += dfs(p);
		return {res, m};
	}
};
//listings:/matching

// Max matching ends -- Example usage (SPOJ - MATCHING)

int main() {
	int N, M, P, A, B;
	scanf("%d%d%d", &N, &M, &P);
	bipartite_matching G(N, M);
	for (int i = 0; i < P; i++) {
		scanf("%d%d", &A, &B);
		G.add_edge(A - 1, B - 1);
	}
	printf("%d\n", G.match().first);
}