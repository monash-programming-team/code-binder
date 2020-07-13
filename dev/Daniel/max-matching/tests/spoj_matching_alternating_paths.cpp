// Test on SPOJ - MATCHING
//
// Verdict: TLE

#include<bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;

class max_matching {
	int L, R, p;
	vi m, vis;
	vvi adj;
	
	bool dfs(int v) {
		vis[v] = p;
		for (auto u : adj[v]) {
			if (m[u] < 0 || (vis[m[u]] != p && dfs(m[u]))) {
				m[u] = v;
				return 1;
			}
		}
		return 0;
	}
	
 public:
	max_matching(int l, int r) : L(l), R(r), adj(r) { }
	
	void add_edge(int u, int v) { adj[v].push_back(u); }

	pair<int, vi> match() {
		m.assign(L, -1), vis.assign(R, -1);
		int res = 0;
		for (p = 0; p < R; p++) res += dfs(p);
		return {res, m};
	}
};

int main() {
	int N, M, P;
	scanf("%d%d%d", &N, &M, &P);
	max_matching G(N, M);
	for (int i = 0; i < P; i++) {
		int A, B;
		scanf("%d%d", &A, &B);
		A--; B--;
		G.add_edge(A, B);
	}
	printf("%d\n", G.match().first);
}