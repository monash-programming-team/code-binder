// Maximum bipartite matching via alternating paths.
// Theoretically the slow algorithm but very fast in practice.
//
// Author      : Daniel Anderson
// Date        : 5/10/2016
// Reliability : 2
// Tested on   : UVA10092, UVA11138
//
// Complexity: O(nm)
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
	
	int T, n, b;
	cin >> T;
	
	for (int t = 1; t <= T; t++) {
    	cin >> b >> n;
    	
    	max_matching G(b, n);
    	
    	int i;
    	for (int c = 0; c < b; c++) {
    		for (int d = 0; d < n; d++) {
    			cin >> i;
    			if (i) G.add_edge(c, d);
    		}
    	}
    	
    	cout << "Case " << t << ": ";
    	cout << "a maximum of " << G.match().first << " nuts and bolts can be fitted together\n";
    }

    return 0;
}