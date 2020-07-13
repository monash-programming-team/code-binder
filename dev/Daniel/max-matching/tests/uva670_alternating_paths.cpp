// Maximum bipartite matching via alternating paths.
// Theoretically the slow algorithm but very fast in practice.
//
// Author      : Daniel Anderson
// Date        : 5/10/2016
// Reliability : 5
// Tested on   : UVA10092, UVA11138, UVA10080, ANZAC2016 2E, UVA670
//
// Usage:
//  max_matching G(l, r)    create a graph with l left and r right nodes
//  add_edge(u, v)          add an edge from left node u to right node v
//  match()                 returns the size of the max matching and the matches
//
// Complexity: O(nm)
#include<bits/stdc++.h>
using namespace std;

typedef long long ll;
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

const double EPS = 1e-8;

typedef pair<int, int> pii;
#define X first
#define Y second

#define sqr(x) (x)*(x)
double dist(pii a, pii b) {
	return sqrt(sqr(a.X - b.X) + sqr(a.Y - b.Y));
}

int main() {
	int L;
	cin >> L;
	while (L--) {
		int N, M;
		cin >> N >> M;
		vector<pii> bob(N);
		for (int i = 0; i < N; i++) {
			cin >> bob[i].X >> bob[i].Y;
		}
		
		vector<pii> place(M);
		for (int i = 0; i < M; i++) {
			cin >> place[i].X >> place[i].Y;
		}
		
		max_matching G(N - 1, M);
		
		for (int i = 0; i < N - 1; i++) {
			for (int j = 0; j < M; j++) {
				if (dist(bob[i], place[j]) + dist(place[j], bob[i+1]) <= 2 * dist(bob[i], bob[i+1])) {
					G.add_edge(i, j);
				}
			}
		}
		
		auto res = G.match();
		
		vector<pii> route;
		for (int i = 0; i < N - 1; i++) {
			route.push_back(bob[i]);
			if (res.second[i] != -1) {
				assert(dist(bob[i], place[res.second[i]]) + dist(place[res.second[i]], bob[i+1]) <= 2 * dist(bob[i], bob[i+1]));
				route.push_back(place[res.second[i]]);
			}
		}
		route.push_back(bob.back());
		assert((int)route.size() == N + res.first);
		
		cout << N + res.first << '\n';
		for (int i = 0; i < (int)route.size(); i++) {
			if (i) cout << ' ';
			cout << route[i].X << ' ' << route[i].Y;
		}
		
		cout << '\n';
		if (L) cout << '\n';
	}
}