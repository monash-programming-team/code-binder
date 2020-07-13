// Alternating paths algorithm test
// 
// Verdict: AC
//
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

const double EPS = 1e-8;

#define SQR(x) (x)*(x)
double dist(double x, double y, pair<double,double> gopher) {
	return sqrt(SQR(x - gopher.first) + SQR(y - gopher.second));
}

int main() {
	int n, m, s, v;
	while (cin >> n >> m >> s >> v) {
		max_matching G(n, m);
		vector<pair<double,double>> gophers(n);
		for (int i = 0; i < n; i++) {
			cin >> gophers[i].first >> gophers[i].second;
		}
		for (int j = 0; j < m; j++) {
			double x, y;
			cin >> x >> y;
			for (int i = 0; i < n; i++) {
				if (dist(x, y, gophers[i]) / v + EPS < s) {
					G.add_edge(i, j);
				}
			}
		}
		cout << (n - G.match().first) << endl;
	}
}