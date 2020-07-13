// Maximum bipartite matching via alternating paths.
// Theoretically the slow algorithm but very fast in practice.
//
// Author      : Daniel Anderson
// Date        : 5/10/2016
// Reliability : 3
// Tested on   : UVA10092, UVA11138, UVA10080
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
	
	int dfs(int v) {
		vis[v] = p;
		for (auto u : adj[v]) 
			if (m[u] < 0 || (vis[m[u]] != p && dfs(m[u])))
				return m[u] = v, 1;
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

string get_operator(ll a, ll b, ll x) {
	if (a + b == x) return " + ";
	else if (a - b == x) return " - ";
	else return " * ";
}

int main() {
	int n, m;
	cin >> n;
	vector<ll> a(n), b(n), RHS;
	for (int i = 0; i < n; i++) {
		cin >> a[i] >> b[i];
		RHS.push_back(a[i] + b[i]);
		RHS.push_back(a[i] - b[i]);
		RHS.push_back(a[i] * b[i]);
	}
	
	sort(RHS.begin(), RHS.end());
	RHS.erase(unique(RHS.begin(), RHS.end()), RHS.end());
	m = (int)RHS.size();
	
	map<ll, int> index;
	for (int i = 0; i < m; i++) {
		index[RHS[i]] = i;
	}
	
	max_matching G(n, m);
	for (int i = 0; i < n; i++) {
		G.add_edge(i, index[a[i] + b[i]]);
		G.add_edge(i, index[a[i] - b[i]]);
		G.add_edge(i, index[a[i] * b[i]]);
	}
	
	auto res = G.match();
	if (res.first < n) cout << "impossible" << endl;
	else {
		auto& X = res.second;
		for (int i = 0; i < n; i++) {
			cout << a[i] << get_operator(a[i], b[i], RHS[X[i]]) << b[i] << " = " << RHS[X[i]] << '\n';
		}
	}
}