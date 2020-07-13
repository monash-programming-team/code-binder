// Test alternating paths algorithm on UVA10092
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

int main() {
	int nk, np;
	while (cin >> nk >> np && nk) {
		vi category;
		for (int i = 0; i < nk; i++) {
			int num; cin >> num;
			for (int j = 0; j < num; j++)
				category.push_back(i);
		}
		int R = (int)category.size();
		
		max_matching G(np, R);
		
		for (int j = 0; j < np; j++) {
			int num; cin >> num;
			for (int i = 0; i < num; i++) {
				int cat; cin >> cat; cat--;
				for (int k = 0; k < R; k++) {
					if (category[k] == cat)
						G.add_edge(j, k);
				}
			}
		}
		
		auto res = G.match();
		if (res.first < R) cout << 0 << endl;
		else {
			vvi allocation(nk);
			for (int j = 0; j < np; j++) {
				if (res.second[j] != -1)
					allocation[category[res.second[j]]].push_back(j);
			}
			
			cout << 1 << endl;
			for (int i = 0; i < nk; i++) {
				for (int j = 0; j < (int)allocation[i].size(); j++) {
					cout << (allocation[i][j] + 1) << " \n"[j == (int)allocation[i].size() - 1];
				}
			}
		}
	}
}