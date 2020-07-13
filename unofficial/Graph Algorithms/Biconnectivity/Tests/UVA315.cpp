#include<bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int, int> pii;

class biconnectivity {
	int n;		// Number of nodes in the graph
	vvi edges;	// Edge adjacency list
	
	int dfs_root, dfs_count, root_children;
	vi dfs_low, dfs_num, dfs_parent;
	vi cut_points, is_cut_point;		
	vector<pii> bridge;
	vector<vector<pii>> bccs;
	vector<pii> comp;		

	void dfs(int u) {
		dfs_low[u] = dfs_num[u] = dfs_count++;
		for (int v : edges[u]) {
			pii edge = minmax(u,v);
			if (dfs_num[v] == -1) {
				dfs_parent[v] = u;
				if (u == dfs_root) root_children++;
				comp.push_back(edge);
				dfs(v);
				if (dfs_low[v] >= dfs_num[u]) {
					is_cut_point[u] = true;
					bccs.emplace_back(comp.rbegin(), find(comp.rbegin(),
										comp.rend(),edge)+1);
					comp.resize(comp.size()-bccs.back().size());
				}
				if (dfs_low[v] > dfs_num[u]) bridge.push_back(edge);
				dfs_low[u] = min(dfs_low[u], dfs_low[v]);
			}
			else if (v != dfs_parent[u]) {
				dfs_low[u] = min(dfs_low[u], dfs_num[v]);
				if (dfs_num[v] < dfs_num[u]) comp.push_back(edge);
			}
		}
	}
	
	void search() {
		cut_points.clear(); bridge.clear(); bccs.clear(); comp.clear();
		is_cut_point.assign(n, 0); dfs_parent.assign(n, -1);
		dfs_num.assign(n, -1);  dfs_low.assign(n, 0); dfs_count = 0;
		for (int v = 0; v < n; v++) {
			if (dfs_num[v] == -1) {
				dfs_root = v; root_children = 0;
				dfs(v);
				is_cut_point[v] = (root_children > 1);
			}
		}
		for (int v = 0; v < n; v++) if (is_cut_point[v]) cut_points.push_back(v);
	}
	
public:
	biconnectivity(int n) : n(n), edges(n) { }
	
	void add_edge(int u, int v) { edges[u].push_back(v); edges[v].push_back(u); }

	vi articulation_points() { search(); return move(cut_points); }
	
	vector<pii> bridges() {	search(); return move(bridge);	}

	vector<vector<pii>> biconnected_components() { search(); return move(bccs);	}
	
	tuple<vi, vector<pii>, vector<vector<pii>>> get_all() {
		search();
		return make_tuple(move(cut_points), move(bridge), move(bccs));
	}
};


// ---------------------------------------------------
// Biconnectivity ends. Tests start.
// --------------------------------------------------

int main() {
	int N, u, v;
	while (cin >> N && N > 0) {
		biconnectivity b(N);
		set<pii> edges;
		string stuff;
		while (cin >> u && u > 0) {
			getline(cin, stuff);
			stringstream ss(stuff);
			while (ss >> v) edges.insert(minmax(u-1,v-1));
		}
		for (auto e : edges) b.add_edge(e.first, e.second);
		auto aps = b.articulation_points();
		cout << aps.size() << endl;
	}
}
