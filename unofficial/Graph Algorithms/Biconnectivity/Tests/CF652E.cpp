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
	vector<pii> bridge;	

	void dfs(int u) {
		dfs_low[u] = dfs_num[u] = dfs_count++;
		for (int v : edges[u]) {
			pii edge = minmax(u,v);
			if (dfs_num[v] == -1) {
				dfs_parent[v] = u;
				if (u == dfs_root) root_children++;
				dfs(v);
				if (dfs_low[v] > dfs_num[u]) bridge.push_back(edge);
				dfs_low[u] = min(dfs_low[u], dfs_low[v]);
			}
			else if (v != dfs_parent[u]) {
				dfs_low[u] = min(dfs_low[u], dfs_num[v]);
			}
		}
	}
	
	void search() {
		bridge.clear();
		dfs_parent.assign(n, -1);
		dfs_num.assign(n, -1);  dfs_low.assign(n, 0); dfs_count = 0;
		for (int v = 0; v < n; v++) {
			if (dfs_num[v] == -1) {
				dfs_root = v; root_children = 0;
				dfs(v);
			}
		}
	}
	
public:
	biconnectivity(int n) : n(n), edges(n) { }
	void add_edge(int u, int v) { edges[u].push_back(v); edges[v].push_back(u); }
	vector<pii> bridges() {	search(); return move(bridge);	}

};



vvi graph;
vvi cc_graph;
map<pii, bool> edge_has_artifact;
vector<bool> cc_has_artifact;
map<pii, bool> cc_edge_has_artifact;
set<pii> bridges;
vvi connected_components;
vi island_cc;

bool is_bridge(int u, int v) {
	return bridges.find(minmax(u,v)) != bridges.end();
}

bool dfs(int v, vi& cc, int cc_num) {
	island_cc[v] = cc_num;
	cc.push_back(v);
	bool has = false;
	for (int child : graph[v]) {
		if (is_bridge(v, child)) continue;
		has = has || edge_has_artifact[minmax(v,child)];
		if (island_cc[child] == -1) has = dfs(child, cc, cc_num) || has;
	}
	return has;
}

bool dfs2(int a, int b, int parent, vi& path) {
	if (a == b) {
		path.push_back(b);
		return true;
	}
	for (int child : cc_graph[a]) {
		if (child != parent && dfs2(child, b, a, path)) {
			path.push_back(a);
			return true;
		}
	}
	return false;
}

int main() {
	ios_base::sync_with_stdio(false);
	cin.tie(0); cout.tie(0);
	
	int n, m, x, y, z, a, b;
	cin >> n >> m;
	
	graph.resize(n);
	biconnectivity bi(n);
	for (int i=0; i<m; i++) {
		cin >> x >> y >> z;
		x--; y--;
		graph[x].push_back(y);
		graph[y].push_back(x);
		bi.add_edge(x,y);
		edge_has_artifact[minmax(x,y)] = z;
	}
	
	cin >> a >> b; a--; b--;
	
	auto bridge = bi.bridges();
	for (auto& br : bridge) {
		if (br.second < br.first) swap(br.first, br.second);
		bridges.insert(br);
	}
	
	int cc_num = 0;
	island_cc.assign(n, -1);
	for (int v=0; v<n; v++) {
		if (island_cc[v] == -1) {
			vi cc;
			bool has = dfs(v, cc, cc_num++);
			connected_components.push_back(cc);
			cc_has_artifact.push_back(has);
		}
	}
	
	cc_graph.resize(cc_num);
	for (const auto& br : bridge) {
		int u = island_cc[br.first];
		int v = island_cc[br.second];
		cc_graph[u].push_back(v);
		cc_graph[v].push_back(u);
		cc_edge_has_artifact[minmax(u,v)] = edge_has_artifact[minmax(br.first, br.second)];
	}
	
	if (island_cc[a] == island_cc[b] && cc_has_artifact[island_cc[a]]) {
		cout << "YES" << endl;
		return 0;
	}
	else {
		vi path;
		dfs2(island_cc[a], island_cc[b], -1, path);
		path.push_back(island_cc[a]);
		for (int i=0; i<(int)path.size()-1; i++) {
			if (cc_has_artifact[path[i]] || cc_has_artifact[path[i+1]] ||
				cc_edge_has_artifact[minmax(path[i], path[i+1])]) {
					cout << "YES" << endl;
					return 0;
				}
		}
		cout << "NO" << endl;
		return 0;
	}
}
