// Dijkstra's Algorithm
//
// Author      : Daniel Anderson
// Date        : April 8, 2016
// Reliability : 5
// Tested On   : SPOJ - SHPATH
//				 UVA - 10986
//				 Codeforces 602C
//				 Codeforces 20C
//				 UVA 10389
//
// Computes the shortest distance between a source vertex and all
// other vertices in a weighted graph (directed or undirected).
// 
// All edge weights must be non-negative.
//
// Complexity: O( E log E )
#include<bits/stdc++.h>
using namespace std;

typedef long long int ll;
typedef vector<int> vi;

//listings:dijkstra
// Shortest paths in non-negative weighted graphs. dist[u] = INF if u is not
// reachable. T is the type of the edge weights / costs. Complexity: O(E log(E))
template<typename T> struct dijkstra {
	typedef pair<T, int> pti;  vector<vector<pti>> edges;
	int n;  const T INF = numeric_limits<T>::max() / 2;
	dijkstra(int n) : n(n), edges(n) { }
	void add_edge(int u, int v, T weight) { edges[u].emplace_back(weight, v); }
	pair<vector<T>, vi> shortest_paths(int src) {
		vector<T> dist(n, INF); vi pred(n, -1); dist[src] = 0;
    priority_queue<pti, vector<pti>, greater<pti>> q; q.emplace(0, src);
		int u, v; T d, w;
		while (!q.empty()) {
			tie(d, u) = q.top(); q.pop();
			if (dist[u] < d) continue;
			for (auto& e : edges[u]) {
				tie(w, v) = e;
				if (dist[u] + w < dist[v]) 
          dist[v] = dist[u] + w, pred[v] = u,	q.emplace(dist[v], v);
			}
		}
		return {dist, pred};
	}
};
//listings:/dijkstra

//listings:get_path
// Reconstruct the path corresponding to pred from Dijkstra and Bellman-Ford
vi get_path(int v, vi& pred) {
	vi p = {v};
	while (pred[v] != -1) p.push_back(v = pred[v]);
	reverse(p.begin(), p.end());
	return p;
}
//listings:/get_path

// ---------------------------------------------------
// Dijkstra ends. Example: Codeforces 20C
// --------------------------------------------------

int main() {
	int n, m;
	cin >> n >> m;
	
	//Initialise the graph
	dijkstra<ll> g(n);
	
	//Add the edges
	for (int i=0; i<m; i++) {
		int a, b, w;
		cin >> a >> b >> w;
		a--; b--;
		g.add_edge(a, b, w);
		g.add_edge(b, a, w);
	}
	
	//Find the shortest paths from 0
	auto ans = g.shortest_paths(0);
	
	//Check for discontinuity
	if (ans.first[n-1] >= g.INF) cout << -1 << endl;
	
	//Recover the path
	else {
		vi path = get_path(n-1, ans.second);
		for (int p : path) cout << p + 1 << ' '; cout << endl;
	}
	
	return 0;
}
