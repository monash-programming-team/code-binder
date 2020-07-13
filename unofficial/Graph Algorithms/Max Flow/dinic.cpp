// Dinic's algorithm for maximum flow
//
// Author      : Daniel Anderson
// Date        : September 1, 2016
// Reliability : 5
// Tested on   : CF - 589F
//				       SPOJ - FASTFLOW
//				       UVA - 820 (Internet Bandwidth)
//				       CF - 653D
//				       NEERC13 Eastern Subregional - H
//
// TLEs On:	 Tuesday Tutorial: Network Flow (TLE on test 57)
//           SPOJ - MATCHING (TLE)
//
// Usage:
//    ll max_flow(int s, int t)
//         Returns the maximum flow s -> t
//
//	  int add_edge(u, v, cap)
//			Adds an edge from u -> v with capacity. Returns the index
//			of the edge. 
//
//	  edge& get_edge(i)
//			Returns a reference to the i'th edge. Use to check flows
//			or update capacities
//
//	Time Complexity: O(V^2 E)

#include<bits/stdc++.h>
using namespace std;

typedef long long int ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:dinic
// Dinic's algorithm for maximum flow. add_edge returns the id of an edge which can be
// used to inspect the final flow value using get_edge(i).flow. Complexity: O(V^2 E)
template<typename T> struct Dinics {
	struct edge { int to; T flow, cap; };  	T INF = numeric_limits<T>::max();
	int n, m;  vi dist, work;  queue<int> q;	vector<edge> edges;  vvi adj;
	bool bfs(int s, int t) {
		dist.assign(n, -1);  dist[s] = 0;  q.push(s);
		while (!q.empty()) {
			int u = q.front(); q.pop();
			for (auto& i : adj[u]) {
				edge& e = edges[i];  int v = e.to;
				if (dist[v] < 0 && e.flow < e.cap) dist[v] = dist[u] + 1, q.push(v);
			}
		}
		return dist[t] >= 0;
	}
	T dfs(int u, int t, T f) {
		if (u == t) return f;
		for (int& i = work[u]; i < (int)adj[u].size(); i++) {
			auto& e = edges[adj[u][i]], &rev = edges[adj[u][i]^1];
			if (e.flow < e.cap && dist[e.to] == dist[u] + 1) {
				T df = dfs(e.to, t, min(f, e.cap - e.flow));
				if (df > 0) { e.flow += df; rev.flow -= df; return df; }
			}
		}
		return 0;
	}
	// Create a flow network with n nodes -- add edges with add_edge(u,v,cap)
	Dinics(int n) : n(n), m(0), adj(n) { }
	int add_edge(int from, int to, T cap) {        // add an edge (from -> to) with
		adj[from].push_back(m++), adj[to].push_back(m++);   // capacity of cap units.
    edges.push_back({to, 0, cap}), edges.push_back({from, 0, 0});
		return m - 2;    // Change {from,0,0} to {from,0,cap} for bidirectional edges
	}
	edge& get_edge(int i) { return edges[i]; }  // get a reference to the i'th edge
	T max_flow(int s, int t) {                 // find the maximum flow from s to t
		T res = 0;  for (auto& e : edges) e.flow = 0;
		while (work.assign(n, 0), bfs(s, t))
			while (T delta = dfs(s, t, INF)) res += delta;
		return res;
	}
};
//listings:/dinic

int main() {
	int n,m,s,t;
	scanf("%d%d%d%d",&n,&m,&s,&t);
	s--; t--;
	Dinics<ll> cc(n);
	for(int i=0;i<m;i++){
	  int u,v,c;
	  scanf("%d%d%d",&u,&v,&c);
	  u--; v--;
	  cc.add_edge(u,v,c);
	}
	cout << cc.max_flow(s, t) << endl;
}
