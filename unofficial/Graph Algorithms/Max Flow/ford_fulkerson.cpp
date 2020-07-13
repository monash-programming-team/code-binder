// Ford-Fulkerson maximum flow using dfs augmenting paths
//
// Author      : Daniel Anderson
// Date        : 28-08-2016
// Reliability : 5
// Tested on   : UVA - 820 (Internet Bandwidth)
//				 CF - 653D
//				 NEERC13 Eastern Subregional H
//				 UVA - 11138 (Nuts and bolts)
//				 UVA - 10092 (The Problem with the Problem Setter)
//
// Tests on large input (AC until TLE)
//  TLE: 	
//          Tuesday Tutorial: Network Flow (TLE on test 21)
//			SPOJ - FASTFLOW	
//			CF - 589F
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
//	Time Complexity: O(Ef), where f is the maximum flow

#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

const ll INF = numeric_limits<ll>::max();

class ford_fulkerson {
    struct edge {
        int to;
        ll flow, cap;
    };	
    int n; vector<edge> edges; vvi g; vi vis;
	
	ll dfs(int u, int t, ll flow) {
		if (u == t) return flow;
		vis[u] = true;
		for (auto id : g[u]) {
			edge& e = edges[id];
			edge& rev = edges[id^1];
			ll residual = e.cap - e.flow, augment = 0;
			if (vis[e.to] || residual <= 0) continue;
			if ((augment = dfs(e.to, t, min(flow, residual))) > 0) {
				e.flow += augment;
				rev.flow -= augment;
				return augment;
			}
		}
		return 0;
	}
	
public:
	// Initialise a flow network with n nodes
    ford_fulkerson(int n) : n(n), g(n) { }

	// Add an edge with capacity cap from node u to node v
	// Returns the index of the edge.
    int add_edge(int u, int v, ll cap){
        g[u].push_back( (int) edges.size() );
        edges.push_back({v, 0, cap});
        g[v].push_back( (int) edges.size() );
        edges.push_back({u, 0, 0});  // {u, 0, cap} for bidirectional edges
		return (int)edges.size() - 2;
    }

	// Get a reference to a specific edge: use to check flows or update capcities
	edge& get_edge(int i) { return edges[i]; }
	
	// Return the max flow from s to t
    ll max_flow(int s, int t) {
        for (auto& e : edges) e.flow = 0;
		ll flow = 0, augment = 0;
		while (vis.assign(n, 0), (augment = dfs(s, t, INF)) != 0) {
			flow += augment;
		}
		return flow;
    }
};

int main(){
    int n,m,s,t;
    scanf("%d%d%d%d",&n,&m,&s,&t);
    ford_fulkerson cc(n);

    for(int i=0;i<m;i++){
		int u,v,cap;
		scanf("%d%d%d",&u,&v,&cap); u--; v--;
		cc.add_edge(u, v, cap);
    }

    ll mf = cc.max_flow(s-1, t-1);
    cout << mf << endl;
    return 0;
}
