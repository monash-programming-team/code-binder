// Minimum cost flow using Leviticus algorithm for cheapest augmenting paths.
//
// The theoretical complexity isn't great but the algorithm is extremely fast
// in practice.
//
// Author      : Daniel Anderson
// Date        : 08-09-2016
// Reliability : 1
// Tested on   : UVA 10594 (Data Flow)
//
// Reference: http://e-maxx.ru/algo/min_cost_flow
//
// Usage:
//    ll flow(int s, int t, ll K)
//   	 Find the minimum cost to send K units of flow from s to t
// 		 If it is not possible to send K units of flow from s to t
// 		 then we will find the minimum cost to send the maximum
// 		 amount of flow possible from s to t. If you just want
// 		 min-cost max-flow then use K = LLONG_MAX
//
// 		Returns: {max_flow, min_cost}
//
//	  int add_edge(u, v, cap, cost)
//			Adds a directed edge from u -> v with capacity and cost. Returns the index
//			of the edge. If you want bidirectional edges, add both (u, v) and (v, u)
//
//	  edge& get_edge(i)
//			Returns a reference to the i'th edge. Use to check flows
//			or update capacities / costs.
//
//	Time Complexity: O(V^2 E^2) theoretically.

//  Runs MUCH faster in practice. Has been tested and runs in 0.1 seconds on
//  multi-testcase problems with V = 100, M = 5000 and cost <= 10^15.

#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<ll> vll;

const ll INF = numeric_limits<ll>::max();

class min_cost_flow {
	struct edge {
        int to;
        ll flow, cap, cost;
    };
	int n; vector<edge> edges; vvi g;
	vi pred, pred_edge;

	bool levit(int s, int t) {
		vi id(n, 0);
		vll dist(n, INF);  dist[s] = 0;
		deque<int> q;  q.push_back(s);
		while (!q.empty()) {
			int v = q.front(); q.pop_front();
			id[v] = 2;
			for (auto i : g[v]) {
				auto& e = edges[i];
				if (e.flow < e.cap && dist[v] + e.cost < dist[e.to]) {
					dist[e.to] = dist[v] + e.cost;
					if (id[e.to] == 0) q.push_back(e.to);
					else if (id[e.to] == 2) q.push_front(e.to);
					id[e.to] = 1;
					pred[e.to] = v;
					pred_edge[e.to] = i;
				}
			}	
		}
		return dist[t] < INF;
	}
	
	pair<ll, ll> augment(int s, int t, ll cap) {
		ll flow = cap, cost = 0;
		for (int v = t; v != s; v = pred[v]) {
			int pe = pred_edge[v];
			flow = min(flow, edges[pe].cap - edges[pe].flow);
		}
		for (int v = t; v != s; v = pred[v]) {
			int pe = pred_edge[v], rev = pe ^ 1;
			edges[pe].flow += flow;
			edges[rev].flow -= flow;
			cost += edges[pe].cost * flow;
		}
		return {flow, cost};
	}
	
 public:
 	// Initialise a flow network with n nodes
 	min_cost_flow(int n) : n(n), g(n), pred(n), pred_edge(n) { }
 	
	// Add a directed edge with capacity and cost from node u to node v
	// Returns the index of the edge.
    int add_edge(int u, int v, ll cap, ll cost) {
        g[u].push_back( (int) edges.size() );
        edges.push_back({v, 0, cap, cost});
        g[v].push_back( (int) edges.size() );
        edges.push_back({u, 0, 0, -cost});
		return (int)edges.size() - 2;
    }
	
	// Returns a reference to the given edge. Useful for modifying
	// capacities and costs or inspecting final flow values.
	edge& get_edge(int i) { return edges[i]; }
	
	// Find the minimum cost to send K units of flow from s to t
	pair<ll, ll> flow(int s, int t, ll K) {
		for (auto& e : edges) e.flow = 0;
		ll flow = 0, cost = 0;
		while (flow < K && levit(s, t)) {
			auto res = augment(s, t, K - flow);
			flow += res.first;
			cost += res.second;
		}
		return {flow, cost};
	}
};


int main(){
    int N, M;
    
    while (scanf("%d %d", &N, &M) != EOF) {
		min_cost_flow mcf(N);
		
		vi us(M), vs(M);
		vll costs(M);
		
		for (int i = 0; i < M; i++) {
			scanf("%d %d %lld", &us[i], &vs[i], &costs[i]);
		}
		
		ll D, K;
		scanf("%lld %lld", &D, &K);
		
		for (int i = 0; i < M; i++) {
			mcf.add_edge(us[i] - 1, vs[i] - 1, K, costs[i]);
			mcf.add_edge(vs[i] - 1, us[i] - 1, K, costs[i]);
		}
		
		auto ans = mcf.flow(0, N - 1, D);
		
		if (ans.first < D) puts("Impossible.");
		else printf("%lld\n", ans.second);
    }
}
