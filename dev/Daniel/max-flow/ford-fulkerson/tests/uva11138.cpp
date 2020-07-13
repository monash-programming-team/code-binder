// Ford-Fulkerson test on UVA 11138
// Expected verdict: Accepted

#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

const ll INF = numeric_limits<ll>::max();

struct ford_fulkerson {
    struct edge {
        int from, to;
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
        edges.push_back({u, v, 0, cap});
        g[v].push_back( (int) edges.size() );
        edges.push_back({v, u, 0, 0});  // {u, 0, cap} for bidirectional edges
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

int SRC = 0;
int TGT = 1;

int n, b;

int bolt(int i) { return 2 + i; }
int nut(int i) { return 2 + b + i; }

int main(){
    
    int T;
    cin >> T;
    
    for (int t = 1; t <= T; t++) {
    	cin >> b >> n;
    	
    	ford_fulkerson cc(2 + n + b);
    	
    	int i;
    	for (int c = 0; c < b; c++) {
    		for (int d = 0; d < n; d++) {
    			cin >> i;
    			if (i) cc.add_edge(bolt(c), nut(d), 1);
    		}
    	}
    	
    	for (int d = 0; d < n; d++) cc.add_edge(nut(d), TGT, 1);
    	for (int c = 0; c < b; c++) cc.add_edge(SRC, bolt(c), 1);

    	cout << "Case " << t << ": ";
    	cout << "a maximum of " << cc.max_flow(SRC, TGT) << " nuts and bolts can be fitted together\n";

    	
    }

    return 0;
}
