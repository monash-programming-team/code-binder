// Ford Fulkerson test on SPOJ - FASTFLOW
//
// Expected verdict: Correct but TLE

#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

const ll INF = numeric_limits<ll>::max();

class ford_fulkerson {
    struct edge {
        int a, b;
        ll cap, flow;
    };
	typedef vector<edge>::iterator edge_it;	
    int n; vector<edge> edges; vvi g; vi vis;
	
	ll dfs(int u, int t, ll flow) {
		if (u == t) return flow;
		vis[u] = true;
		for (auto id : g[u]) {
			edge& e = edges[id];
			edge& rev = edges[id^1];
			ll residual = e.cap - e.flow;
			if (vis[e.b] || residual <= 0) continue;
			if (ll augment = dfs(e.b, t, min(flow, residual)) > 0) {
				e.flow += augment;
				rev.flow -= augment;
				return augment;
			}
		}
		return 0;
	}
	
public:
	// Add an edge with capacity cap from node u to node v
	// Return value can be removed if you don't need to examine specific edges
    edge_it add_edge(int u, int v, ll cap){
        edge e1 = {u, v, cap, 0};
        edge e2 = {v, u, 0, 0};	// {v, u, cap, 0} for bidirectional edges
        g[u].push_back( (int) edges.size() );
        edges.push_back(e1);
        g[v].push_back( (int) edges.size() );
        edges.push_back(e2);
		return edges.end() - 2;
    }

	// Initialise a flow network with n nodes
    ford_fulkerson(int n) : n(n), g(n) { }
	
	// Return the max flow from s to t
    ll max_flow(int s, int t) {
        for (auto& e : edges) e.flow = 0;
		ll flow = 0, augment = 0;
		while (vis.assign(n, 0), augment = dfs(s, t, INF) != 0) {
			flow += augment;
		}
		return flow;
    }
};

int main(){
    int n,m;
    scanf("%d%d",&n,&m);
    ford_fulkerson cc(n);

    for(int i=0;i<m;i++){
		int u,v,cap;
		scanf("%d%d%d",&u,&v,&cap); u--; v--;
		cc.add_edge(u, v, cap);
		cc.add_edge(v, u, cap);
    }

    ll mf = cc.max_flow(0, n - 1);
    cout << mf << endl;
    return 0;
}
