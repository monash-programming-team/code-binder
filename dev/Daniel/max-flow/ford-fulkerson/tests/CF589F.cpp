// Ford-Fulkerson test on Codeforces 589F
//
// Expected verdict: TLE
//
// Good demonstration on how to use the edge iterators
// to modify / check edge flows / capacities.

#include <bits/stdc++.h>
using namespace std;

#define debug(x) cerr << #x << " = " << x << endl;

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

const ll INF = numeric_limits<ll>::max();

class ford_fulkerson {
    struct edge {
        int a, b;
        ll cap, flow;
    };	
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
    int add_edge(int u, int v, ll cap){
        edge e1 = {u, v, cap, 0};
        edge e2 = {v, u, 0, 0};	// {v, u, cap, 0} for bidirectional edges
        g[u].push_back( (int) edges.size() );
        edges.push_back(e1);
        g[v].push_back( (int) edges.size() );
        edges.push_back(e2);
		return (int)edges.size() - 2;
    }

	// Initialise a flow network with n nodes
    ford_fulkerson(int n) : n(n), g(n) { }

	// Get a reference to a specific edge
	edge& get_edge(int i) { return edges[i]; }
	
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

int n, num_nodes, max_t;
vi a, b;

const int SRC = 0;
const int TGT = 1;

int food_node(int i) { return 2 + i; }
int time_node(int t) { return 2 + n + t; }

vi food_edges;

// True if all dishes can be eaten for T seconds
bool check(int T, ford_fulkerson& cc) {
	for (auto i : food_edges) cc.get_edge(i).cap = T;
	ll flow = cc.max_flow(SRC, TGT);
	return (flow == n * T);
}

int main(){
	cin >> n;
	a.resize(n); b.resize(n);
	for (int i = 0; i < n; i++) cin >> a[i] >> b[i];
	max_t = *max_element(b.begin(), b.end());
	num_nodes = 2 + n + max_t;
	
    ford_fulkerson cc(num_nodes);
	
	// Connect source to food
	for (int i = 0; i < n; i++) {
		food_edges.push_back(cc.add_edge(SRC, food_node(i), 0));
	}
	
	// Connect times to sink
	for (int t = 0; t < max_t; t++) {
		cc.add_edge(time_node(t), TGT, 1);
	}
    
	// Connect food to times
	for (int i = 0; i < n; i++) {
		for (int t = a[i]; t < b[i]; t++) {
			cc.add_edge(food_node(i), time_node(t), 1);
		}
	}
	
	// Binary search on maximum time
	int lo = 0, hi = 20000;
	while (lo != hi - 1) {
		int T = (lo + hi) / 2;
		if (check(T, cc)) lo = T;
		else hi = T;
	}
	
	cout << lo * n << endl;
	return 0;
}
