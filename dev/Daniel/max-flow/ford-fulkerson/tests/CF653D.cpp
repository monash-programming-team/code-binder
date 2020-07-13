// Test for Ford-Fulkerson on Codeforces 653D
//
// Expected verdict: Accepted

#include<bits/stdc++.h>
using namespace std;

typedef vector<vector<pair<int, int>>> graph;

typedef long double ld;
typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

const ld EPS = 1e-10;
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

int n, m, x;
graph adj;

bool check(ld weight) {
  ld wpb = weight / x;  // weight carried per bear
  ford_fulkerson mf(n);
  for (int v=0; v < n; v++) {
    for (auto e : adj[v]) {
      ll cap = (ll) (e.second / wpb);
      mf.add_edge(v, e.first, cap);
    }
  }
  ll flow = mf.max_flow(0, n-1);
  return flow >= x;
}

int main() {
  int a, b, c;
  cin >> n >> m >> x;

  adj.resize(n);
  for (int i=0; i<m; i++) {
    cin >> a >> b >> c; a--; b--;
    adj[a].emplace_back(b, c);
  }

  ld lo = 0;
  ld hi = 1e12;
  int count = 0;
  while (hi - lo > EPS && count++ < 200) {
    ld mid = (hi + lo) / 2.0;
    if (check(mid)) lo = mid;
    else hi = mid;
  }
  
  cout << fixed << setprecision(10) << lo << endl;
}