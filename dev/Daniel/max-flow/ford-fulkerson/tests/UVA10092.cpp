// Ford-Fulkerson test on UVA10092
// Expected Verdict: Accepted
//

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

int nk, np;
const int SRC = 0;
const int TGT = 1;

int category(int i) { return 2 + i; }
int problem(int i) { return 2 + nk + i; }

int main(){
	
	while (cin >> nk >> np, nk != 0) {
		
		ford_fulkerson cc(2 + nk + np);

		vi requirements(nk);
		for (int i = 0; i < nk; i++)
			cin >> requirements[i];
		
		vvi edges(nk);
		
		// Add edges from source to categories
		for (int i = 0; i < nk; i++) {
			cc.add_edge(SRC, category(i), requirements[i]);
		}
		
		// Add edges from problems to sink
		for (int i = 0; i < np; i++) {
			cc.add_edge(problem(i), TGT, 1);
		}
		
		// Add edges from categories to problems
		for (int i = 0; i < np; i++) {
			int n, cat;
			cin >> n;
			for (int j = 0; j < n; j++) {
				cin >> cat; cat--;
				edges[cat].push_back(cc.add_edge(category(cat), problem(i), 1));
			}
		}
		
		// Do flow
		auto flow = cc.max_flow(SRC, TGT);
		
		// Output the solution
		if (flow == accumulate(requirements.begin(), requirements.end(), 0)) {
			cout << 1 << '\n';
			for (int i = 0; i < nk; i++) {
				int num_printed = 0;
				for (int j = 0; j < (int)edges[i].size(); j++) {
					auto& e = cc.get_edge(edges[i][j]);
					if (e.flow == 1) {
						num_printed++;
						cout << (e.to - 2 - nk + 1) << " \n"[num_printed == requirements[i]];
					}
				}
			}
		}
		else {
			cout << 0 << '\n';
		}
	}
}
