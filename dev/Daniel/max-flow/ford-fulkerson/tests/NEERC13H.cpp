// Ford Fulkerson test on NEERC13 Eastern Subregional - Problem H
// Expected verdict: Accepted
//

#include <bits/stdc++.h>
using namespace std;;

#define debug(x) cerr << #x << " = " << x << endl;

typedef long long int ll;
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<vi> vvi;

const ll INF = numeric_limits<ll>::max();

class ford_fulkerson {
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
			ll residual = e.cap - e.flow;
			if (vis[e.to] || residual <= 0) continue;
			if (ll augment = dfs(e.to, t, min(flow, residual)) > 0) {
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
		while (vis.assign(n, 0), augment = dfs(s, t, INF) != 0) {
			flow += augment;
		}
		return flow;
    }
};


int a, b, n;
vector<pii> input;
vi t, d;

int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);
	
	cin >> a >> b >> n;
	t.resize(n); d.resize(n);
	input.resize(n);
	
	for (int i = 0; i < n; i++) {
		cin >> input[i].first >> input[i].second;
	}
	
	sort(input.begin(), input.end());
	
	for (int i = 0; i < n; i++) {
		t[i] = input[i].first;
		d[i] = input[i].second;
	}
	
	vi opening_index, closing_index;
	vi opening_times, closing_times;
	
	for (int i = 0; i < n; i++) {
		if (d[i] == 0) {
			opening_index.push_back(i);
			opening_times.push_back(t[i]);
		}
		else {
			closing_index.push_back(i);
			closing_times.push_back(t[i]);
		}
	}

	int u = opening_times.size();
	int v = closing_times.size();

	if (u != v) {
		cout << "Liar" << endl;
		return 0;
	}

	// Build the graph
	ford_fulkerson G(2 + u + v);
	
	for (int i = 0; i < u; i++) {
		G.add_edge(0, 2 + i, 1);
	}
	
	for (int i = 0; i < v; i++) {
		G.add_edge(2 + u + i, 1, 1);
	}

	vi pairs;

	for (int i = 0; i < v; i++) {
		int time = closing_times[i];
		
		// Add smuggler edges
		auto smuggler_it = upper_bound(opening_times.begin(), opening_times.end(), time - a);
		for (auto it = opening_times.begin(); it != smuggler_it; it++) {
			int j = distance(opening_times.begin(), it);
			pairs.push_back(G.add_edge(2 + j, 2 + u + i, 1));
		}
		
		// Add cargo edges
		auto cargo_it = lower_bound(opening_times.begin(), opening_times.end(), time - b);
		auto cargo_end = upper_bound(opening_times.begin(), opening_times.end(), time);
		for (auto it = cargo_it; it != cargo_end; it++) {
			int j = distance(opening_times.begin(), it);
			pairs.push_back(G.add_edge(2 + j, 2 + u + i, 1));
		}
	}
	
	// Do matching
	int numMatches = G.max_flow(0,1);
	
	// Answer
	if (numMatches != u) {
		cout << "Liar" << endl;
		return 0;
	}
	else {
		cout << "No reason" << endl;
		for (auto i : pairs) {
			auto& e = G.get_edge(i);
			if (e.flow == 1) {
				int s = e.from - 2;
				int t = e.to - 2 - u;
				int start = opening_times[s];
				int end = closing_times[t];
				cout << start << ' ' << end << '\n';
			}
		}
	}

}