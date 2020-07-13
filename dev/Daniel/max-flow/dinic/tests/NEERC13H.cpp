// Dinic test on NEERC13 H
//
// Verdict: Accepted

#include <bits/stdc++.h>
using namespace std;;

#define debug(x) cerr << #x << " = " << x << endl;

typedef long long int ll;
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<vi> vvi;

// Dinic's algorithm for maximum flow
//
// Parameters:
// 	T is the type of your flows / capacities
template<typename T>
class dinic {
	struct edge {
		int from, to;
		T flow, cap;
	};
	
	T INF = numeric_limits<T>::max();

	int n;
	vi dist, work;
	vector<edge> edges;
	vvi adj;
	queue<int> q;

	bool bfs(int s, int t) {
		dist.assign(n, -1);
		dist[s] = 0;
		q.push(s);
		while (!q.empty()) {
			int u = q.front(); q.pop();
			for (auto& i : adj[u]) {
				edge& e = edges[i];
				int v = e.to;
				if (dist[v] < 0 && e.flow < e.cap) {
					dist[v] = dist[u] + 1;
					q.push(v);
				}
			}
		}
		return dist[t] >= 0;
	}

	T dfs(int u, int t, T f) {
		if (u == t) return f;
		for (int& i = work[u]; i < (int)adj[u].size(); i++) {
			edge& e = edges[adj[u][i]];
			edge& rev = edges[adj[u][i]^1];
			if (e.cap <= e.flow) continue;
			int v = e.to;
			if (dist[v] == dist[u] + 1) {
				T df = dfs(v, t, min(f, e.cap - e.flow));
				if (df > 0) {
					e.flow += df;
					rev.flow -= df;
					return df;
				}
			}
		}
		return 0;
	}
	

public:
	// Create a flow network with n vertices
	dinic(int nn) : n(nn), adj(n) { }

	// Add an edge with capacity cap between the given edges
	// Can be modified to use bidirectional edges
	// Returns the index of the edge
	int add_edge(int from, int to, T cap) {
		adj[from].push_back((int)edges.size());
		edges.push_back({from, to, 0, cap});
		adj[to].push_back((int)edges.size());
		edges.push_back({to, from, 0, 0});	// change to (from, 0, cap) for bidirectional
		return (int)edges.size() - 2;
	}
	
	// Returns a reference to the given edge. Useful for modifying
	// capacities or inspecting final flow values.
	edge& get_edge(int i) { return edges[i]; }

	// Finds the max flow from s to t
	T max_flow(int s, int t) {
		T res = 0;
		for (auto& e : edges) e.flow = 0;
		while (bfs(s, t)) {
			work.assign(n, 0);
			while (T delta = dfs(s, t, INF))
				res += delta;
		}
		return res;
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
	dinic<int> G(2 + u + v);
	
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
