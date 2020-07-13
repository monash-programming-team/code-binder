// Dinic test on Codeforces 589F
//
// Expected verdict: Accepted
//
// Good demonstration on how to use the edge indices
// to modify capacities

#include <bits/stdc++.h>
using namespace std;

#define debug(x) cerr << #x << " = " << x << endl;

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

template<typename T>
class dinic {
	struct edge {
		int to;
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
		edges.push_back({to, 0, cap});
		adj[to].push_back((int)edges.size());
		edges.push_back({from, 0, 0});	// change to (from, 0, cap) for bidirectional
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

int n, num_nodes, max_t;
vi a, b;

const int SRC = 0;
const int TGT = 1;

vi food_edges;

int food_node(int i) { return 2 + i; }
int time_node(int t) { return 2 + n + t; }

// True if all dishes can be eaten for T seconds
bool check(int T, dinic<int>& cc) {
	for (auto i : food_edges) cc.get_edge(i).cap = T;
	int flow = cc.max_flow(SRC, TGT);
	return (flow == n * T);
}

int main(){
	cin >> n;
	a.resize(n); b.resize(n);
	for (int i = 0; i < n; i++) cin >> a[i] >> b[i];
	max_t = *max_element(b.begin(), b.end());
	num_nodes = 2 + n + max_t;
	
    dinic<int> cc(num_nodes);
	
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
