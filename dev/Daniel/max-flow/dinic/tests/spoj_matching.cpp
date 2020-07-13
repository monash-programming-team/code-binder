// Dinic's test on spoj matching
// verdict: TLE

#include<bits/stdc++.h>
using namespace std;

typedef long long int ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

// Dinic's algorithm for maximum flow
//
// Parameters:
// 	T is the type of your flows / capacities
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

int N, M, P, A, B;

int cow(int i) {
	return 2 + i;
}

int bull(int i) {
	return 2 + N + i;
}

int main() {

	scanf("%d%d%d", &N, &M, &P);
	dinic<int> cc(N + M + 2);
	
	for (int i = 0; i < N; i++) {
		cc.add_edge(0, cow(i), 1);
	}
	
	for (int i = 0; i < M; i++) {
		cc.add_edge(bull(i), 1, 1);
	}
	
	for (int i = 0; i < P; i++) {
		scanf("%d%d", &A, &B);
		A--; B--;
		cc.add_edge(cow(A), bull(B), 1);
	}
	
	cout << cc.max_flow(0, 1) << endl;
}
