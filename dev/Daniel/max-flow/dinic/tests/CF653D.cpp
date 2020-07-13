// Dinic test on Codeforces 653D
//
// Expected verdict: Accepted
//

#include<bits/stdc++.h>
using namespace std;

typedef vector<vector<pair<int, int>>> graph;

typedef long double ld;
typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

const ld EPS = 1e-10;
const ll INF = numeric_limits<ll>::max();

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

int n, m, x;
graph adj;

bool check(ld weight) {
  ld wpb = weight / x;  // weight carried per bear
  dinic<ll> mf(n);
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
