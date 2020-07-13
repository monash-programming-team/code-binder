// Bellman-Ford Algorithm
//
// Author      : Daniel Anderson
// Date        : April 9, 2016
// Reliability : 5
// Tested On   : Codeforces 20C, UVA558, UVA10449, UVA11721, UVA515
//
// Computes the shortest distance between a source vertex and all
// other vertices in a weighted graph (directed or undirected) that may
// contain negative edge weights. Can also be used to detect and find
// negative cycles.
//
// Complexity: O(VE)
#include <bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef pair<int,int> pii;
typedef vector<vi> vvi;

//listings:bellman_ford
// Shortest paths and negative cycle finding in graphs with any weights.
// dist[u] = INF if u is not reachable. dist[u] = -INF if u is reachable via a
// negative cycle. T is the type of the edge weights / costs. Complexity: O(VE)
template <typename T> struct BellmanFord {
  typedef pair<T, int> pti;  vector<vector<pti>> adj;
  int n, last = -1; const T INF = numeric_limits<T>::max() / 2;
  BellmanFord(int n) : adj(n), n(n) {}
  void add_edge(int u, int v, T weight) { adj[u].emplace_back(weight, v); }
  pair<vector<T>, vi> shortest_paths(int src) {
    vector<T> dist(n, INF);  dist[src] = 0;  vi pred(n, -1);  last = 0;
    for (int k = 0; k < n && last != -1; k++) {  last = -1;
      for (int u = 0; u < n; u++) if (dist[u] < INF) for (auto &e : adj[u]) {
        int v = e.Y; T len = dist[u] + e.X;
        if (len < dist[v]) dist[v] = len, pred[v] = u, last = v;
      }
    }
    if (last == -1) return {dist, pred};  // there were no negative cycles
    for (int k = 0, upd = 1; k < n && upd; k++) {  upd = 0;
      for (int u = 0; u < n; u++) if (dist[u] < INF) for (auto &e : adj[u]) {
        int v = e.Y; T len = dist[u] + e.X;
        if (len < dist[v]) dist[v] = -INF, upd = 1;
      }
    }
    return {dist, pred};  // there was a negative cycle
  }  // Returns true if the most recent invocation encountered a negative cycle
  bool had_negative_cycle() { return last != -1; }
  // OPTIONAL: Find a negative cycle in the graph
  vi find_negative_cycle() {
    n++; adj.resize(n);  // add a new temp vertex
    for (int v = 0; v < n - 1; v++) add_edge(n-1, v, 0);
    vi C, pred = shortest_paths(n-1).Y;
    n--; adj.resize(n);  // delete the temp vertex
    if (!had_negative_cycle()) return C; // no negative cycle found
    for (int i = 0; i < n; i++) last = pred[last];
    for (int u = last; u != last || C.empty(); u = pred[u]) C.push_back(u);
    reverse(C.begin(), C.end());
    return C;
  }
};
//listings:/bellman_ford

vi get_path(int v, vi &pred) {
	vi p = {v};
	while (pred[v] != -1) p.push_back(v = pred[v]);
	reverse(p.begin(), p.end());
	return p;
}

// ---------------------------------------------------
// Bellman-Ford ends. Examples
// --------------------------------------------------

// SCC. Needed for UVA11721
struct SCC {
	int n, comp;  vvi g, gt;  vi seq, vis;
	void dfs(int u, const vvi &adj) {
		for (int v : adj[u]) if (vis[v] == -1) { vis[v] = comp; dfs(v, adj); }
		seq.push_back(u);
	}

	SCC(int n) : n(n), g(n), gt(n) { }
	void add_edge(int u, int v) { g[u].push_back(v); gt[v].push_back(u); }
	pair<int, vi> find_SCC() {
		vis.assign(n, -1); comp = 0;
		for (int i = 0; i < n; i++) if (vis[i] == -1) { vis[i] = comp; dfs(i, g); }
		vis.assign(n, -1); comp = 0;
		for (int i = n-1; i >= 0; i--) {
			int u = seq[i];
			if (vis[u] == -1) { vis[u] = comp; dfs(u, gt); comp++; }
		}
		return {comp, vis};
	}
	vvi get_dag() { // OPTIONAL: find_SCC() must be called first
		map<pii, int> mmap;  vvi dag(comp, vi());
		for (int u = 0; u < n; u++) for (int v : g[u]) {
			if (vis[u] == vis[v]) continue;
			if (!mmap.count(pii(vis[u], vis[v]))){
				dag[vis[u]].push_back(vis[v]);
				mmap[pii(vis[u], vis[v])] = 1;
			}
		}
		return dag;
	}
};

map<pii, int> edges;

// Verdict: AC
void solve_uva558() {
  int c;  cin >> c;
  while (c--) {
    int n, m;
    cin >> n >> m;
    BellmanFord<int> G(n + 1);
    for (int i = 0; i < m; i++) {
      int x, y, t;
      cin >> x >> y >> t;
      G.add_edge(x, y, t);
      edges[{x,y}] = t;
    }
    for (int i = 0; i < n; i++) G.add_edge(n, i, 0);
    auto d = G.shortest_paths(n);
    if (G.had_negative_cycle()) {
      cout << "possible" << endl;
      // Verify cycle finding
      auto c = G.find_negative_cycle();
      assert(!c.empty());
      int w = edges[{c.back(), c.front()}];
      for (int i = 0; i < (int)c.size() - 1; i++) w += edges[{c[i],c[i+1]}];
      assert(w < 0);
    }
    else {
      cout << "not possible" << endl;
      // Verify cycle finding
      auto c = G.find_negative_cycle();
      assert(c.empty());
    }
  }
}

// Verdict: AC
int dfs(int c, vvi& dag, vi& comp_cycles, vi& vis) {
  if (vis[c]) return comp_cycles[c];
  vis[c] = 1;
  for (auto v : dag[c])
    if (dfs(v, dag, comp_cycles, vis))
      return (comp_cycles[c] = 1);
  return 0;
}
void solve_uva11721() {
  int cn = 1;
  int T; cin >> T;
  while (T--) {
    int n, m;
    cin >> n >> m; n++;
    SCC scc(n);
    vector<vector<pii>> adj(n);
    for (int i = 0; i < m; i++) {
      int x, y, t; cin >> x >> y >> t;
      scc.add_edge(x, y);
      adj[x].emplace_back(y, t);
    }
    
    // Find SCCs
    int num_scc; vi comp;
    tie(num_scc, comp) = scc.find_SCC();
    vvi comps(num_scc); vi comp_cycles(num_scc), vis(num_scc), ans;
    for (int v = 0; v < n; v++) comps[comp[v]].push_back(v);
    
    // Find which SCCs contain negative cycles
    for (int c = 0; c < num_scc; c++) {
      sort(comps[c].begin(), comps[c].end());
      map<int,int> id;
      for (int i = 0; i < (int)comps[c].size(); i++)
        id[comps[c][i]] = i;
      BellmanFord<int> G((int)comps[c].size());
      for (int u : comps[c]) for (auto v : adj[u])
        if (binary_search(comps[c].begin(), comps[c].end(), v.X))
          G.add_edge(id[u], id[v.X], v.Y);
      auto dist = G.shortest_paths(0).X;
      if (dist[0] == -G.INF) comp_cycles[c] = 1, vis[c] = 1;
    }
    
    // Find which SCCs can reach the negative cycles
    auto dag = scc.get_dag();
    for (int c = 0; c < num_scc; c++) if (dfs(c, dag, comp_cycles, vis)) 
        for (auto& x : comps[c]) ans.push_back(x);
     
    // Output the answer
    sort(ans.begin(), ans.end());
    cout << "Case " << cn++ << ":";
    if (ans.empty()) cout << " impossible" << '\n';
    else { for (auto x : ans) cout << ' ' << x;  cout << '\n'; }
  }
}

// Verdict: AC
template<typename T> T cube(T x) { return x*x*x; }
void solve_uva10449() {
  int n, cn = 1;
  while (cin >> n) {
    cout << "Set #" << cn++ << '\n';
    if (n == 0) { cin >> n >> n; continue; }
    vi b(n);
    for (int i = 0; i < n; i++) cin >> b[i];
    BellmanFord<ll> G(n);
    int r; cin >> r;
    for (int i = 0; i < r; i++) {
      int u, v; cin >> u >> v; u--; v--;
      G.add_edge(u, v, cube(b[v] - b[u]));
    }
    int q; cin >> q;
    auto res = G.shortest_paths(0);
    for (int i = 0; i < q; i++) {
      int x; cin >> x; x--;
      if (res.X[x] < 3 || res.X[x] == G.INF) cout << '?' << '\n';
      else cout << res.X[x] << '\n';
    }
  }
}

// Verdict: TLE on test 28
void solve_CF20C() {
  int n, m;
	cin >> n >> m;
	BellmanFord<int> g(n);
	
	for (int i=0; i<m; i++) {
		int a, b, w;	cin >> a >> b >> w;	a--; b--;
		g.add_edge(a, b, w);
		g.add_edge(b, a, w);
	}
	
	auto ans = g.shortest_paths(0);
	if (ans.first[n-1] == g.INF) cout << -1 << endl;
  else {
    vi path = get_path(n-1, ans.second);
    for (int p : path) cout << p + 1 << ' '; cout << endl;
  }
}

// Verdict: AC
void solve_uva515() {
  int n, m;
  while (cin >> n, n != 0) {
    cin >> m;
    BellmanFord<int> G(2*n+1);
    while (m--) {
      int ss, nn, kk; string oo;
      cin >> ss >> nn >> oo >> kk;
      nn += ss; ss--;
      if (oo == "lt") { kk--; G.add_edge(ss, nn, kk); }
      else { kk++; G.add_edge(nn, ss, -kk); }
    }
    for (int v = 0; v < 2*n; v++) G.add_edge(2*n, v, 0);
    if (G.shortest_paths(2*n), G.had_negative_cycle()) cout << "successful conspiracy" << '\n';
    else cout << "lamentable kingdom" << '\n';
  }
}

int main() {
  solve_uva558();
  //solve_uva10449();
  //solve_uva11721();
  //solve_CF20C();
  //solve_uva515();
}
