// Test for vertex biconnectivity -- identifying biconnected components
//
// Verdict: AC
//
#include <bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;

struct biconnectivity {
  struct edge {
    int u, v, vcomp;  bool used, bridge;
    edge(int a, int b) : u(a), v(b) { }
    int other(int w) { return w == u ? v : u; }
  };
  int n, m, n_bcomps, n_vcomps, dfs_root, dfs_count, root_children;
  vi dfs_num, dfs_low, cut_point, vcur, bcur, bcomp;  vvi bccs, adj;
  vector<edge> edges;
  void make_vcomp(int i) {
    bccs.emplace_back(vcur.rbegin(), find(vcur.rbegin(), vcur.rend(), i) + 1);
    vcur.resize(vcur.size() - bccs.back().size());
    for (auto j : bccs.back()) edges[j].vcomp = n_vcomps; n_vcomps++;
  }
  void make_bcomp(int v) {
    int u = -1; n_bcomps++;
	while (u != v) { u = bcur.back(); bcur.pop_back(); bcomp[u] = n_bcomps - 1; }
  }
  void dfs(int u) {
    dfs_low[u] = dfs_num[u] = dfs_count++;
    for (auto i : adj[u]) {
      auto& e = edges[i]; int v = e.other(u);
      if (e.used) continue; e.used = true;
      if (dfs_num[v] == -1) {
        if (u == dfs_root) root_children++;
        vcur.push_back(i); bcur.push_back(v);
        dfs(v);
		if (dfs_low[v] > dfs_num[u]) { e.bridge = true; make_bcomp(v); }
        if (dfs_low[v] >= dfs_num[u]) { cut_point[u] = true; make_vcomp(i); }
        dfs_low[u] = min(dfs_low[u], dfs_low[v]);
      } else {
        dfs_low[u] = min(dfs_low[u], dfs_num[v]);
        if (dfs_num[v] < dfs_num[u]) vcur.push_back(i);
      }
    }
  }

  biconnectivity(int n) : n(n), m(0), adj(n) { }
  edge& get_edge(int i) { return edges[i]; }
  int add_edge(int u, int v) {
    adj[u].push_back(m);  adj[v].push_back(m);
    edges.emplace_back(u, v);
	return m++;
  }
  void find_components() {
    dfs_num.assign(n, -1); dfs_low.assign(n, 0); dfs_count = 0;
    vcur.clear(); bcur.clear(); bccs.clear(); cut_point.assign(n, 0);
	bcomp.assign(n, -1); n_vcomps = 0, n_bcomps = 0;
    for (auto& e : edges) e.used = false, e.bridge = false;
    for (int v = 0; v < n; v++) if (dfs_num[v] == -1) {
      bcur = {v}; dfs_root = v;  root_children = 0;  dfs(v);
      cut_point[v] = (root_children > 1); make_bcomp(v);
    }
  }
};


int main() {

  int n, m;
  while (cin >> n >> m, n > 0) {
    biconnectivity G(n);
    int u, v;
    while (cin >> u >> v, u != -1) {
      G.add_edge(u, v);
    }
  
    G.find_components();

    vi num_comps(n);
    for (auto comp : G.bccs) {
      set<int> contains;
      for (int w : comp) {
        auto& e = G.get_edge(w);
        contains.insert(e.u);
        contains.insert(e.v);
      }
      for (auto w : contains) {
        num_comps[w]++;
      }
    }
    
    vector<pair<int,int>> ans;
    for (int w = 0; w < n; w++) {
      ans.emplace_back(-num_comps[w], w);
    }
    
    sort(ans.begin(), ans.end());
    
    for (int w = 0; w < m; w++) {
      cout << ans[w].second << ' ' << (-ans[w].first) << endl;
    }
    cout << endl;
  
  }

  return 0;
}
