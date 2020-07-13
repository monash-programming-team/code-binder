// Test edge biconnectivity - Identifying bridges
//
// Verdict: AC
//
#include <bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;

class edge_biconnectivity {
  struct edge {
    int u, v;
    bool used, bridge;
    edge(int a, int b) : u(a), v(b) { }
    int other(int w) { return w == u ? v : u; }
  };
  int n, n_comps; vi comp;
  vector<edge> edges; vvi adj;

  int dfs_count;
  vi dfs_num, dfs_low, cur;

  void new_comp(int u) {
    bool done = false;    
    while (!done) {
      comp[cur.back()] = n_comps;
      if (cur.back() == u) done = true;
      cur.pop_back();
    }
    n_comps++;
  }

  void dfs(int u) {
    dfs_low[u] = dfs_num[u] = dfs_count++;
    for (auto i : adj[u]) {
      auto& e = edges[i]; int v = e.other(u);
      if (e.used) continue; e.used = true;
      if (dfs_num[v] == -1) {
        cur.push_back(v);
        dfs(v);
        if (dfs_low[v] > dfs_num[u]) { 
          e.bridge = true;
          new_comp(v);
        }
        dfs_low[u] = min(dfs_low[u], dfs_low[v]);
      } else {
        dfs_low[u] = min(dfs_low[u], dfs_num[v]);
      }
    }
  }

 public:
  // Add an undirected edge from vertex u to vertex v
  // Returns the index of the edge
  int add_edge(int u, int v) {
    adj[u].push_back((int)edges.size());
    adj[v].push_back((int)edges.size());
    edges.emplace_back(u, v);
    return (int)edges.size() - 1;
  }

  // Create an undirected graph with n vertices
  edge_biconnectivity(int n) : n(n), adj(n) { }

  // Construct the tree of bridge-connected components
  // Must run find_bridges() first
  vvi get_tree() {
    vvi tree(n_comps);
    for (auto& e : edges) {
      if (e.bridge) {
        tree[comp[e.u]].push_back(comp[e.v]);
        tree[comp[e.v]].push_back(comp[e.u]);
      }
    }
    return tree;
  }

  // Returns the number of edges in the graph
  int num_edges() { return (int)edges.size(); }

  // Return a reference to the i'th edge. Use to check
  // if a particular edge was marked as a bridge
  edge& get_edge(int i) { return edges[i]; }

  // Search the graph and return the bridge-connected components
  pair<int, vi> find_components() {
    dfs_num.assign(n, -1); dfs_low.assign(n, 0); dfs_count = 0;
    for (auto& e : edges) e.used = false, e.bridge = false;
    comp.assign(n, -1); n_comps = 0;
    for (int v = 0; v < n; v++) {
      if (dfs_num[v] == -1) {
        cur = {v};
        dfs(v);
        new_comp(v);
      }
    }
    return {n_comps, comp};
  }

};


int main() {

  int n;
  while (cin >> n) {
    edge_biconnectivity bcc(n);
    int num_edges = -1;
    int u, v; string num;
    for (int i = 0; i < n; i++) {
      cin >> u >> num;
      num = num.substr(1, num.size() - 2);
      int con = stoi(num);
      for (int j = 0; j < con; j++) {
        cin >> v;
        if (u < v) num_edges = max(num_edges, bcc.add_edge(u, v));
      }
    }
    
    auto res = bcc.find_components();
    
    vector<pair<int,int>> ans;
    for (int i = 0; i <= num_edges; i++) {
      auto& e = bcc.get_edge(i);
      if (e.bridge) {
        ans.push_back(minmax(e.u, e.v));
      }
    }
    
    sort(ans.begin(), ans.end());
    
    cout << ans.size() << " critical links" << endl;
    for (auto a : ans) {
      cout << a.first << " - " << a.second << endl;
    }
    cout << endl;  
  
  }

  return 0;
}
