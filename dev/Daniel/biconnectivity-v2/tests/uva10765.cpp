// Test for vertex biconnectivity -- identifying biconnected components
//
// Verdict: AC
//
#include <bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;

class vertex_biconnectivity {
  struct edge {
    int u, v, comp;
    bool used;
    edge(int a, int b) : u(a), v(b) { }
    int other(int w) { return w == u ? v : u; }
  };
  int n; vvi comps;
  vector<edge> edges; vvi adj;

  int dfs_root, dfs_count, root_children;
  vi dfs_num, dfs_low, cut_point, cur;

  void dfs(int u) {
    dfs_low[u] = dfs_num[u] = dfs_count++;
    for (auto i : adj[u]) {
      auto& e = edges[i]; int v = e.other(u);
      if (e.used) continue; e.used = true;
      if (dfs_num[v] == -1) {
        if (u == dfs_root) root_children++;
        cur.push_back(i);
        dfs(v);
        if (dfs_low[v] >= dfs_num[u]) {
          cut_point[u] = true;
          comps.emplace_back(cur.rbegin(), find(cur.rbegin(), cur.rend(), i) + 1);
          cur.resize(cur.size() - comps.back().size());
          for (auto j : comps.back()) edges[j].comp = (int)comps.size() - 1;
        }
        dfs_low[u] = min(dfs_low[u], dfs_low[v]);
      } else {
        dfs_low[u] = min(dfs_low[u], dfs_num[v]);
        if (dfs_num[v] < dfs_num[u]) cur.push_back(i);
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
  vertex_biconnectivity(int n) : n(n), adj(n) { }

  // OPTIONAL: Construct the block-cut tree of the graph
  // Must run find_bccs() first
  // WARNING: Complexity is O(n^2 + m)
  vvi get_tree() {
    vvi mat(n, vi(comps.size()));
    for (int c = 0; c < (int)comps.size(); c++) {
      for (auto i : comps[c]) {
        auto& e = edges[i];
        mat[e.u][c] = 1;
        mat[e.v][c] = 1;
      }
    }
    vvi tree(comps.size());
    for (int u = 0; u < n; u++) {
      vi my_comps;
      for (int c = 0; c < (int)comps.size(); c++) {
        if (mat[u][c]) my_comps.push_back(c);
      }
      for (auto c1 : my_comps) {
        for (auto c2 : my_comps) {
          if (c1 != c2) tree[c1].push_back(c2);
        }
      }
    }
    return tree;
  }

  // Return a reference to the i'th edge. Use to check
  // if a particular edge was marked as a bridge
  edge& get_edge(int i) { return edges[i]; }

  // Search the graph and return the articulation points and 
  // biconnected components.
  pair<vi, vvi> find_components() {
    dfs_num.assign(n, -1); dfs_low.assign(n, 0); dfs_count = 0;
    cur.clear(); comps.clear(); cut_point.assign(n, 0);
    for (auto& e : edges) e.used = false;
    for (int v = 0; v < n; v++) {
      if (dfs_num[v] == -1) {
        dfs_root = v; root_children = 0;
        dfs(v);
        cut_point[v] = (root_children > 1);
      }
    }
    return {cut_point, comps};
  }

};


int main() {

  int n, m;
  while (cin >> n >> m, n > 0) {
    vertex_biconnectivity b(n);
    int u, v;
    while (cin >> u >> v, u != -1) {
      b.add_edge(u, v);
    }
  
    auto comps = b.find_components();

    vi num_comps(n);
    for (auto comp : comps.second) {
      set<int> contains;
      for (int w : comp) {
        u = b.get_edge(w).u;
        v = b.get_edge(w).v;
        contains.insert(u);
        contains.insert(v);
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
