// Biconnectivity Algorithms for Vertex-Biconnectivity
//
// Finds articulation points and biconnected components.
//
// Author      : Daniel Anderson
// Date        : 15/09/2016
// Reliability : 4
// Tested on   : UVA315, UVA10199, SPOJ-SUBMERGE, UVA10765
//
// Usage:
//  pair<vi, vvi> find_components()
//    Returns two vectors. The first is of length n and contains true if
//    vertex i is an articulation point, false otherwise. The second contains
//    the biconnected components, represented as vectors of edge ids.
//    Complexity O(n + m)  
//
//  edge& get_edge(int i)
//    Returns a reference to the i'th edge. Useful to check the component
//    number of a specific edge.
//    Complexity O(1)
//
//  vvi get_block_graph() [Not tested!]
//    Find the Block-Graph of the graph. This is a graph where each
//    vertex represents a biconnected component, and contains an edge
//    to any other biconnected component that shares an articulation point.
//    Complexity O(n^2 + m)
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

  // OPTIONAL: Construct the block-graph of the graph
  // Must run find_bccs() first
  // WARNING: Not tested thoroughly
  // WARNING: Complexity is O(n^2 + m)
  vvi get_block_graph() {
    vvi mat(n, vi(comps.size()));
    for (int c = 0; c < (int)comps.size(); c++) {
      for (auto i : comps[c]) {
        auto& e = edges[i];
        mat[e.u][c] = 1;
        mat[e.v][c] = 1;
      }
    }
    vvi G(comps.size());
    for (int u = 0; u < n; u++) {
      vi my_comps;
      for (int c = 0; c < (int)comps.size(); c++) {
        if (mat[u][c]) my_comps.push_back(c);
      }
      for (auto c1 : my_comps) {
        for (auto c2 : my_comps) {
          if (c1 != c2) G[c1].push_back(c2);
        }
      }
    }
    return G;
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

  vertex_biconnectivity b(6);
  b.add_edge(0,1);
  b.add_edge(1,2);
  b.add_edge(2,3);
  b.add_edge(3,4);
  b.add_edge(4,2);
  b.add_edge(2,5);
  b.add_edge(5,0);
  
  auto res = b.find_components();

  cout << "Number of components: " << res.second.size() << endl;
  for (auto& comp : res.second) {
    for (auto i : comp) cout << i << ' ';
    cout << endl;
  }
  
  cout << "Articulation points: ";
  for (auto b : res.first) cout << b << ' ';
  cout << endl;
  
  auto G = b.get_block_graph();
  
  cout << "BLOCK GRAPH: " << endl;
  for (int i = 0; i < (int)G.size(); i++) {
    cout << i << ": ";
    for (auto j : G[i]) cout << j << ' ';
    cout << endl;
  }

  return 0;
}
