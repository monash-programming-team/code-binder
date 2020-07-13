// Biconnectivity Algorithms for Edge-Biconnectivity
//
// Finds bridges and bridge-connected components.
//
// Author      : Daniel Anderson
// Date        : 15/09/2016
// Reliability : 2
// Tested on   : UVA796, CF711D
//
// Usage:
//  pair<int, vi> find_components()
//    Returns the number of bridge connected components and a vector
//    containing for each vertex the id of its component.
//    Complexity: O(n + m)
//    
//  edge& get_edge(int i)
//    Returns a reference to the i'th edge. Useful to check whether
//    this edge is a bridge or not.
//    Complexity: O(1)
//
//  vvi get_tree()
//    Creates the adjacency list of the bridge-connected component tree.
//    Must call find_bridges() first.
//    Complexity: O(m)
//
#include <bits/stdc++.h>
using namespace std;

#define DEBUG_ON 1
#include "../code-template/debug.h"

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


  return 0;
}
