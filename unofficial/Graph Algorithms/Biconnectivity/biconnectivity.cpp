// Biconnectivity Algorithms
//
// Finds articulation points, brdiges, bridge-connected components
// and biconnected components.
//
// Author      : Daniel Anderson
// Date        : 5/12/2016
// Reliability : 5
// Tested on   : CF711D - Identifying bridge-connected components
//               SPOJ: Submerge - Counting articulation points
//               UVA315 - Counting articulation points
//               UVA796 - Identifying bridges
//               UVA10199 - Identifying articulation points
//               UVA10765 - Identifying biconnected components
//
// Usage:
//  void find_components()
//    Computes all of the articulation points, brdiges, bridge-connected components
//    and biconnected components. See comment above class for details.
//    Complexity O(V + E)  
//
//  edge& get_edge(int i)
//    Returns a reference to the i'th edge. Useful to check the component
//    number of a specific edge.
//    Complexity O(1)
//
#include <bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:biconnectivity
// Find articulation points, bridges, biconnected components and bridge-connected
// components. cut_point[v] = true if v is an articulation point. e.bridge = true
// if e is a bridge. n_vcomps is the number of biconnected components, n_bcomps
// is the number of bridge-connected components. bccs contains biconnected
// components specified by edge indices. bcomp[v] is the index of the
// bridge-connected component containing vertex v.  Complexity: O(V + E)
struct Biconnectivity {
  struct edge {
    int u, v, vcomp;  bool used, bridge;
    edge(int a, int b) : u(a), v(b) { }
    int other(int w) { return w == u ? v : u; }
  };
  int n, m, n_bcomps, n_vcomps, dfs_root, dfs_count, root_children;
  vi dfs_num, dfs_low, cut_point, vcur, bcur, bcomp; vvi bccs, adj; vector<edge> edges;
  void make_vcomp(int i) {  // omit if biconnected components are not required
    bccs.emplace_back(vcur.rbegin(), find(vcur.rbegin(), vcur.rend(), i) + 1);
    vcur.resize(vcur.size() - bccs.back().size());
    for (auto j : bccs.back()) edges[j].vcomp = n_vcomps; n_vcomps++;
  }
  void make_bcomp(int v) {  // omit if bridge-connected components are not required
    int u = -1; n_bcomps++;
    while (u != v) { u = bcur.back(); bcur.pop_back(); bcomp[u] = n_bcomps - 1; }
  }
  void dfs(int u) {
    dfs_low[u] = dfs_num[u] = dfs_count++;
    for (auto i : adj[u]) if (!edges[i].used) {
      auto& e = edges[i]; int v = e.other(u);  e.used = true;
      if (dfs_num[v] == -1) {
        if (u == dfs_root) root_children++;
        vcur.push_back(i), bcur.push_back(v), dfs(v);
        if (dfs_low[v] > dfs_num[u]) { e.bridge = true; make_bcomp(v); }
        if (dfs_low[v] >= dfs_num[u]) { cut_point[u] = true; make_vcomp(i); }
        dfs_low[u] = min(dfs_low[u], dfs_low[v]);
      } else {
        dfs_low[u] = min(dfs_low[u], dfs_num[v]);
        if (dfs_num[v] < dfs_num[u]) vcur.push_back(i);
      }
    }
  }
  Biconnectivity(int n) : n(n), m(0), adj(n) { }
  edge& get_edge(int i) { return edges[i]; }
  int add_edge(int u, int v) {
    adj[u].push_back(m), adj[v].push_back(m), edges.emplace_back(u, v);
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
//listings:/biconnectivity

// Graph with 7 vertices
Biconnectivity G(7);

string edge(int i) {
	return "(" + to_string(G.get_edge(i).u) + "," + to_string(G.get_edge(i).v) + ")";
}

int main() {

  G.add_edge(0,1);
  G.add_edge(1,2);
  G.add_edge(2,3);
  G.add_edge(3,4);
  G.add_edge(4,2);
  G.add_edge(2,5);
  G.add_edge(5,0);
  G.add_edge(4, 6);
  
  G.find_components();

  cout << "Number of bicomponents: " << G.n_vcomps << endl;
  for (auto& comp : G.bccs) {
    for (auto i : comp) cout << edge(i) << ' ';
    cout << endl;
  }
  
  cout << "Articulation points: ";
  for (int v = 0; v < 6; v++) if (G.cut_point[v]) cout << v << ' ';
  cout << endl;
  
  cout << "Number of bridge-connected components: " << G.n_bcomps << endl;
  cout << "Bridges: ";
  for (int i = 0; i < G.m; i++) if (G.get_edge(i).bridge) cout << edge(i) << ' ';

  return 0;
}
