// Edmonds' blossom matching algorithm for general graph matchings.
//
// Author      : Daniel Anderson (based on indy256 and Darcy binder)
// Date        : 03-12-2016
// Reliability : 1
// Tested On   : COJ1192
//
// Finds maximum matchings in general, unweighted graphs.
// 
// Complexity: O(V^3)

#include<bits/stdc++.h>
using namespace std;

#include "../code-template/debug.h"

#define X first
#define Y second

typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:matching
// Maximum matching in a general, unweighted graph. Returns the number of matches
// and a vector containing each node's match, or -1 if no match.  Complexity: O(V^3)
struct GraphMatching {
  int n, m;  vi match, p, base;  vvi adj;
  int lca(int a, int b) {
    vi used(n);
    while (1) { a=base[a], used[a]=1; if (match[a] == -1) break; a = p[match[a]]; }
    while (1) { b = base[b]; if (used[b]) return b; b = p[match[b]]; }
  }
  void mark_path(vi& blossom, int v, int b, int c) {
    for (; base[v] != b; v = p[match[v]])
      blossom[base[v]] = blossom[base[match[v]]] = 1,  p[v] = c,  c = match[v];
  }
  int find_path(int root) {
    vi used(n); iota(base.begin(), base.end(), 0);  p.assign(n, -1);
    used[root] = 1;  queue<int> q;  q.push(root);
    while (!q.empty()) {
      int v = q.front(); q.pop();
      for (int u : adj[v]) {
        if (base[v] == base[u] || match[v] == u) continue;
        if (u == root || (match[u] != -1 && p[match[u]] != -1)) {
          int cb = lca(v, u);  vi blossom(n);
          mark_path(blossom, u, cb, v), mark_path(blossom, v, cb, u);
          for (int i=0; i<n; i++) if (blossom[base[i]]) {
            base[i] = cb; if (!used[i]) used[i] = 1, q.push(i);
          }
        } else if (p[u] == -1) {
          p[u] = v; if (match[u] == -1) return u;
          u = match[u], used[u] = 1, q.push(u);
        }
      }
    }
    return -1;
  }
  // Create a graph on n vertices
  GraphMatching(int n) : n(n), m(0), base(n), adj(n) { }
  void add_edge(int u, int v) { adj[u].push_back(v); adj[v].push_back(u); }
  pair<int,vi> max_matching() {  // Returns the number of matches and each node's match
    p.assign(n, -1), match.assign(n, -1);
    for (int i=0; i<n; i++) {
      if (match[i] != -1) continue;
      int v = find_path(i), ppv = -1;
      while (v != -1) ppv = match[p[v]], match[v] = p[v], match[p[v]] = v, v = ppv;
    }
    return {(n - count(match.begin(), match.end(), -1)) / 2, match};
  }
};
//listings:/matching

// ----------------------------------------------------------------------------
//                            TEST PROBLEMS
// ----------------------------------------------------------------------------

void solve_COJ1192() {
  int N; cin >> N;
  GraphMatching G(N);
  int u,v ;
  while (cin >> u >> v) G.add_edge(u-1,v-1);
  cout << G.max_matching().X * 2 << endl;
}

int main() {
  solve_COJ1192();
}