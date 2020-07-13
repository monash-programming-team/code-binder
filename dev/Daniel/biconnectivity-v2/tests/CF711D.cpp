// Test for edge biconnectivity -- identifying bridge connected components
//
// Currently WA

#include <bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef long long int ll;

struct IOS { IOS() { ios::sync_with_stdio(0); cin.tie(0); } } IOS;

const ll MOD = 1000000007LL;

ll expmod(ll a,ll b,ll mod = MOD) {ll res=1;a%=mod;for(;b;b>>=1){if(b&1)res=res*a%mod;a=a*a%mod;}return res;}

class edge_biconnectivity {
  struct edge {
    int u, v;
    bool used, bridge = false;
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

void dfs(int u, int& sz, vvi& adj, vi& vis, edge_biconnectivity& b) {
  sz++;
  vis[u] = true;
  for (auto i : adj[u]) {
    auto& e = b.get_edge(i);
    int v = e.other(u);
    if (!vis[v] && !e.bridge)
      dfs(v, sz, adj, vis, b);
  }
}

int main() {

  int n;
  cin >> n;
  
  vvi adj;
  adj.resize(n);
  
  edge_biconnectivity bcc(n);
  
  int a;
  for (int i = 0; i < n; i++) {
    cin >> a; a--;
    int e = bcc.add_edge(a, i);
    adj[a].push_back(e);
    adj[i].push_back(e);
  }
  
  int num_comps;
  vi comp;

  tie(num_comps, comp) = bcc.find_components();
  
  vi sizes(num_comps);
  for (int i = 0; i < n; i++) {
    assert(comp[i] >= 0 && comp[i] < num_comps);
    sizes[comp[i]]++;
  }
  
  int num_bridges = 0;
  for (int i = 0; i < n; i++)
    if (bcc.get_edge(i).bridge)
      num_bridges++;
  
  cerr << "num_bridges = " << num_bridges << endl;
  
  vi cc_sizes;
  vi vis(n);
  for (int i = 0; i < n; i++) {
    if (!vis[i]) {
      int sz = 0;
      dfs(i, sz, adj, vis, bcc);
      cc_sizes.push_back(sz);
    }
  }
  
  sort(sizes.begin(), sizes.end());
  sort(cc_sizes.begin(), cc_sizes.end());
  
  assert(sizes.size() == cc_sizes.size());
  assert(sizes.size() == num_comps);
  
  assert(accumulate(sizes.begin(), sizes.end(), 0) == accumulate(cc_sizes.begin(), cc_sizes.end(), 0));
  
  for (int i = 0; i < num_comps; i++) {
    if (sizes[i] != cc_sizes[i]) {
      cout << i << ": " << sizes[i] << " != " << cc_sizes[i] << endl;  
    }
  }
  
  ll ans = expmod(2, num_bridges);
  for (auto sz : sizes) {
    if (sz != 1)
      ans = (ans * (expmod(2, sz) - 2)) % MOD;
  }
  cout << ans << endl;

  return 0;
}
