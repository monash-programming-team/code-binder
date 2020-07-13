// Test for edge biconnectivity -- identifying bridge connected components
//
// Expected Verdict: AC

#include <bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef long long int ll;

struct IOS { IOS() { ios::sync_with_stdio(0); cin.tie(0); } } IOS;

const ll MOD = 1000000007LL;

ll expmod(ll a,ll b,ll mod = MOD) {ll res=1;a%=mod;for(;b;b>>=1){if(b&1)res=res*a%mod;a=a*a%mod;}return res;}

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

void dfs(int u, int& sz, vvi& adj, vi& vis, biconnectivity& G) {
  sz++;
  vis[u] = true;
  for (auto i : adj[u]) {
    auto& e = G.get_edge(i);
    int v = e.other(u);
    if (!vis[v] && !e.bridge)
      dfs(v, sz, adj, vis, G);
  }
}

int main() {

  int n;
  cin >> n;
  
  vvi adj;
  adj.resize(n);
  
  biconnectivity G(n);
  
  int a;
  for (int i = 0; i < n; i++) {
    cin >> a; a--;
    int e = G.add_edge(a, i);
    adj[a].push_back(e);
    adj[i].push_back(e);
  }
  
  G.find_components();
  
  vi sizes(G.n_bcomps);
  for (int i = 0; i < n; i++) {
    assert(G.bcomp[i] >= 0 && G.bcomp[i] < G.n_bcomps);
    sizes[G.bcomp[i]]++;
  }
  
  int num_bridges = 0;
  for (int i = 0; i < n; i++)
    if (G.get_edge(i).bridge)
      num_bridges++;
  
  cerr << "num_bridges = " << num_bridges << endl;
  
  vi cc_sizes;
  vi vis(n);
  for (int i = 0; i < n; i++) {
    if (!vis[i]) {
      int sz = 0;
      dfs(i, sz, adj, vis, G);
      cc_sizes.push_back(sz);
    }
  }
  
  sort(sizes.begin(), sizes.end());
  sort(cc_sizes.begin(), cc_sizes.end());
  
  assert(sizes.size() == cc_sizes.size());
  assert((int)sizes.size() == G.n_bcomps);
  
  assert(accumulate(sizes.begin(), sizes.end(), 0) == accumulate(cc_sizes.begin(), cc_sizes.end(), 0));
  
  for (int i = 0; i < G.n_bcomps; i++) {
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
