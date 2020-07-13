// Centroid Decomposition (AKA Divide-and-conquer on trees)
//
// Author: Daniel Anderson
// Date: 11-01-2017
// Reliability: 5
// Tested on: CF342E, HackerRank - BST maintenance, CF321C, CF100570F (Gym), SPOJ - QTREE5
//
#include<bits/stdc++.h>
using namespace std;

#include "../code-template/debug.h"

#define X first
#define Y second

#define all(x) begin(x), end(x)

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> pii;

//listings:centroid_decomposition
// Centroid Decomposition. Constructs a valid centroid tree of the given tree.
// croot -- the root of the centroid tree
// cadj -- downward adjacency list of the centroid tree
// par -- parent in the centroid tree (-1 for the root)
struct CentroidDecomposition {
  int n, cnt = 0, croot;  vvi adj, cadj;  vi par, mark, size;
  int dfs(int u, int p) {
    size[u] = 1;
    for (int v : adj[u]) if (v != p && !mark[v]) dfs(v, u), size[u] += size[v];
    return size[u];
  }
  int find_centroid(int u, int p, int sz) {
    for (int v : adj[u]) if (v != p && !mark[v])
      if (size[v] * 2 > sz) return find_centroid(v, u, sz);
    return u;
  }
  int find_centroid(int src) { return find_centroid(src, -1, dfs(src, -1)); }
  // Create a tree on n vertices -- add edges using add_edge(u, v)
  CentroidDecomposition(int n) : n(n), adj(n), cadj(n), par(n,-1), mark(n), size(n) { }
  void add_edge(int u, int v) { adj[u].push_back(v), adj[v].push_back(u); }
  int decompose_tree(int src = 0) {
    int c = find_centroid(src);  mark[c] = 1;
    for (int u : adj[c]) if (!mark[u]) {
      int v = decompose_tree(u);
      cadj[c].push_back(v), par[v] = c;
    }
    return croot = c;
  }
};
//listings:/centroid_decomposition

// Lowest common ancestor for distance queries
struct LCA {
	const int LOGN = 20;	// works for n <= 10^6. Change appropriately.
	const vvi& adj;  int n;  vi par, lvl;  vvi anc;
	void dfs(int u, int p, int d) {
		par[u] = p; lvl[u] = d;
		for (int v : adj[u]) if (v != p) dfs(v, u, d+1);
	}
	// Create the LCA structure from the given adjacency list at the given root
	LCA(const vvi& adj, int root = 0) : adj(adj), n(adj.size()), par(n), lvl(n) {
		dfs(root,-1,0),  anc.assign(n, vi(LOGN, -1));
		for (int i = 0; i < n; i++) anc[i][0] = par[i];
		for (int k = 1; k < LOGN; k++) for (int i = 0; i < n; i++)
			if (anc[i][k-1] != -1) anc[i][k] = anc[anc[i][k-1]][k-1];
	}
	int query(int u, int v) {  // LCA with respect to original root
		if (lvl[u] > lvl[v]) swap(u, v);
		for (int k = LOGN - 1; k >= 0; k--) 
			if (lvl[v] - (1 << k) >= lvl[u]) v = anc[v][k];
		if (u == v) return u;
		for (int k = LOGN - 1; k >= 0; k--) {
			if (anc[u][k] == anc[v][k]) continue;
			u = anc[u][k]; v = anc[v][k];
		}
		return par[u];
	}
  int query(int u, int v, int root) {  // OPTIONAL: LCA with respect to any root
    int a = query(u, v), b = query(u, root), c = query(v, root);
    if (a == c && c != b) return b;
    else if (a == b && c != b) return c;
    else return a;
  }
	int dist(int u, int v) { return lvl[u] + lvl[v] - 2 * lvl[query(u,v)]; }
};

// ------------------------------------------------------------------
//                      CODEFORCES 342E. Verdict: AC
// ------------------------------------------------------------------
const int INF = INT_MAX / 2;
vi dist;
void paint(int u, LCA& lca, CentroidDecomposition& cd) {
  int x = u;
  while (x != -1) {
    dist[x] = min(dist[x], lca.dist(x,u));
    x = cd.par[x];
  }
}
int get_dist(int u, LCA& lca, CentroidDecomposition& cd) {
  int ans = INT_MAX, x = u;
  while (x != -1) {
    ans = min(ans, lca.dist(x,u) + dist[x]);
    x = cd.par[x];
  }
  return ans;
}  
void solve_CF342E() {
  ios_base::sync_with_stdio(0); cin.tie(0);
  int n, m; cin >> n >> m;
  CentroidDecomposition cd(n);
  for (int i=0; i<n-1; i++) {
    int a, b; cin >> a >> b; a--; b--;
    cd.add_edge(a,b);
  }
  LCA lca(cd.adj);  cd.decompose_tree();  dist.assign(n, INF);  paint(0, lca, cd);
  while (m--) {
    int t, v;
    cin >> t >> v; v--;
    if (t == 1) paint(v, lca, cd);
    else cout << get_dist(v, lca, cd) << '\n';
  }
}

// ------------------------------------------------------------------
//              HackerRank BST-MAINTENANCE. Verdict: AC
// ------------------------------------------------------------------
void solve_hackerrank_bst_maintenance() {
  int N; cin >> N;
  CentroidDecomposition cd(N);
  vi a(N);  for (auto& x : a) { cin >> x; x--; }
  // Build the BST
  set<int> numbers = {a[0]};
  map<int,int> left, right;
  vi parent(N, -1);
  for (int i=1; i<N; i++) {
    auto it = numbers.upper_bound(a[i]);  int result;
    if (it != numbers.end() && left.count(*it) == 0) 
      left[*it] = a[i], result = *it;
    else it--, right[*it] = a[i], result = *it;
    numbers.insert(a[i]), parent[a[i]] = result;
    cd.add_edge(a[i], result);
  }
  // Do Centroid Decomposition
  cd.decompose_tree();  ll ans = 0;
  vector<ll> total_length(N), num_paths(N); vector<bool> active(N);
  vector<unordered_map<int, ll>> length_in_subtree(N), num_paths_in_subtree(N);
  LCA lca(cd.adj);
  // Answer queries
  for (int u : a) {
    active[u] = true;
    ans += total_length[u];
    int x = cd.par[u], pred = u;
    while (x != -1) {
      int dist = lca.dist(x,u);
      // Add to answer -------------------------
      if (active[x])
        ans += (total_length[x] - length_in_subtree[x][pred])
          + (num_paths[x] - num_paths_in_subtree[x][pred] + 1) * dist;
      // Update values -------------------------
      total_length[x] += dist;
      length_in_subtree[x][pred] += dist;
      num_paths[x]++;
      num_paths_in_subtree[x][pred]++;
      // Move up the tree ---------------------
      pred = x, x = cd.par[x];
    }
    cout << ans << '\n';
  }
}

// ------------------------------------------------------------------
//                      CODEFORCES 321C. Verdict: AC
// ------------------------------------------------------------------
void dfs(int u, char rank, vector<char>& ans, const vvi& cadj) {
  ans[u] = rank;
  for (int v : cadj[u]) dfs(v, rank + 1, ans, cadj);
}
void solve_CF321C() {
  int n; cin >> n;
  CentroidDecomposition cd(n);
  for (int i=0; i<n-1; i++) {
    int a, b; cin >> a >> b; a--; b--;
    cd.add_edge(a, b);
  }
  int root = cd.decompose_tree();
  vector<char> assignment(n);
  dfs(root, 'A', assignment, cd.cadj);
  for (int i=0; i<n; i++) cout << assignment[i] << " \n"[i == n-1];
}


// ------------------------------------------------------------------
//             CODEFORCES 100570F (Gym). Verdict: AC
// ------------------------------------------------------------------
void dfs2(int u, int p, ll d, vector<ll>& dist, const vvi& adj, map<pii, ll>& weight) {
  dist[u] = d;
  for (int v : adj[u]) if (v != p)
    dfs2(v, u, d + weight[minmax(u,v)], dist, adj, weight);
}
void solve_CF100570F() {
  int n, q; cin >> n >> q;
  // Build tree and centroid decomposition
  CentroidDecomposition cd(n);
  map<pii, ll> weight;
  for (int i=0; i<n-1; i++) {
    int a, b; ll w; cin >> a >> b >> w;
    a--; b--;
    cd.add_edge(a, b);
    weight[minmax(a,b)] = w;
  }
  cd.decompose_tree();
  // Pre-compute tree distances
  LCA lca(cd.adj); vector<ll> d(n); dfs2(0, -1, 0, d, cd.adj, weight);
  auto get_dist = [&](int u, int v) { return d[u] + d[v] - 2 * d[lca.query(u,v)]; };
  // Pre-compute path lengths
  vector<vector<ll>> to_me(n), to_par(n);
  for (int u = 0; u < n; u++) {
    to_me[u].push_back(0);
    for (int p = cd.par[u], c = u; p != -1; c = p, p = cd.par[p]) {
      ll dist = get_dist(u, p);
      to_me[p].push_back(dist);
      to_par[c].push_back(dist);
    }
  }
  // Sort path lengths for easy counting
  for (int u = 0; u < n; u++) {
    sort(all(to_me[u]));
    sort(all(to_par[u]));
  }
  // Answer queries
  while (q--) {
    int v; ll l; cin >> v >> l; v--;
    int ans = upper_bound(all(to_me[v]), l) - to_me[v].begin();
    for (int p = cd.par[v], c = v; p != -1; c = p, p = cd.par[p]) {
      ll dist = get_dist(v, p);
      int add = upper_bound(all(to_me[p]), l - dist) - to_me[p].begin();
      int overlap = upper_bound(all(to_par[c]), l - dist) - to_par[c].begin();
      ans += add - overlap;
    }
    cout << ans << '\n';
  }
}

// ------------------------------------------------------------------
//                  SPOJ - QTREE5. Verdict: AC
// ------------------------------------------------------------------
void solve_SPOJ_QTREE5() {
  int N; cin >> N;
  CentroidDecomposition cd(N);
  for (int i=0; i<N-1; i++) {
    int a, b; cin >> a >> b; a--, b--;
    cd.add_edge(a, b);
  }
  cd.decompose_tree();  LCA lca(cd.adj);
  vector<multiset<int>> distances(N);  vi colour(N);
  int Q; cin >> Q;
  while (Q--) {
    int t; cin >> t;
    if (t == 0) {
      int i; cin >> i; i--;
      colour[i] ^= 1;
      for (int u = i; u != -1; u = cd.par[u]) {
        if (colour[i] == 1) distances[u].insert(lca.dist(i, u));
        else distances[u].erase(distances[u].find(lca.dist(i, u)));
      }
    } else {
      int v; cin >> v; v--;
      int res = INT_MAX;
      for (int u = v; u != -1; u = cd.par[u])
        if (!distances[u].empty()) res = min(res, *distances[u].begin() + lca.dist(u, v));
      if (res == INT_MAX) cout << -1 << '\n';
      else cout << res << '\n';
    }
  }
}

int main() {
  //solve_CF342E();
  //solve_hackerrank_bst_maintenance();
  //solve_CF321C();
  //solve_CF100570F();
  solve_SPOJ_QTREE5();
}

