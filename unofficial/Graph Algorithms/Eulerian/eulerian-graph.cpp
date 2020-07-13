// Find Euler Path/Tour
//
// Author      : Daniel Anderson
// Date        : 12-12-16
// Reliability : 5
// Tested On   : UVA10054, NEERC_Moscow9I, CF267B, TIMMUS1137, UVA10441
//
// Finds an Eulerian path or tour in a graph if one exists.
//
// Complexity: O( N + M )
#include <bits/stdc++.h>
using namespace std;

//#include "../code-template/debug.h"

#define X first
#define Y second

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> pii;

//listings:eulerian
// Find an Eulerian path or tour in a given graph if one exists. For a connected,
// undirected graph, an Euler tour exists if every vertex has an even degree.
// An Euler path exists if all but two vertices have an even degree, these will
// be the endpoints. A connected directed graph has an Euler tour if all
// vertices have indegree == outdegree, or an Euler path if one vertex has
// outdegree - indegree = 1 and one vertex has indegree - outdegree = 1, these
// will be the start and endpoints respectively. You must check existence yourself.
// Call find(start) where start is the first vertex of the path / tour.
// NOTE: Both the start and end-point are included in a tour. Complexity: O(V + E)
struct Eulerian {
  struct edge { int u, v; bool used; int opp(int x) { return x == u ? v : u; } };
  int n, m; vector<edge> edges;  vvi adj;  vi cnt, tour;
  void dfs(int u) {
    while (cnt[u] < (int)adj[u].size()) {
      auto& e = edges[adj[u][cnt[u]++]];
      if (!e.used) e.used = 1, dfs(e.opp(u)), tour.push_back(u);
    }
    if (tour.empty()) tour.push_back(u);
  }
  Eulerian(int n) : n(n), m(0), adj(n) { }
  void add_edge(int u, int v, bool dir) {  // dir = true if the edge is directed
    edges.push_back({u,v,0}), adj[u].push_back(m++);       // or false otherwise
    if (!dir) adj[v].push_back(m-1);
  }
  vi find(int start=0) {
    tour.clear(); cnt.assign(n, 0);  for (auto& e : edges) e.used = 0;
    dfs(start), reverse(tour.begin(), tour.end());
    return tour;
  }
};
//listings:/eulerian

// ----------------------------------------------------------------------------
//                                   TEST PROBLEMS
// ----------------------------------------------------------------------------

// Undirected Euler Tour
// Verdict: AC
void solve_UVA10054() {
  int T, u, v, cn=1; cin >> T;
  while (T--) {
    if (cn > 1) cout << '\n';
    cout << "Case #" << cn++ << '\n';
    int N; cin >> N;
    Eulerian G(51);
    vi degree(51);
    for (int i=0; i<N; i++) {
      cin >> u >> v;
      G.add_edge(u, v, false);
      degree[u]++; degree[v]++;
    }
    bool good = true;
    for (int c=0; c<=50; c++) if (degree[c] % 2 != 0) good = false;
    auto res = G.find(u);
    if ((int)res.size() != N + 1 || !good) cout << "some beads may be lost" << '\n';
    else for (int i=0; i<N; i++) cout << res[i] << ' ' << res[i+1] << '\n';
  }
}

// Directed Euler Tour
// Verdict: AC if you increase the time limit slightly. Good enough.
void solve_NEERC_moscow9I() {
#ifdef ONLINE_JUDGE
	freopen("infinite.in", "r", stdin);
	freopen("infinite.out", "w", stdout);
#endif
  ios::sync_with_stdio(0); cin.tie(0);
  int N; cin >> N;
  if (N == 1) { cout << "01" << endl; return; }
  int V = 1 << (N - 1), mask = (1 << (N - 2)) - 1;
  Eulerian G(V);
  for (int u=0; u<V; u++) {
    G.add_edge(u, (u & mask) << 1, true);
    G.add_edge(u, ((u & mask) << 1) + 1, true);
  }
  auto res = G.find(0);
  for (int i = 0; i < N - 1; i++) cout << ((res[0] & (1 << i)) > 0);
  for (int i = 1; i < (int)res.size(); i++) cout << (res[i] & 1);
	cout << endl;
}

// Undirected Euler Path or Euler Tour
// Verdict: AC
void solve_CF267B() {
  int n, u, v; cin >> n;
  if (n == 1) { cout << "1 +" << endl; return; }
  Eulerian G(7);
  vi degree(7);
  map<pii, vector<pair<int, pii>>> dominoes;
  for (int i=0; i<n; i++) {
    cin >> u >> v;
    degree[u]++; degree[v]++;
    G.add_edge(u, v, false);
    dominoes[minmax(u,v)].emplace_back(i, pii{u,v});
  }
  // Look for a starting vertex
  vi odds;
  for (int v=0; v<=6; v++) if (degree[v] % 2 == 1) odds.push_back(v);
  if (odds.size() != 0 && odds.size() != 2) {
    cout << "No solution" << endl; return;
  } else {
    auto res = odds.empty() ? G.find(u) : G.find(odds.front());
    if ((int)res.size() != n + 1) { cout << "No solution" << endl; return; }
    for (int i=0; i<n; i++) {
      pii dom = {res[i], res[i+1]}, ref = minmax(res[i], res[i+1]);
      auto use = dominoes[ref].back();  dominoes[ref].pop_back();
      if (dom == use.Y) cout << use.X + 1 << " +" << '\n';
      else cout << use.X + 1 << " -" << '\n';
    }
  }
}

// Directed Euler Tour
// Verdict: AC
const int MAXN = 10001;
void solve_timmus1137() {
  int n, m, tot = 0, u, v; cin >> n;
  // Build the graph
  Eulerian G(MAXN);
  vi indegree(MAXN), outdegree(MAXN);
  for (int i=0; i<n; i++) {
    cin >> m >> u;
    for (int j=0; j<m; j++) {
      cin >> v;
      G.add_edge(u, v, true);
      indegree[v]++; outdegree[u]++;
      u = v;  tot++;
    }
  }
  // Check feasibility
  for (int x=0; x<MAXN; x++) if (indegree[x] != outdegree[x]) { cout << 0 << endl; return; }
  // Find a route
  auto res = G.find(u);
  if ((int)res.size() != tot + 1) { cout << 0 << endl; return; }
  else {
    cout << res.size() - 1;
    for (auto& x : res) cout << ' ' << x;
    cout << endl;
  }
}

// Lexicographically least directed Euler path or tour
// Verdict: AC
void solve_UVA10441() {
  int t; cin >> t;
  while (t--) {
    int n; cin >> n;
    vector<string> dict;
    for (int i=0; i<n; i++) {
      string s; cin >> s;
      dict.push_back(s);
    }
    sort(dict.begin(), dict.end());
    vi indegree(26), outdegree(26);
    vector<vector<string>> edges(26);
    map<pii, vector<string>> words;
    // Build graph
    for (int i=0; i<n; i++) {
      edges[dict[i].front()-'a'].push_back(dict[i]);
      words[pii{dict[i].front()-'a', dict[i].back()-'a'}].push_back(dict[i]);
    }
    // Add edges in lexicographical order
    Eulerian G(26);
    for (int c=0;c<26;c++) {
      sort(edges[c].begin(), edges[c].end());
      for (auto& s : edges[c]) {
        G.add_edge(c, s.back()-'a',true);
        indegree[s.back()-'a']++;
        outdegree[c]++;
      }
    }
    // Sort the edges
    for (auto& E : words) sort(E.Y.rbegin(), E.Y.rend());
    
    // Check feasibility
    vi one_more, one_less;
    bool good = true;
    for (int c=0; c<26; c++) {
      if (indegree[c] != outdegree[c]) {
        if (indegree[c] - outdegree[c] == 1) one_more.push_back(c);
        else if (outdegree[c] - indegree[c] == 1) one_less.push_back(c);
        else good = false;
      }
    }
    bool tour = good && (one_more.empty() && one_less.empty());
    bool path = good && (one_more.size() == 1 && one_less.size() == 1);
    
    if (!tour && !path) { cerr << "Degrees not feasible" << endl; cout << "***" << endl; continue; }
    
    vi res;
    
    // Find a Path
    if (path) {
      cerr << "Seeking path" << endl;
      int start = one_less.back();
      res = G.find(start);
      if ((int)res.size() != n + 1) { cerr << "Path not viable" << endl; cout << "***" << endl; continue; }
    }
    // Find a Tour
    else {
      cerr << "Seeking tour" << endl;
      int start = dict[0].front()-'a';
      res = G.find(start);
      if ((int)res.size() != n + 1) { cerr << "Tour not viable" << endl; cout << "***" << endl; continue; }
    }
    
    // Print the answer
    cerr << "Printing answer" << endl;
    bool first = true;
    for (int v=0; v<n; v++) {
      if (!first) cout << '.';
      first = false;
      auto& es = words[pii{res[v], res[v+1]}];
      cout << es.back(); es.pop_back();
    }
    cout << endl;
  }
}

int main() {
  //solve_UVA10054();
	//solve_NEERC_moscow9I();
  //solve_CF267B();
  //solve_timmus1137();
  solve_UVA10441();
}

