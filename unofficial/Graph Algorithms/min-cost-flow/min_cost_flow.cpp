// Minimum cost flow using successive shortest paths.
//
// Author      : Daniel Anderson
// Date        : 06-12-2016
// Reliability : 5
// Tested on   : UVA 10594, SPOJ - GREED, UVA 10746, SPOJ - BABY, CF 491C
//
// Usage:
//    ll flow(int s, int t, ll cap)
//   	 Find the minimum cost to send cap units of flow from s to t
// 		 If it is not possible to send cap units of flow from s to t
// 		 then we will find the minimum cost to send the maximum
// 		 amount of flow possible from s to t. If you just want
// 		 min-cost max-flow then use cap = LLONG_MAX
//
// 		Returns: {max_flow, min_cost}
//
//	  int add_edge(u, v, cap, cost)
//			Adds a directed edge from u -> v with capacity and cost.
//      Returns the index of the edge. If you want bidirectional
//      edges, add both (u,v) and (v, u)
//
//	  edge& get_edge(i)
//			Returns a reference to the i'th edge. Use to check flows
//			or update capacities / costs.
//
//	Time Complexity: O(VE + E log(V)*FLOW) theoretically.
//
//  Runs MUCH faster in practice. Has been tested and runs in 0.1 seconds on
//  multi-testcase problems with V = 100, E = 5000 and cost <= 10^15.

#include <bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<ll> vll;

//listings:min_cost_flow
// Minimum cost flow using successive shortest paths. Finds the minimum cost
// to send cap units of flow from s to t. If you want max flow, use cap = INF.
// F = Flow type, C = Cost type.  Complexity: O(VE + E log(V) * FLOW)
template <typename F, typename C> struct MinCostFlow {
  struct edge { int from, to;  F flow, cap;  C cost; };
  const C INF = numeric_limits<C>::max();  vector<C> pi, dist;
  int n, m;  vi pred, pe;  vvi g;  vector<edge> edges;
  typedef pair<C, int> pci;  priority_queue<pci, vector<pci>, greater<pci>> q;
  void bellman_ford() {  // Omit bellman_ford if using Leviticus instead of Dijkstra
    pi.assign(n, 0);
    for (int i = 0; i < n - 1; i++) for (auto &e : edges) if (e.flow < e.cap)
        pi[e.to] = min(pi[e.to], pi[e.from] + e.cost);
  }
  bool dijkstra(int s, int t) {  // Swap this for levit(s, t) for random data
    dist.assign(n, INF); pred.assign(n, -1); dist[s] = 0; q.emplace(0, s);
    while (!q.empty()) {
      C d; int u;  tie(d, u) = q.top(); q.pop();
      if (dist[u] == d) for (int i : g[u]) {
        auto &e = edges[i], v = e.to; C rcost = e.cost + pi[u] - pi[v];
        if (e.flow < e.cap && dist[u] + rcost < dist[v]) 
          pred[v]=u, pe[v]=i, dist[v]=dist[u]+rcost, q.emplace(dist[u]+rcost, v);
      }
    }
    for (int v = 0; v < n; v++) if (pred[v] != -1) pi[v] += dist[v];
    return dist[t] < INF;
  }
  pair<F, C> augment(int s, int t, F cap) {
    F flow = cap;  C cost = 0;
    for (int v = t; v != s; v = pred[v])
      flow = min(flow, edges[pe[v]].cap - edges[pe[v]].flow);
    for (int v = t; v != s; v = pred[v])
      edges[pe[v]].flow += flow, edges[pe[v]^1].flow -= flow,
        cost += edges[pe[v]].cost * flow;
    return {flow, cost};
  }
  // Create a flow network on n vertices
  MinCostFlow(int n) : n(n), m(0), pred(n), pe(n), g(n) {}
  int add_edge(int u, int v, F cap, C cost) {
    edges.push_back({u, v, 0, cap, cost}), g[u].push_back(m++);
    edges.push_back({v, u, 0, 0, -cost}), g[v].push_back(m++);
    return m - 2;
  }
  edge &get_edge(int i) { return edges[i]; }
  pair<F, C> flow(int s, int t, F cap) {
    for (auto &e : edges) e.flow = 0;
    F flow = 0; C cost = 0;  bellman_ford();
    while (flow < cap && dijkstra(s, t)) {
      auto res = augment(s, t, cap - flow); flow += res.X, cost += res.Y;
    }
    return {flow, cost};
  }
};
//listings:/min_cost_flow

/*
// TESTED ON: UVA11301 - Dijkstra TLEd at 6 seconds. Leviticus passed in 0.14 seconds.
//listings:leviticus
// If Dijkstra's is too slow for min-cost flow, it can be substituted with
// the Leviticus algorithm. This has average complexity O(E * flow) on random
// graphs which is better than Dijkstra but is O(VE * flow) in the worst case.
bool levit(int s, int t) {
  vi id(n, 0);  dist.assign(n, INF);  dist[s] = 0; deque<int> q;  q.push_back(s);
  while (!q.empty()) {
    int v = q.front();  q.pop_front();	 id[v] = 2;
    for (auto i : g[v]) {
      auto& e = edges[i];
      if (e.flow < e.cap && dist[v] + e.cost < dist[e.to]) {
        dist[e.to] = dist[v] + e.cost;
        if (id[e.to] == 0) q.push_back(e.to);
        else if (id[e.to] == 2) q.push_front(e.to);
        id[e.to] = 1, pred[e.to] = v,	pe[e.to] = i;
      }
    }	
  }
  return dist[t] < INF;
}
//listings:/leviticus
*/

void solve_spoj_greedy() {
  int t;
  cin >> t;
  while (t--) {
    int N; cin >> N;
    int S = N, T = N + 1;
    MinCostFlow<int, int> G(N + 2);
    vi quantity(N);
    for (int i = 0; i < N; i++) {
      int x; cin >> x; x--;
      quantity[x]++;
    }
    for (int i = 0; i < N; i++) {
      G.add_edge(i, T, 1, 0);
      if (quantity[i] != 0) G.add_edge(S, i, quantity[i], 0);
    }
    int e; cin >> e;
    for (int i = 0; i < e; i++) {
      int x, y; cin >> x >> y; x--; y--;
      G.add_edge(x, y, N, 1);
      G.add_edge(y, x, N, 1);
    }
    auto res = G.flow(S, T, N);
    assert(res.first == N);
    cout << res.second << '\n';
  }
}

void solve_uva10594() {
  int N, M;
  while (scanf("%d %d", &N, &M) != EOF) {
    MinCostFlow<ll, ll> mcf(N);
    vi us(M), vs(M);
    vll costs(M);

    for (int i = 0; i < M; i++) scanf("%d %d %lld", &us[i], &vs[i], &costs[i]);

    ll D, K;
    scanf("%lld %lld", &D, &K);

    for (int i = 0; i < M; i++) {
      mcf.add_edge(us[i] - 1, vs[i] - 1, K, costs[i]);
      mcf.add_edge(vs[i] - 1, us[i] - 1, K, costs[i]);
    }

    auto ans = mcf.flow(0, N - 1, D);

    if (ans.first < D) puts("Impossible.");
    else printf("%lld\n", ans.second);
  }
}

const int src = 0, sink = 1;
int n, m;
int bank(int i) { return 2 + i; }
int popo(int j) { return 2 + n + j; }

const double EPS = 1e-8;
void solve_UVA10746() {
	while (cin >> n >> m, n > 0) {
    MinCostFlow<int, double> G(n + m + 2);
    for (int i = 0; i < n; i++) G.add_edge(src, bank(i), 1, 0);
    for (int j = 0; j < m; j++) G.add_edge(popo(j), sink, 1, 0);
    for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) {
      double t; cin >> t;
      G.add_edge(bank(i), popo(j), 1, t);
    }
    auto res = G.flow(src, sink, n);
		cout << fixed << setprecision(2) << (res.second + EPS) / n << '\n';
	}
}

int dist(int x1, int y1, int x2, int y2) {
	return abs(x1 - x2) + abs(y1 - y2);
}
void solve_spoj_baby() {
	while (cin >> n, n > 0) {
		vi a(n), b(n);
		
		for (int i = 0; i < n; i++) cin >> a[i];
		for (int i = 0; i < n; i++) cin >> b[i];
		
    MinCostFlow<int,int> G(2*n + 2);
    for (int i = 0; i < n; i++) {
        G.add_edge(src, 2 + i, 1, 0);
        G.add_edge(2 + n + i, sink, 1, 0);
    }
    
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
        G.add_edge(2 + i, 2 + n + j, 1, dist(i, a[i], j, b[j]));
			}
		}
		
		auto res = G.flow(src, sink, n);
    assert(res.first == n);
		cout << res.second << endl;
	}
}

char c(int k) {
	if (k < 26) return 'a' + k;
	else return 'A' + k - 26;
}
int idx(char a) {
	if ('a' <= a && a <= 'z') return a - 'a';
	else return (a - 'A') + 26;
}
void solve_CF491C() {
	int N, K;
	string enc, ans;
	cin >> N >> K >> enc >> ans;
	
	vvi A(K, vi(K, 0));
	for (int i = 0; i < N; i++)	A[idx(enc[i])][idx(ans[i])]--;
  
  MinCostFlow<int,int> G(2*K + 2);
  for (int i = 0; i < K; i++) {
    G.add_edge(src, 2 + i, 1, 0);
    G.add_edge(2 + K + i, sink, 1, 0);
  }
  
  for (int i = 0; i < K; i++) for (int j = 0; j < K; j++) 
      G.add_edge(2 + i, 2 + K + j, 1, A[i][j]);
	
	auto res = G.flow(src, sink, K);
	vi match(K);
  for (int i = 0; i < G.m; i++) 
    if (G.edges[i].flow == 1 && G.edges[i].from != src && G.edges[i].to != sink) 
      match[G.edges[i].from - 2] = G.edges[i].to - 2 - K;
  
	cout << -res.second << endl;
	for (int k = 0; k < K; k++) cout << c(match[k]); cout << endl;
}

int main() {
  //solve_uva10594();
  solve_spoj_greedy();
  //solve_UVA10746();
  //solve_spoj_baby();
  //solve_CF491C();
}