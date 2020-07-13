// Maximum bipartite matching using the Hopcroft-Karp algorithm.
// This is asymptotically faster than the alternating paths algorithm.
//
// Author      : Daniel Anderson
// Date        : 5/10/2016
// Reliability : 5
// Tested on   : UVA10092, UVA11138, UVA10080, UVA670, SPOJ-MATCHING
//
// Usage:
//  BipartiteMatching G(l, r)  create a graph with l left and r right nodes
//  add_edge(u, v)              add an edge from left node u to right node v
//  match()                     returns the size of the matching and the matches
//                              the matches indicate each left node's mate, or -1.
//
// Complexity: O(Sqrt(V) E)
#include<bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:matching
// Maximum unweighted bipartite matching using the Hopcroft-Karp algorithm.
// Returns the number of matches and a vector of each left node's match,
// or -1 if the node had no match.  Complexity: O(Sqrt(V) E)
struct BipartiteMatching {
	int L, R, p;  vi m, used, d; vvi adj; queue<int> q;
  bool bfs() {
    for (int v=0; v<R; v++) if (!used[v]) d[v] = p, q.push(v);
    while (!q.empty()) {
      int v = q.front(); q.pop();
      if (d[v] != d[R]) for (int u : adj[v]) if (d[m[u]] < p)
        d[m[u]] = d[v] + 1, q.push(m[u]);
    }
    return d[R] >= p;
  }
	int dfs(int v) {
    if (v == R) return 1;
    for (int u : adj[v]) if (d[m[u]] == d[v] + 1 && dfs(m[u])) return m[u] = v, 1;
    d[v] = d[R];  return 0;
	}
  // Create a Bipartite graph with L and R vertices in the left and right part
	BipartiteMatching(int L, int R) : L(L), R(R), d(R+1), adj(R) { }
	void add_edge(int u, int v) { adj[v].push_back(u); } // Add edge left(u) -> right(v)
	pair<int, vi> match() {   // Returns the number of matches and the matches for each
		int res = 0;  m.assign(L, R), used.assign(R+1, 0);       // node in the left part
		for (p=0; bfs(); p = d[R]+1) for (int v=0; v<R; v++)
      if (!used[v] && dfs(v)) used[v] = 1, res++;
    replace(m.begin(), m.end(), R, -1); return {res, m};
	}
};
//listings:/matching

// SPOJ - MATCHING
// Verdict: AC
void solve_SPOJ_MATCHING() {
	int N, M, P, A, B;
	scanf("%d%d%d", &N, &M, &P);
	BipartiteMatching G(N, M);
	for (int i = 0; i < P; i++) {
		scanf("%d%d", &A, &B);
		G.add_edge(A - 1, B - 1);
	}
	printf("%d\n", G.match().first);
}

// Verdict: AC
void solve_UVA10092() {
  int nk, np;
	while (cin >> nk >> np && nk) {
		vi category;
		for (int i = 0; i < nk; i++) {
			int num; cin >> num;
			for (int j = 0; j < num; j++)
				category.push_back(i);
		}
		int R = (int)category.size();
		
		BipartiteMatching G(np, R);
		
		for (int j = 0; j < np; j++) {
			int num; cin >> num;
			for (int i = 0; i < num; i++) {
				int cat; cin >> cat; cat--;
				for (int k = 0; k < R; k++) {
					if (category[k] == cat)
						G.add_edge(j, k);
				}
			}
		}
		
		auto res = G.match();
		if (res.first < R) cout << 0 << endl;
		else {
			vvi allocation(nk);
			for (int j = 0; j < np; j++) {
				if (res.second[j] != -1)
					allocation[category[res.second[j]]].push_back(j);
			}
			
			cout << 1 << endl;
			for (int i = 0; i < nk; i++) {
				for (int j = 0; j < (int)allocation[i].size(); j++) {
					cout << (allocation[i][j] + 1) << " \n"[j == (int)allocation[i].size() - 1];
				}
			}
		}
	}
}

// Verdict: AC
void solve_UVA11138() {
  int T, n, b;
	cin >> T;
	for (int t = 1; t <= T; t++) {
    	cin >> b >> n;
    	BipartiteMatching G(b, n);
    	int i;
    	for (int c = 0; c < b; c++) {
    		for (int d = 0; d < n; d++) {
    			cin >> i;
    			if (i) G.add_edge(c, d);
    		}
    	}
    	cout << "Case " << t << ": ";
    	cout << "a maximum of " << G.match().first << " nuts and bolts can be fitted together\n";
    }
}

// Verdict: AC
const double EPS = 1e-8;
#define SQR(x) (x)*(x)
double dist(double x, double y, pair<double,double> gopher) {
	return sqrt(SQR(x - gopher.first) + SQR(y - gopher.second));
}
void solve_UVA10080() {
	int n, m, s, v;
	while (cin >> n >> m >> s >> v) {
		BipartiteMatching G(n, m);
		vector<pair<double,double>> gophers(n);
		for (int i = 0; i < n; i++) {
			cin >> gophers[i].first >> gophers[i].second;
		}
		for (int j = 0; j < m; j++) {
			double x, y;
			cin >> x >> y;
			for (int i = 0; i < n; i++) {
				if (dist(x, y, gophers[i]) / v + EPS < s) {
					G.add_edge(i, j);
				}
			}
		}
		cout << (n - G.match().first) << endl;
	}
}

// Verdict: 
typedef pair<int, int> pii;
#define X first
#define Y second
#define sqr(x) (x)*(x)
double dist(pii a, pii b) {
	return sqrt(sqr(a.X - b.X) + sqr(a.Y - b.Y));
}
void solve_UVA670() {
	int L;  cin >> L;
	while (L--) {
		int N, M;
		cin >> N >> M;
		vector<pii> bob(N);
		for (int i = 0; i < N; i++) {
			cin >> bob[i].X >> bob[i].Y;
		}
		
		vector<pii> place(M);
		for (int i = 0; i < M; i++) {
			cin >> place[i].X >> place[i].Y;
		}
		
		BipartiteMatching G(N - 1, M);
		
		for (int i = 0; i < N - 1; i++) {
			for (int j = 0; j < M; j++) {
				if (dist(bob[i], place[j]) + dist(place[j], bob[i+1]) <= 2 * dist(bob[i], bob[i+1])) {
					G.add_edge(i, j);
				}
			}
		}
		
		auto res = G.match();
		
		vector<pii> route;
		for (int i = 0; i < N - 1; i++) {
			route.push_back(bob[i]);
			if (res.second[i] != -1) {
				assert(dist(bob[i], place[res.second[i]]) + dist(place[res.second[i]], bob[i+1]) <= 2 * dist(bob[i], bob[i+1]));
				route.push_back(place[res.second[i]]);
			}
		}
		route.push_back(bob.back());
		assert((int)route.size() == N + res.first);
		
		cout << N + res.first << '\n';
		for (int i = 0; i < (int)route.size(); i++) {
			if (i) cout << ' ';
			cout << route[i].X << ' ' << route[i].Y;
		}
		
		cout << '\n';
		if (L) cout << '\n';
	}
}

int main() {
  //solve_SPOJ_MATCHING();
  //solve_UVA10092();
  //solve_UVA11138();
  //solve_UVA10080();
  solve_UVA670();
}
