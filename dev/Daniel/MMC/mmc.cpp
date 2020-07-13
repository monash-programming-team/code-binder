// Minimum mean cycle in a directed graph
//
// Author: Daniel
// Date: 20-01-2017
// Reliability: 1
// Tested on: UVA11090
//
// Complexity: O(VE)
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:mmc
// Find the mean edge weight of the minimum mean cycle in a directed graph.
// If the graph contains no directed cycle, returns INF.  Complexity: O(VE)
struct MinimumMeanCycle {
	int n; vector<vector<pair<int,double>>> adj; const double INF = DBL_MAX / 2.0;
	MinimumMeanCycle(int N) : n(N), adj(n) { }
	void add_edge(int u, int v, double w) { adj[u].emplace_back(v,w); }
	double find_weight() {
		vector<vector<double>> DP(n+1, vector<double>(n, INF));
		fill(DP[0].begin(), DP[0].end(), 0);
		for (int i=0; i<n; i++) for (int u=0; u<n; u++) for (auto& e : adj[u])
			DP[i+1][e.X] = min(DP[i+1][e.X], DP[i][u] + e.Y);
		double res = INF;
		for (int i=0; i<n; i++) if (DP[n][i] < INF) {
			double hi = -INF;
			for (int j=0; j<n; j++) hi = max(hi, (DP[n][i]-DP[j][i]) / (n-j));
			res = min(res, hi);
		}
		return res;
	}
};
//listings:/mmc

// ----------------------------------------------------------------------------
//                                TEST PROBLEMS
// ----------------------------------------------------------------------------
namespace problems {
  // Verdict: AC
  namespace UVA_11090 {
    void solve() {
      int N; cin >> N;
      for (int t=1; t<=N; t++) {
        int n, m; cin >> n >> m;
        MinimumMeanCycle mmc(n);
        while (m--) {
          int a, b; double c; cin >> a >> b >> c;
          mmc.add_edge(a-1,b-1,c);
        }
        cout << "Case #" << t << ": ";
        double res = mmc.find_weight();
        if (res >= mmc.INF) cout << "No cycle found." << '\n';
        else cout << fixed << setprecision(2) << res << '\n';
      }
    }
  }
}

int main() {
	problems::UVA_11090::solve();
}
