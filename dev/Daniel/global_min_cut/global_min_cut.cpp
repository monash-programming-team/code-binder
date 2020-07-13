// Minimum cut in a weighted, undirected graph.
//
// Author: Daniel Anderson (Based on Stanford's ICPC notebook)
// Date: 21-01-2017
// Reliability: 0
// Tested on: UVA10989
//
// Complexity: O(V^3)
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:mincut
// Find the minimum cut in a weighted, undirected graph. Complexity: O(V^3)
template<typename T> struct MinCut {
  int n; vector<vector<T>> adj; const T INF = numeric_limits<T>::max();
  MinCut(int N) : n(N), adj(n, vector<T>(n)) { }
  void add_edge(int u, int v, T w) { adj[u][v] = adj[v][u] += w; }
  pair<T,vi> cut() {  // Returns the weight and the contents of one side of the cut
    T best = INF;  vi used(n), cut, best_cut;  auto weights = adj;
    for (int p=n-1; p >= 0; p--) {
      int prev, last = 0;  vi add = used, w = weights[0];
      for (int i=0; i<p; i++) {
        prev = last, last = -1;
        for (int j=1; j<n; j++) if (!add[j] && (last==-1 || w[j]>w[last])) last = j;
        if (i == p-1) {
          for (int j=0; j<n; j++) weights[prev][j] += weights[last][j];
          for (int j=0; j<n; j++) weights[j][prev] = weights[prev][j];
          used[last] = 1, cut.push_back(last);
          if (w[last] < best) best = w[last], best_cut = cut;
        } else {
          for (int j=0; j<n; j++) w[j] += weights[last][j];
          add[last] = 1;
        }
      }
    }
    return {best, best_cut};
  }
};
//listings:/mincut

// ----------------------------------------------------------------------------
//                                TEST PROBLEMS
// ----------------------------------------------------------------------------
namespace problems {
  // Verdict: 
  namespace UVA_10989 {
    void solve() {
      int N; cin >> N;
      for (int t=1; t<=N; t++) {
        cout << "Case #" << t << ": ";
        int n, m; cin >> n >> m;
        MinCut<int> G(n);
        while (m--) {
          int a, b, c; cin >> a >> b >> c;
          G.add_edge(a-1,b-1,c);
        }
        cout << G.cut().X << '\n';
      }
    }
  }
}

int main() {
  problems::UVA_10989::solve();
}
