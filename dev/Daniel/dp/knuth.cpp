// Knuth optimisation
//
// Author: Daniel
// Date: 01-02-2016
// Reliability: 1
// Tested on: UVA10304
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long int ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

const int INF = 1e8;
vector<int> psums;

ll cost(int i, int j) { return psums[j] - (i ? psums[i-1] : 0); }

//listings:knuth
// Knuth optimisation for problems of the form:
//   DP[i][j] = min(DP[i][k-1] + DP[k+1][j]) + C[i][j] for i <= k <= j
// The minimiser must be monotonic (satisfy opt[i][j-1] <= opt[i][j] <= opt[i+1][j])
// Alternatively, also applicable if instead the following conditions are met:
//   1. C[a][c] + C[b][d] <= C[a][d] + C[b][c] (quadrangle inequality)
//   2. C[b][c] <= C[a][d]                     (monotonicity)
// for all a <= b <= c <= d
// To use: -- define cost function cost(i,j) = C[i][j]
//         -- Compute base cases DP[i][i] for all 0 < i < n
//         -- Fill the rest of DP[i][j] = INF
//         -- Compute knuth(DP);
// Returns the optimal split points k = opt[i][j].  Complexity: O(N^2)
template<typename T> vvi knuth(vector<vector<T>>& DP) {
  int n = (int)DP.size();  vvi opt(n, vi(n));
  for (int i=0; i<n; i++) opt[i][i] = i;
  for (int len=1; len<n; len++) for (int i=0; i+len<n; i++) {
    int j = i + len;
    for (int k=opt[i][j-1]; k <= opt[i+1][j]; k++) {
      T new_cost = (k-1>=i ? DP[i][k-1] : 0) + (k+1<=j ? DP[k+1][j] : 0) + cost(i,j);
      if (new_cost < DP[i][j]) DP[i][j] = new_cost, opt[i][j] = k;
    }
  }
  return opt;
}
//listings:/knuth

namespace problems {
  namespace UVA10304 {
    void solve() {
      int n;
      while (cin >> n) {
        vi f(n); for (int i=0; i<n; i++) cin >> f[i]; psums.clear();
        partial_sum(f.begin(), f.end(), back_inserter(psums));
        vvi DP(n, vi(n, INF));
        for (int i=0; i<n; i++) DP[i][i] = f[i];
        auto opt = knuth(DP);
        cout << DP[0][n-1] - psums[n-1] << '\n';
      }
    }
  }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  problems::UVA10304::solve();  
}
