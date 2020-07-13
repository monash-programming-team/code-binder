// Divide and conquer DP
//
// Author: Daniel
// Date: 01-02-2016
// Reliability: 1
// Tested on: HackerRank - Guardians of the Lunatics
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long int ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

const ll INF = 1e18;

int L, G;
vector<ll> C, psum;
vector<vector<ll>> DP;

ll cost(int i, int j) { return (j - i) * (psum[j-1] - (i ? psum[i-1] : 0)); }

//listings:dc
// Divide and conquer optimisation for dynamic programs of the form
//   DP[i][j] = min(DP[i-1][k] + C[k][j]) for k < j.  i <= K, j <= N
// The minimiser must be monotonic (satisfy opt[i][j] <= opt[i][j+1]).
// To use: -- define cost function cost(k,j) = C[k][j]
//         -- fill base cases DP[0][0], DP[0][j], DP[i][0]
//         -- fill the rest DP[i][j] = INF
//         -- compute each row: for(int i=1; i<=K; i++) compute(i,1,N,0,N)
// Complexity: O(KN log(N))
void compute(int i, int l, int r, int optL, int optR) {
  if (r < l) return;  int mid = (l + r) / 2, opt = optL;
  for (int k=optL; k<=min(mid-1,optR); k++) {
    ll new_cost = DP[i-1][k] + cost(k,mid);
    if (new_cost < DP[i][mid]) DP[i][mid] = new_cost, opt = k;
  }
  compute(i, l, mid-1, optL, opt), compute(i, mid+1, r, opt, optR);
}
//listings:/dc

int main() {
  cin >> L >> G;
  C.resize(L);
  for (auto& x : C) cin >> x;
  partial_sum(C.begin(), C.end(), back_inserter(psum));
  DP.assign(G+1,vector<ll>(L+1, INF)); 
  DP[0][0] = 0; for (int g=1; g<=G; g++) DP[g][0] = 0;
  for (int g=1; g<=G; g++) compute(g,1,L,0,L);
  cout << DP[G][L] << endl;
}
