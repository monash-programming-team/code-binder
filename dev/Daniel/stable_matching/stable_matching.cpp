// Stable matching problem
//
// Author: Daniel Anderson
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

//listings:stable_matching
// Finds a stable matching with the given preferences. mpref[i] lists male i's
// preferred matches in order (highest first). fpref lists female preferences
// in the same format. Returns a list of each male's match. Complexity: O(n^2)
vi stable_matching(const vvi& mpref, const vvi& fpref) {
  int n = (int)mpref.size(); vi mpair(n,-1), fpair(n,-1), p(n); vvi forder(n, vi(n));
  for (int i=0;i<n; i++) for (int j=0; j<n; j++) forder[i][fpref[i][j]] = j;
  for (int i=0; i<n; i++) {
    while (mpair[i] < 0) {
      int w = mpref[i][p[i]++], m = fpair[w];
      if (m == -1) mpair[i] = w, fpair[w] = i;
      else if (forder[w][i] < forder[w][m])
        mpair[m] = -1, mpair[i] = w, fpair[w] = i, i = m;
    }
  }
  return mpair;
};
//listings:/stable_matching

// ----------------------------------------------------------------------------
//                                TEST PROBLEMS
// ----------------------------------------------------------------------------
namespace problems {
    namespace SPOJ_STABLEMP {
      void solve() {
        int T; cin >> T;
        while (T--) {
          int n; cin >> n;
          vvi wp(n, vi(n)), mp(n, vi(n));
          for (int i=0; i<n; i++) {
            int w; cin >> w; w--;
            for (int j=0; j<n; j++) { cin >> wp[w][j]; wp[w][j]--; }
          }
          for (int i=0; i<n; i++) {
            int m; cin >> m; m--;
            for (int j=0; j<n; j++) { cin >> mp[m][j]; mp[m][j]--; }
          }
          auto sol = stable_matching(mp, wp);
          for (int i=0; i<n; i++) cout << (i+1) << ' ' << (sol[i]+1) << '\n';
        }
      }
    }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  problems::SPOJ_STABLEMP::solve();
}
