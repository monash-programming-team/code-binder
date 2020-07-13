// Josephus problem
//
// Author: Daniel Anderson
// Date: 18-01-2017
// Reliability: 5
// Tested on: UVA11351, UVA10774, CodeChef-IGNUS15B, SPOJ-ANARC08H,
//            SPOJ-DANGER, IEEEXtreme10-Flower Games
// 
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;

//listings:josephus
// Josephus Problem (0-based): k=2 special case. Complexity: O(1)
ll survivor(ll n) { return (n - (1LL << (63 - __builtin_clzll(n)))) * 2; }

// Josephus Problem (0-based): Determine the survivor. Complexity: O(n)
int survivor(int n, int k) {
  vi A(n+1);  // A[i] is the survivor with i people, killing every k'th
  for (int i=2; i<=n; i++) A[i] = (A[i-1]+(k%i))%i;
  return A[n];  // OPTIONAL: Return entire array if multiple values needed
}
//listings:/josephus

#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;

typedef tree<int, null_type, less<int>, rb_tree_tag,
		tree_order_statistics_node_update> ordered_set;

// Determine the killing order in O(n log(n))
vi fast_killing_order(int n, int k) {
  ordered_set ppl;  vi ans;
  for (int i=0; i<n; i++) ppl.insert(i);
  int pos = 0;
  while (!ppl.empty()) {
    pos = (pos + k - 1) % ppl.size();
    auto it = ppl.find_by_order(pos);
    ans.push_back(*it);
    ppl.erase(it);
  }
  return ans;
}

// ----------------------------------------------------------------------------
//                               TEST PROBLEMS
// ----------------------------------------------------------------------------

// Verdict: AC
void solve_UVA11351() {
  int k; cin >> k;
  for (int t=1; t<=k; t++) {
    int n, k; cin >> n >> k;
    cout << "Case " << t << ": " << survivor(n,k)+1 << '\n';
  }
}

// Verdict: AC
void solve_UVA10774() {
  int T; cin >> T;
  for (int t=1; t<=T; t++) {
    int n; cin >> n;
    int winner = -1, rounds = 0;
    while (1) {
      rounds++;
      winner = survivor(n);
      if (winner == n - 1) break;
      n = winner + 1;
    }
    cout << "Case " << t << ": " << (rounds - 1) << ' ' << (winner + 1) << '\n'; 
  }
}

// Verdict: AC
void solve_CC_IGNUS15B() {
  int T; cin >> T;
  while (T--) {
    ll n; cin >> n;
    cout << survivor(n)+1 << '\n';
  }
}

// Verdict: AC
void solve_SPOJ_ANARC08H() {
  int N, D;
  while (cin >> N >> D && N && D) {
    cout << N << ' ' << D << ' ' << survivor(N,D)+1 << '\n';
  }
}

// Verdict: AC
void solve_SPOJ_DANGER() {
  int d,z;
  while (scanf("%de%d", &d, &z) && d) {
    int n = d; while (z--) n *= 10;
    cout << survivor(n)+1 << '\n';
  }
}

// Verdict: AC
void solve_IEEE10_flower() {
  int T; cin >> T;
  while (T--) {
    ll N; cin >> N;
    cout << survivor(N)+1 << '\n';
  }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  //solve_UVA11351();
  //solve_UVA10774();
  //solve_CC_IGNUS15B();
  //solve_SPOJ_ANARC08H();
  //solve_SPOJ_DANGER();
  solve_IEEE10_flower();
}

