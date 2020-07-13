// Misof Tree - a lightweight order statistics interval tree
//
// Source: http://codeforces.com/contest/374/submission/5568759 (modified by Daniel)
// Reliability: 1
// Tested on: ASC18-E
//
// Complexity: O(log(N)) time per query, O(4*N) space.

#include<bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;

//listings:misof
// Maintains minimal order statistics on a multiset of elements in 0 to N-1.
// If you need to store larger elements, compress the input first.
// Complexity: O(log(N)) time per query, O(4*N) space.
struct MisofTree {
  int n, leaf, size=0; vi a, s, tree;
  MisofTree(int N) : n(N), leaf(1<<(32-__builtin_clz(n))), a(N), s(N), tree(100+4*N) {}
  void insert(int x) { for (int i=leaf+x; i; i /= 2) ++tree[i]; size++; }
  void erase(int x) { for (int i=leaf+x; i; i /= 2) --tree[i]; size--; }
  int find_by_order(int k) {  // Return the k'th smallest element in the multiset
    if (++k > tree[1] || k <= 0) return -1;  // indicates no result
    int i; for (i=1; i<leaf; i*=2) if (k > tree[i]) k -= tree[i++];
    return i - leaf + (k > tree[i]);
  }
};
//listings:/misof

namespace problems {
  // ASC18: Problem E: Infinity Sect
  namespace infinity_sect {
    ll m, n, k;

    int id(int loc) {
      if (loc <= m) return loc;
      else if (loc == m+1) return 0;
      else return -(loc - (m+1));
    }

    void solve() {
    #ifdef ONLINE_JUDGE
      freopen("infinity.in", "r", stdin);
      freopen("infinity.out", "w", stdout);
    #endif
      ios::sync_with_stdio(0); cin.tie(0);
      
      cin >> m >> n >> k;
      ll total_size = m + n + 2;
      
      MisofTree ppl(m+n+2+10);
      for (int i=0; i <= m + n + 1; i++) ppl.insert(i);
      
      int loc = 0;
      while ((total_size > 1 && ppl.find_by_order(0) != 0) || total_size > 2) {
        loc = (loc + k - 1) % total_size;
        int dead_guy = ppl.find_by_order(loc);
        if (dead_guy == 0 || dead_guy == m+1){
          ppl.erase(0), ppl.erase(m+1);
          total_size -= 2;
          if (dead_guy == m+1) loc--;
        }
        else {
          ppl.erase(dead_guy);
          total_size--;	
        }
      }
      int survivor = ppl.find_by_order(0);
      cout << id(survivor) << endl;
    }
  }
}

int main() {
  problems::infinity_sect::solve();
}