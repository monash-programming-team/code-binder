#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef vector<int> vi;

// Fenwick Tree
//
// Author      : Peter Whalan (modified by Daniel)
// Date        : February 29, 2016
// Reliability : 5
// Tested On   : CF:634/C, CERC 12 G, Monash group, ASC18:E, Brute-force
//
// Performs range sum queries on array. Allows individual values to be
// modified.
// 
// Treats arguments as 0-based
//
// Complexity:
//		O( N ) to build
//		O( log N ) to update and query

//listings:ft
// Fenwick Tree with ranged queries and point updates. Complexity: O(log(n))
template<typename T> struct FenwickTree {
	int N;  vector<T> A;
	FenwickTree(int n): N(n+1), A(N) {}                    // Create tree with n elements
	T rq(int b) { T r=0; for (;b;b-=b&-b) r+=A[b]; return r; }          // Get sum A[0,b)
	T rq(int a,int b) { return rq(b)-rq(a); }                           // Get sum A[a,b)
	void adjust(int i, T v) { for (i++;i<N;i+=i&-i) A[i]+=v; }               // A[i] += v
  int lower_bound(T sum) {                  // find min i such that sum(A[0..i]) >= sum
    int i = 0;                                       // Returns n if there is no such i
    for (int b = 1 << (31-__builtin_clz(N)); b; b /= 2)     // (Only works if A[i] >= 0
      if (i+b < N && sum > A[i+b]) sum -= A[i+b], i+=b;                   // for all i)
    return i;
  }
};
//listings:/ft

void test_lower_bound() {
  srand(time(0));
  int testno = 1;
  while(1) {
    cout << "Test " << testno++ << "           \r" << flush;
    int N = 1000;
    FenwickTree<int> a(N);
    for (int i=0; i<N; i++) a.adjust(i, rand()%5);
    int total_sum = a.rq(0,N+1);
    for (int sum=0; sum<total_sum+1; sum++) {
      int lb = 0;
      while (a.rq(0,lb+1) < sum) lb++;
      int ans = a.lower_bound(sum);
      assert(ans == lb);
    }
  }
}

namespace problems {
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
      
      FenwickTree<int> ppl(m+n+2+10);
      auto insert = [&](int i) { ppl.adjust(i, 1); };
      auto erase = [&](int i) { ppl.adjust(i, -1); };
      auto find_by_order = [&](int k) { return ppl.lower_bound(k+1); };
      
      for (int i=0; i <= m + n + 1; i++) insert(i);
      
      int loc = 0;
      while ((total_size > 1 && ppl.lower_bound(1) != 0) || total_size > 2) {
        loc = (loc + k - 1) % total_size;
        int dead_guy = find_by_order(loc);
        if (dead_guy == 0 || dead_guy == m+1){
          erase(0); erase(m+1);
          total_size -= 2;
          if (dead_guy == m+1) loc--;
        }
        else {
          erase(dead_guy);
          total_size--;	
        }
      }
      int survivor = find_by_order(0);
      cout << id(survivor) << endl;
    }
  }
  
  namespace monash {
    void solve_monash() {
      ios::sync_with_stdio(0); cin.tie(0);
      
      int n; cin >> n;
      FenwickTree<int> ft(n);
      for (int i=0; i<n; i++) {
        int a; cin >> a; ft.adjust(i, a);
      }
      int m; cin >> m;
      int caseno = 1;
      while (m--) {
        int type; cin >> type;
        if (type == 1) {
          int i, v; cin >> i >> v;
          ft.adjust(i-1, v);
          cout << caseno++ << '\n';
        }
        else {
          int i, j; cin >> i >> j;
          cout << caseno++ << ' ' << ft.rq(i-1,j) << '\n';
        }
      }
    }
  }
}

int main() {
  //test_lower_bound();
  //problems::infinity_sect::solve();
  
  FenwickTree<int> ft(10);
  for (int i=0; i<10; i++) ft.adjust(i,1);
  cout << ft.lower_bound(9) << ' ' << ft.lower_bound(10) << ' ' << ft.lower_bound(100) << endl;
}

