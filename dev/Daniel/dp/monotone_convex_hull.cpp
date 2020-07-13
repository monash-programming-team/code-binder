// Monotone upper hull data structure
//
// Author: Daniel
// Date: 01-02-2016
// Reliability: 2
// Tested on: SPOJ-ACQUIRE, ICPC-WF11:F
//
#include<bits/stdc++.h>
using namespace std;

#define L first
#define W second

typedef long long int ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:monotone_hull
// Monotonic convex hull trick. Find the maximum value on the upper envelope
// of a dynamic set of lines such that the queries and slopes are monotonically
// non-decreasing. To maintain a lower hull, just negate the values.
// Complexity: amortized O(1) per operation. T = type of slope/intercept.
template<typename T> struct MonotoneHull {
  struct Line { T a, b; double x; };  deque<Line> lines;
  void add_line(T a, T b) { // Add a line of the form y = ax + b
    double x = -1e200;
    while (!lines.empty()) {
      if (a == lines.back().a) x *= b < lines.back().b ? -1 : 1;
      else x = 1.0 * (lines.back().b - b) / (a - lines.back().a);
      if (x < lines.back().x) lines.pop_back();
      else break;
    }
    lines.push_back({a,b,x});
  }
  T query(T x) {     // Find min A[i]x + B[i]. Can alter this to binary search if the
    while (lines.size() > 1 && lines[1].x <= x) lines.pop_front();    // query points
    return lines[0].a * x + lines[0].b;  // are not monotone but the slopes still are
  }
};
//listings:/monotone_hull

// Code for comparison
const int MAXN = 100100;
struct Hull {
    long long a[MAXN], b[MAXN];    double x[MAXN];    int head, tail;
    Hull(): head(1), tail(0) {}
    long long get(long long xQuery) {
        if (head > tail) return 0;
        while (head < tail && x[head + 1] <= xQuery) head++;
        x[head] = xQuery;
        return a[head] * xQuery + b[head];
    }
    void add(long long aNew, long long bNew) {
        double xNew = -1e18;
        while (head <= tail) {
            if (aNew == a[tail]) return;
            xNew = 1.0 * (b[tail] - bNew) / (aNew - a[tail]);
            if (head == tail || xNew >= x[tail]) break;
            tail--;
        }
        a[++tail] = aNew; b[tail] = bNew; x[tail] = xNew;
    }
};


namespace problems {
  namespace SPOJ_ACQUIRE {
    void solve() {
      int N; cin >> N;
      vector<pair<ll,ll>> plot(N);
      for (int i=0; i<N; i++) cin >> plot[i].L >> plot[i].W;
      sort(plot.begin(), plot.end());
      vi smax(N); iota(smax.begin(), smax.end(), 0);
      for (int i=N-2; i >= 0; i--)
        if (plot[i].W <= plot[smax[i+1]].W) smax[i] = smax[i+1];
      vector<pair<ll,ll>> key_plots;
      for (int i=0; i<N; i++) if (smax[i] == i) key_plots.push_back(plot[i]);
      N = (int)key_plots.size();
      vector<ll> DP(N); DP[0] = key_plots[0].L * key_plots[0].W;
      MonotoneHull<ll> hull; hull.add_line(key_plots[0].W, 0);
      for (int i=1; i<N; i++) {
        hull.add_line(key_plots[i].W, DP[i-1]);
        DP[i] = hull.query(key_plots[i].L);
      }
      cout << DP[N-1] << endl;
    }
  }
}

// Test against other code binder code
namespace brute_force_tests {
  const int T = 1000000;
  const int Q = 1000;
  const int MAXA = 10;
  const int MAXB = 10;
  const int MAXADD = 2;
  void test() {
    srand(time(0));
    for (int t=1; t<=T; t++) {
      cout << "Test " << t << "/" << T << "        \r" << flush;
      MonotoneHull<ll> hull1;
      Hull hull2;
      vector<ll> queries; vector<pair<ll,ll>> ab;
      ll x = 0, a = -1;
      for (int q=1; q<=Q; q++) {
        if (rand()&1 || ab.empty()) {  // insert
          a -= (rand() % MAXA + 1);
          ll b = rand() % MAXB;
          ab.emplace_back(a,b);
          hull1.add_line(a,b);
          hull2.add(a,b);
        } else {         // query
          x += rand() % MAXADD;
          queries.push_back(x);
          ll val1 = hull1.query(x);
          ll val2 = hull2.get(x);
          if (val1 != val2) {
            assert(val1 == val2);
          }
        }
      }
    }
  }
}

int main() {
  //problems::SPOJ_ACQUIRE::solve();
  brute_force_tests::test();
}
