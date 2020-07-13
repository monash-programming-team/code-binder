// Dynamic upper hull data structure
//
// Author: Daniel
// Date: 01-02-2016
// Reliability: 2
// Tested on: SPOJ-ACQUIRE, ICPC-WF11:F
//
#include<bits/stdc++.h>
using namespace std;

#include "../code-template/debug.h"

#define L first
#define W second
#define X first
#define Y second

typedef long long int ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:dynamic_hull
// Dynamic upper convex hull trick. Maintains the upper hull of a dynamic 
// set of lines. To maintain a lower hull, just negate the values.
// Complexity: O(log(N)) per operation. T = type of slope / intercept.
template<typename T> struct DynamicHull {
  struct Line {
    typedef typename multiset<Line>::iterator It;
    T a, b; mutable It me, endit, none;
    Line(T a, T b, It endit) : a(a), b(b), endit(endit) { }
    bool operator<(const Line& rhs) const {
      if (rhs.endit != none) return a < rhs.a;
      if (next(me) == endit) return 0;
      return (b - next(me)->b) < (next(me)->a - a) * rhs.a;
    }
  };
  multiset<Line> lines;
  void add_line(T a, T b) {
    auto bad = [&](auto y) {
      auto z = next(y);
      if (y == lines.begin()) {
        if (z == lines.end()) return false;
        return y->a == z->a && z->b >= y->b;
      }
      auto x = prev(y);
      if (z == lines.end()) return y->a == x->a && x->b >= y->b;
      return (x->b-y->b)*(z->a-y->a) >= (y->b-z->b)*(y->a-x->a);
    }; // WARNING: Change above comparison to doubles if you fear overflow
    auto it = lines.emplace(a, b, lines.end()); it->me = it;
    if (bad(it)) { lines.erase(it); return; }
    while (next(it) != lines.end() && bad(next(it))) lines.erase(next(it));
    while (it != lines.begin() && bad(prev(it))) lines.erase(prev(it));
  }
  T query(T x) {
    auto it = lines.lower_bound(Line{x,0,{}});
    return it->a * x + it->b;
  }
};
//listings:/dynamic_hull

ostream& operator<<(ostream& o, const DynamicHull<int>::Line& l) {
  if (l.endit == l.none) return o << "[ Query for x=" << l.a << "]";
  else return o << "[" << l.a << ", " << l.b << "]";
}

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
      DynamicHull<ll> hull; hull.add_line(-key_plots[0].W, 0);
      for (int i=1; i<N; i++) {
        hull.add_line(-key_plots[i].W, -DP[i-1]);
        DP[i] = -hull.query(key_plots[i].L);
      }
      cout << DP[N-1] << endl;
    }
  }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  problems::SPOJ_ACQUIRE::solve();
}
