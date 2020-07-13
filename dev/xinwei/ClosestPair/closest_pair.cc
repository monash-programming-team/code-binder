#include <bits/stdc++.h>

using namespace std;
const double EPS = 1e-8;
const double INF = 1e16;
bool dEqual(double a, double b) { return fabs(a - b) < EPS; }
#define sq(x) (x)*(x)

struct Point {
  double x, y;
  bool operator<(const Point &other) const {
    return y < other.y || (y == other.y && x < other.x);
  }
};

vector<Point> P;
double dist(Point p, Point q){
  return sqrt( sq(p.x - q.x) + sq(p.y - q.y) );
}

bool cmp_x(const Point &p, const Point &q){
  return p.x < q.x || (p.x == q.x && p.y < q.y);
}

int main(){
  std::ios_base::sync_with_stdio(false);

  int n; 
  scanf("%d", &n);
  P.resize(n);
  for (int i = 0; i < n; i++)
    scanf("%lf%lf", &P[i].x, &P[i].y);

  sort(P.begin(), P.end(), cmp_x);
  set<Point> S;
  double d = INF;
  int ptr = 0;
  for (int i = 0; i < n; i++){
    while (ptr < i && (P[i].x - P[ptr].x) >= d){
      S.erase(P[ptr]);
      ptr++;
    }
    auto lb = S.lower_bound( {INF, P[i].y - d} );
    auto ub = S.upper_bound( {INF, P[i].y + d} );
    for (auto it = lb; it != ub; ++it){
      double cur = dist(P[i], *it);
      if (cur < d)
	d = cur;
    }
    S.insert(P[i]);
  }
  cout << fixed << setprecision(10) << d << endl;

  return 0;
}
