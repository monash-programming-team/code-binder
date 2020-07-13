// Geometry attempt again
//
// Author: Daniel (algorithms from Darcy's binder and elsewhere)
// Reliability: 5
// Tested on: SPOJ-QCJ4 (Minimum enclosing circle)
//            SPOJ-ALIENS (Minimum enclosing circle)
//            SWERC15J (Convex hull and point in convex polygon)
//            SPOJ-VCIRCLES (Area of union of circles)
//            SPOJ-CIRU (Area of union of circles)
//            CF600D (Area of intersection of two circles)
//            MONASH-C (Closest pair of points)
//            SPOJ-CLOPPAIR (Closest pair of points)
//            Brute-force closest pair of points
#include<bits/stdc++.h>
using namespace std;

#include "../code-template/debug.h"

#define X first
#define Y second

typedef long long int ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

/*
%%% NEEDED ::
% -- 2D GEOMETRY --
% * Common tangents to two circles
% * Convex cut
% * Convex polygon and circle intersection
% * Area of convex polygon intersection
% * Diameter of a convex polygon
% * Area of union of rectangles
% * KD-Tree
% * Delaunay triangulation
% * Voronori diagrams
% * Great circle distance
% * Lat-long conversion
% * Rotating calipers technique
% * Bentley Ottman algorithm (all intersection points of N lines)

% -- 3D GEOMETRY --
% * Point-to-line segment distance
% * Point-to-line distance
% * Point-to-triangle distance
% * Point-to-plane distance
% * Distance between line segments
% * Distance between lines
% * Rotate a point about a line
% * Ray-sphere intersection
% * Ray-plane intersection
% * 3D convex hull
% * Volume of convex polyhedron
*/

//listings:point
// Geometric constants
const double INF = 1e100;
const double EPS = 1e-8;
const double PI = acos(-1);

// Point geometry ------------------------------------------------------------------
#define x real()
#define y imag()
typedef complex<double> cd;
struct Point : public cd {
  Point() = default;  using cd::cd;
  Point(const cd& p) : cd(p) { }
  double& real() const { return (double&)*this; }
  double& imag() const { return ((double*)this)[1]; }
};
bool operator==(const Point& a, const Point& b) { return abs(a-b) <= EPS; }
bool operator<(const Point& a, const Point& b) { return a.x<b.x||(a.x==b.x&&a.y<b.y); }

double dot(const Point& a, const Point& b) { return (conj(a)*b).x; }    // Dot product
double det(const Point& a, const Point& b) { return (conj(a)*b).y; }  // Cross product
double dist(const Point& a, const Point& b) { return abs(a-b); }  // Euclidean distance
double dist2(const Point& a, const Point& b) { return norm(a-b); }  // Distance squared

// Rotate p anticlockwise by PI/2 radians (90 degrees) about (0,0)
Point inv(const Point& p) { return {-p.y, p.x}; }
// Rotate p anticlockwise by t radians about (0,0)
Point rotate(const Point& p, double t) { return p*polar(1.0, t); } 
// Rotate a anticlockwise by t radians about p
Point rotate(const Point& a, const Point& p, double t) { return rotate(a-p,t)+p; }
// Project the vector p onto the vector v
Point project(const Point& p, const Point& v) { return v * dot(p, v) / norm(v); }
// Orientation test (1: anticlockwise, -1: clockwise, 0: colinear)
int orientation(const Point& a, const Point& b, const Point& c) {
  double d=det(b-a,c-b); return d > EPS ? 1 : d < -EPS ? -1 : 0;
}
// Returns true if the points a, b, c are colinear
bool colinear(const Point& a, const Point& b, const Point& c) {
  return orientation(a,b,c) == 0;
}
// Returns the interior angle of ABC in [-PI,PI]
double angle(const Point& a, const Point& b, const Point& c) {
  return arg((a-b)/(c-b));
}
// Closest pair of points. Complexity: O(N log(N))
pair<Point,Point> closest_pair(vector<Point> P) {
  sort(P.begin(), P.end(), [](auto p1, auto p2) { return p1.y < p2.y; });
  set<Point> a; double d = INF; pair<Point,Point> cp{{0,0},{INF,0}};
  for (auto p = P.begin(), c = p; c != P.end(); c++) {
    while (p != c && p->y < c->y-d) a.erase(*p++);
    for (auto i=a.lower_bound(Point{c->x-d,0}); i != a.end() && i->x < c->x+d; i++)
      if (dist(*c, *i) < d) d = dist(*c, *i), cp = {*c, *i};
    a.insert(*c);
  }
  return cp;
}
//listings:/point

//listings:line
// Line geometry -------------------------------------------------------------------
struct Line { Point a, b; };

// The angle made by the line L with the x-axis in [-PI, PI]
double angle(const Line& L) { return arg(L.b - L.a); }

// Project the point p onto the line L
Point project(const Point& p, const Line& L) { return L.a + project(p-L.a, L.b-L.a); }
// Reflect the point p across the line L
Point reflect(const Point& p, const Line& L) {
  return L.a + conj((p-L.a)/(L.b-L.a)) * (L.b-L.a);
}
// Returns true if L1 and L2 are parallel
bool parallel(const Line& L1, const Line& L2) {
  return abs(det(L1.b-L1.a,L2.b-L2.a)) < EPS;
}
// Intersection point of two lines. Don't use for parallel lines.
Point intersection(const Line& L1, const Line& L2) {
  Point ab = L1.b-L1.a, qp = L2.a-L2.b, ap = L2.a-L1.a;
  double t = det(ab,ap)/det(ab,qp);
  return L1.a+t*ab;
}
// Distance between a point and a line
double dist(const Point& p, const Line& L) {
  return abs(det(L.b-L.a,p-L.a)/abs(L.b-L.a));
}
// Returns true if p is on the line L
bool point_on_line(const Point& p, const Line& L) { return dist(p, L) < EPS; }
//listings:/line

//listings:segment
// Segment geometry ----------------------------------------------------------------
struct Segment { Point a, b; };

// Returns true if the point p lies on the line segment seg
bool point_on_segment(const Point& p, const Segment& seg) {
  return abs(dist(seg.a,p) + dist(p,seg.b) - dist(seg.a,seg.b)) < EPS;
}
// Returns the closest point to p on the line segment seg
Point closest_point(const Point& p, const Segment& seg) {
  double t = dot(seg.b-seg.a,p-seg.a) / norm(seg.b);
  if (t > 1) return seg.b;
  else if (t < 0) return seg.a;
  else return seg.a + t * (seg.b - seg.a);
}
// Distance between a point and a line segment
double dist(const Point& p, const Segment& seg) {
  return dist(p, closest_point(p, seg));
}
// Returns the intersection of two non-parallel line segments. Returns {NAN,NAN}
// if there is no intersection. Check for NAN with isnan(), NOT ==.
Point intersection(const Segment& seg1, const Segment& seg2) {
  Point ab = seg1.b-seg1.a, qp = seg2.a-seg2.b, ap = seg2.a-seg1.a;
  double t = det(ab,ap)/det(ab,qp), s=det(ap,qp)/det(ab,qp);
  Point p = seg1.a + t * ab;
  if (t >= -EPS && t <= 1+EPS && s >= -EPS && s <= 1+EPS) return p;
  else return {NAN,NAN};
}
// Distance between two line segments. Can exclude conditional if segments are
double dist(const Segment& seg1, const Segment& seg2) {  //  guaranteed to not
  if (abs(det(seg1.b-seg1.a,seg2.b-seg2.a)) > EPS        //          intersect
    && !isnan(intersection(seg1,seg2).x)) return 0.0;
  return min({dist(seg1.a,seg2),dist(seg1.b,seg2),dist(seg2.a,seg1),dist(seg2.b,seg1)});
}
//listings:/segment

//listings:circle
// Circle geometry -----------------------------------------------------------------
struct Circle { Point c; double r; };
double area(const Circle& cir) { return PI*cir.r*cir.r; } // Area of a circle

// Returns true if the given point is contained within the given circle
bool point_in_circle(const Point& p, const Circle& c) {
  return dist(p,c.c) <= c.r + EPS;
}
// Returns true if the given point is on the boundary of the given circle
bool point_on_circle(const Point& p, const Circle& c) {
  return abs(dist(p,c.c)-c.r) <= EPS;
}
// Finds all intersection points of two non-coincident circles
vector<Point> intersection(const Circle& c1, const Circle& c2) {
  double d = dist(c2.c, c1.c);
  if (d > c1.r + c2.r + EPS || d < abs(c1.r - c2.r) - EPS) return {};  // None
  double a = (c1.r*c1.r - c2.r*c2.r + d*d) / (2*d);
  double h = sqrt(abs(c1.r*c1.r - a*a));
  Point p = c1.c + a/d * (c2.c - c1.c), q = h/d * inv(c2.c-c1.c);
  if (abs(h) < EPS) return {p};  // One intersection point
  else return {p + q, p - q};    // Two intersection points
}
// Compute the area of the intersection of two circles
double area_of_intersection(Circle c1, Circle c2) {
  if (c1.r < c2.r) swap(c1,c2);  double d = dist(c1.c,c2.c);
  if (d + c2.r <= c1.r + EPS) return area(c2);  // c2 is in c1
  if (d >= c1.r + c2.r - EPS) return 0.0;       // no overlap
  double alpha = acos((c1.r*c1.r+d*d-c2.r*c2.r)/(2*c1.r*d))*2;
  double beta = acos((c2.r*c2.r+d*d-c1.r*c1.r)/(2*c2.r*d))*2;
  double a1 = 0.5*beta*c2.r*c2.r-0.5*c2.r*c2.r*sin(beta);
  double a2 = 0.5*alpha*c1.r*c1.r-0.5*c1.r*c1.r*sin(alpha);
  return a1 + a2;
}
// Finds all intersection points of a circle and an infinite line
vector<Point> intersection(const Circle& c, const Line& L) {
  Point a = L.a - c.c, b = L.b - c.c, d = b - a;
  double dr = norm(d), D = det(a,b), desc = c.r*c.r * dr - D*D;
  if (desc < 0) return {};  // No intersection points
  if (abs(desc) < EPS) return { c.c-D/dr*inv(d) };  // One intersection
  double sgn = (d.y < -EPS ? -1 : 1);
  Point f = sgn*sqrt(desc)/dr * d;  d = c.c - D/dr * inv(d);
  return {d + f, d - f};  // Two intersection points
}
// Construct the circumcircle of the (non-colinear) points a, b, c
Circle circumcircle(const Point& a, const Point& b, const Point& c) {
  double g = 2*det(b-a,c-b); if (abs(g) < EPS) return {0,0};  // Colinear points
  double e = dot(b-a,b+a)/g, f = dot(c-a,c+a)/g;
  Point center = inv(f*(b-a)- e*(c-a));
  return {center, dist(a,center)};
}
// Construct a circle from two antipodal points on the boundary
Circle circle_from_diameter(const Point& a, const Point& b) {
  return {0.5*(a+b), dist(0.5*(a+b), a)};
}
// Find the smallest circle that encloses all of the given points. Complexity: O(N)
Circle minimum_enclosing_circle(vector<Point> P) {
  int N = (int)P.size(); random_shuffle(P.begin(), P.end());
  Circle c{P[0], 0};
  for (int i=1; i<N; i++) if (!point_in_circle(P[i], c)) {
    c = Circle{P[i],0};
    for (int j=0; j<i; j++) if (!point_in_circle(P[j], c)) {
      c = circle_from_diameter(P[i],P[j]);
      for (int k=0; k<j; k++) if (!point_in_circle(P[k], c))
        c = circumcircle(P[i],P[j],P[k]);
    }
  }
  return c;
}
// Find the area of the union of the given circles. Complexity: O(n^2 log(n))
double circle_union_area(const vector<Circle>& cir) {
  int n = (int)cir.size(); vector<bool> ok(n, 1); double ans = 0.0;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) if (i != j && ok[j]) 
    if (dist(cir[i].c,cir[j].c)+cir[i].r-cir[j].r < EPS) { ok[i] = false; break; }
  for (int i=0; i<n; i++) if (ok[i]) {
    bool flag = false;  vector<pair<double,double>> reg;
    for (int j=0; j<n; j++) if (i != j && ok[j]) {
      auto p = intersection(cir[i], cir[j]);
      if (p.size() < 2) continue; else flag = true;
      auto ang1 = arg(p[1] - cir[i].c), ang2 = arg(p[0] - cir[i].c);
      if (ang1 < 0) ang1 += 2*PI;
      if (ang2 < 0) ang2 += 2*PI;
      if (ang1 > ang2) reg.emplace_back(ang1, 2*PI), reg.emplace_back(0, ang2);
      else reg.emplace_back(ang1, ang2);
    }
    if (!flag) { ans += area(cir[i]); continue; }
    int cnt = 1; sort(reg.begin(), reg.end());
    for (int j=1; j<(int)reg.size(); j++)
      if (reg[cnt-1].Y >= reg[j].X) reg[cnt-1].Y = max(reg[cnt-1].Y, reg[j].Y);
      else reg[cnt++] = reg[j];
    reg.emplace_back(0,0);  reg[cnt] = reg[0];
    for (int j=0; j<cnt; j++) {
      auto p1 = cir[i].c + polar(cir[i].r, reg[j].Y);
      auto p2 = cir[i].c + polar(cir[i].r, reg[j+1].X);
      ans += det(p1, p2) / 2.0;
      double ang = reg[j+1].X - reg[j].Y;
      if (ang < 0) ang += 2*PI;
      ans += 0.5 * cir[i].r*cir[i].r * (ang - sin(ang));
    }
  }
  return ans;
}
//listings:/circle

//listings:polygon
// Polygon geometry ----------------------------------------------------------------
typedef vector<Point> Polygon;  // vertices of a polygon must be ordered

// Signed area of a polygon (positive if counter-clockwise, negative if clockwise)
double signed_area(const Polygon& P) {
  int n = (int)P.size();  double sum = 0.0;
  for (int i=0, j=n-1; i<n; j=i++) sum += det(P[i], P[j]);
  return sum / 2.0;
}

// Convex hull of the points P. colinear = true to include colinear points. Can
// return the upper and lower hulls separately if required. Complexity: O(N log(N))
Polygon convex_hull(vector<Point> P, bool colinear) {
  int n = (int)P.size(); auto eps = colinear ? EPS : -EPS;  vector<Point> L, U;
  sort(P.begin(), P.end()); P.erase(unique(P.begin(), P.end()), P.end());
  for (int i=0; i<n; i++) {
    while (L.size()>1 && det(P[i]-L.back(),L[L.size()-2]-P[i]) <= -eps) L.pop_back();
    while (U.size()>1 && det(P[i]-U.back(),U[U.size()-2]-P[i]) >= eps) U.pop_back();
    L.push_back(P[i]), U.push_back(P[i]);
  } // Return upper and lower hull here if desired
  auto hull = L; hull.insert(hull.end(), U.rbegin()+1, U.rend()-1);
  return hull;
}

// Point in polygon test. boundary = true if points on the boundary are considered
// to be in the polygon. Complexity: O(N)
bool point_in_polygon(const Polygon& poly, const Point& p, bool boundary) {
  int n = (int)poly.size(), i, j, c = 0;
  for (i=0, j=n-1; i<n; j=i++)
    if (poly[i] == p || colinear(poly[i],poly[j],p)) return boundary;
  for (i=0, j=n-1; i<n; j=i++) if (((poly[i].y <= p.y && p.y < poly[j].y) ||
    (poly[j].y <= p.y && p.y < poly[i].y)) && (p.x < (poly[j].x - poly[i].x)
    * (p.y - poly[i].y) / (poly[j].y - poly[i].y) + poly[i].x)); c = !c;
  return c;
}

// Point in polygon test for convex polygons. P must not contain colinear points.
// boundary = true if points on the boundary are considered to be in the polygon.
// Complexity: O(log(N))
bool point_in_convex_polygon(const Polygon& P, const Point& p, bool boundary) {
  int a = 1, b = (int)P.size()-1;
  if (point_on_segment(p,{P[a],P[0]}) || point_on_segment(p,{P[b],P[0]}))
    return boundary; else if (orientation(P[a],P[0],P[b]) > 0) swap(a,b);
  if (orientation(P[a],P[0],p) > 0 || orientation(P[b],P[0],p) < 0) return false;
  while (abs(a-b) > 1) {
    int c = (a+b)/2;
    if (orientation(P[c],P[0],p) > 0) b = c; else a = c;
  }
  return orientation(P[b],P[a],p) < 0 || (orientation(P[b],P[a],p)==0 && boundary);
}
//listings:/polygon

namespace problems {
  // Verdict: AC
  // Minimum enclosing circle problem
  namespace SPOJ_QCJ4 {
    void solve() {
      int n; cin >> n;
      vector<Point> points(n);
      for (auto& p : points) cin >> p.x >> p.y;
      auto mec = minimum_enclosing_circle(points);
      for (auto& p : points) assert(point_in_circle(p, mec));
      cout << setprecision(2) << fixed << mec.r * 2 << endl;
    }
  }
  
  // Verdict: AC
  // Minimum enclosing circle problem
  namespace SPOJ_ALIENS {
    void solve() {
      int c; cin >> c;
      while (c--) {
        int n; cin >> n;
        vector<Point> points(n);
        for (auto& p : points) cin >> p.x >> p.y;
        auto mec = minimum_enclosing_circle(points);
        for (auto& p : points) assert(point_in_circle(p, mec));
        cout << setprecision(2) << fixed << mec.r << '\n' << mec.c.x << ' ' << mec.c.y << '\n';
      }
    }
  }
  // Verdict: AC
  // Convex hull and point-in-convex-polygon problem
  namespace SWERC15J {
    void solve() {
      int L; cin >> L; vector<Point> large(L);
      for (auto& p : large) cin >> p.x >> p.y;
      int S; cin >> S; vector<Point> small(S);
      for (auto& p : small) cin >> p.x >> p.y;
      auto hull = convex_hull(large, false);
      int ans = 0;
      for (const auto& p : small) 
        if (point_in_convex_polygon(hull, p, true)) ans++;
      cout << ans << endl;
    }
  }
  // Verdict: AC
  // Area of union of circles problem
  namespace SPOJ_VCIRCLES {
    void solve() {
      int n; cin >> n;
      vector<Circle> circ(n);
      for (auto& c : circ) cin >> c.c.x >> c.c.y >> c.r;
      cout << fixed << setprecision(5) << circle_union_area(circ) << endl;
    }
  }
  // Verdict: AC
  // Area of union of circles problem
  namespace SPOJ_CIRU {
    void solve() {
      int n; cin >> n;
      vector<Circle> circ(n);
      for (auto& c : circ) cin >> c.c.x >> c.c.y >> c.r;
      cout << fixed << setprecision(3) << circle_union_area(circ) << endl;
    }
  }
  // Verdict: AC
  namespace CF600D {
    void solve() {
      Circle c1, c2;
      cin >> c1.c.x >> c1.c.y >> c1.r >> c2.c.x >> c2.c.y >> c2.r;
      cout << fixed << setprecision(20) << area_of_intersection(c1,c2) << endl;
    }
  }
  // Verdict: AC
  namespace MONASH_CPOP {
    void solve() {
      int n; cin >> n;
      vector<Point> pts(n);
      for (auto& p : pts) {
        int xx, yy; cin >> xx >> yy;
        p.x = xx, p.y = yy;
      }
      auto cpop = closest_pair(pts);
      cout << fixed << setprecision(20) << dist(cpop.X, cpop.Y) << endl;
    }
  }
  // Verdict: AC
  namespace SPOJ_CLOPPAIR {
    void solve() {
      int N; cin >> N;
      vector<Point> P(N);
      for (int i=0; i<N; i++) {
        double xx, yy; cin >> xx >> yy;
        P[i] = Point{(double)xx,(double)yy};
      }
      auto clop = closest_pair(P);
      vi idx;
      for (int i=0; i<N; i++) if (P[i] == clop.X || P[i] == clop.Y) idx.push_back(i);
      cout << idx[0] << ' ' << idx[1] << ' ' <<
        fixed << setprecision(6) << dist(clop.X, clop.Y) << endl;
    }
  }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  //problems::SPOJ_QCJ4::solve();
  //problems::SPOJ_ALIENS::solve();
  //problems::SWERC15J::solve();
  //problems::SPOJ_VCIRCLES::solve();
  //problems::SPOJ_CIRU::solve();
  //problems::CF600D::solve();
  problems::MONASH_CPOP::solve();
}

