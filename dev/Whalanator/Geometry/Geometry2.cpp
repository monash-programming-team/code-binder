#include<bits/stdc++.h>

using namespace std;

typedef long long ll;

//listings:geometry
#define x real()
#define y imag()
#define X first
#define Y second
#define cpt const Pt&

const double EPS = 1e-9;
const double INF = 1e15;
const double pi=acos(-1);
bool deq(double a,double b) {return abs(a-b)<EPS;}

typedef complex<double> cpx;

struct Pt : public cpx {
    Pt() = default;
	using cpx::cpx;
	Pt(cpx a) : cpx(a) {}
    double& x const { return (double&)*this; }
    double& y const { return ((double*)this)[1]; }
	bool operator ==(cpt b) const {return abs(*this-b) < EPS; }
	bool operator <(cpt b) const {return x<b.x || (x==b.x && y<b.y); }
};
//listings:/geometry

//Allow points to be read in by input streams
istream& operator >>(istream& is, Pt& p) {
	return is >> p.x >> p.y;
}

double dot(cpt a, cpt b) {return (conj(a) * b).x;}// Dot product
double det(cpt a, cpt b) {return (conj(a) * b).y;}//Determinant/"Cross Product"
double angle(cpt a, cpt b) {return arg(b - a);}// [-pi,pi] a to b with x axis
double angle (cpt a, cpt b, cpt c) {return arg((a-b)/(c-b));}//[-pi,pi]
//double slope(cpt a, cpt b) {return tan(arg(b - a));}// m for line segment (a,b)
double slope(cpt a, cpt b) {return (b.y-a.y)/(b.x-a.x);}

Pt rotate(cpt a, double theta) {return a * polar((double)1.0, theta);}//anticlockwise
//around p by theta anticlockwise
Pt rotate(cpt a, cpt p, double theta) {return rotate(a - p,theta) + p;}
Pt project(cpt p, cpt v) {return v * dot(p, v) / norm(v);}// p onto v
Pt project(cpt p, cpt a, cpt b) {return a+project(p-a,b-a);}//p onto line (a,b)
//reflect p across the line (a,b)
Pt reflect(cpt p, cpt a, cpt b) {return a + conj((p - a) / (b - a)) * (b - a);}

//bool colinear(Pt a, Pt b, Pt c) {return deq(det(b-a,c-b),0);}

// Orientation test (1 anticlockwise, -1 clockwise, 0 colinear)
int orient(cpt a, cpt b, cpt c) {
	double d=det(b-a,c-b);
	return d>EPS?1:d<-EPS?-1:0;
}

//Compare points by principal argument (-pi,pi] breaking ties by norm.
//0 is considered less than everything else.
//untested
bool argcompapprox(cpt a,cpt b) {
	if (b==0) return 0;
	if (a==0) return 1;
	double a1=arg(a),a2=arg(b);
	if (a1<-pi+EPS/2) a1+=2*pi;
	if (a2<-pi+EPS/2) a2+=2*pi;
	return a1+EPS<a2 || (deq(a1,a2) && norm(a)<norm(b));
}

//untested
bool argcompexact(cpt a,cpt b) {
	if (b==0) return 0;
	if (a==0) return 1;
	bool r1=a.y>0 || (a.y==0 && a.x<0);
	bool r2=b.y>0 || (b.y==0 && b.x<0);
	ll d=det(a,b);
	return r1<r2 || (r1==r2 && (d>0 || (d==0 && norm(a)<norm(b))));
}

// Point on line segment (including endpoints)
bool ptonseg(cpt a, cpt b, cpt p) {
	Pt u=b-a,v=p-a;
	return a==p || b==p ||
		((0 < dot(u,v) && dot(u,v) < norm(u)) && deq(det(u,v),0));
}

// Signed area of polygon
// Positive for anticlockwise orientation
double polygonarea(const vector<Pt>& p) {
	double r=0;
	int n=p.size();
	for (int j=0,i=n-1;j<n;i=j++) r+=det(p[i],p[j]);
	return r/2;
}

// Convex hull O(NlogN)
// if all points are colinear the middle points come up twice forwards and
// backwards e.g. a-b-c-d becomes a-b-c-d-c-b
// To remove colinear points change <-EPS and >EPS to <EPS and >-EPS.
vector<Pt> convexhull(vector<Pt> p) {
  sort(p.begin(),p.end()); p.resize(unique(p.begin(),p.end())-p.begin());
  int l=0,u=0;
  vector<Pt> L(p),U(p);
  if (p.size()<=2) return p;
  for (Pt& i:p) {
    while (l>1 && det(i-L[l-1],L[l-2]-i)<-EPS) l--;
    while (u>1 && det(i-U[u-1],U[u-2]-i)>EPS) u--;
    L[l++]=U[u++]=i;
  }
  L.resize(l+u-2);
  copy(U.rend()-u+1,U.rend()-1,L.begin()+l);

  return L;
}

//Point in polygon test O(N)
//Returns:
// 0 if not in polygon
// 1 if on boundary
// 2 if in interior
int ptinpoly(const vector<Pt>& p, cpt q) {
	int n=p.size();
	int i,j,r=0;
	for (j=0,i=n-1;j<n;i=j++) { //that trick to avoid modding
		if (ptonseg(p[i],p[j],q)) return 1;
		if (((p[i].y <= q.y && q.y < p[j].y)
			|| (p[j].y <= q.y && q.y < p[i].y))
			&& q.x < (p[j].x-p[i].x) * (q.y-p[i].y)/(p[j].y-p[i].y) + p[i].x)
			r^=2;
	}
	return r;
}

// Point in polygon test for convex polygons. P must not contain colinear points.
// boundary = true if points on the boundary are considered to be in the polygon.
// Complexity: O(log(N))
bool point_in_convex_polygon(const vector<Pt>& P, const Pt& p, bool boundary) {
  int a = 1, b = (int)P.size()-1;
  if (ptonseg(P[a],P[0],p) || ptonseg(P[b],P[0],p))
    return boundary; else if (orient(P[a],P[0],P[b]) > 0) swap(a,b);
  if (orient(P[a],P[0],p) > 0 || orient(P[b],P[0],p) < 0) return false;
  while (abs(a-b) > 1) {
    int c = (a+b)/2;
    if (orient(P[c],P[0],p) > 0) b = c; else a = c;
  }
  return orient(P[b],P[a],p) < 0 || (orient(P[b],P[a],p)==0 && boundary);
}

Pt solve(cpt a, cpt b, cpt v) {// solves [a b]x==v with Cramer's rule.
	return Pt(det(v,b)/det(a,b),det(a,v)/det(a,b));
}

//Intersection of 2 line segments. Divides by 0 if they are parallel.
//Returns {nan,nan} if they don't intersect.
//Uncomment if statements below to get infinite lines.
Pt intersectline(Pt a, Pt b, Pt p, Pt q) {
	Pt ab=b-a,qp=p-q,ap=p-a;
	double s=det(ap,qp)/det(ab,qp),t=det(ab,ap)/det(ab,qp);
	//double s,t;
	//tie(s,t)=solve(b-a,p-q,p-a);
	//a+t(b-a)=p+s(q-p)
	
	//Can also just use ptonseg.
	if (-EPS<s && s<1+EPS //Answer is on ab
		&& -EPS<t && t<1+EPS) //Answer is on pq 
		return a+s*ab;
	return Pt(NAN,NAN);
}

//Distance between infinite line and point.
double distlinept(cpt a, cpt b, cpt p) {
	return abs(det(b-a,p-a)/abs(b-a));
}

//Distance between finite line and point
double distfinitelinept(Pt a, Pt b, Pt p) {
	b-=a;p-=a;
	double sp=(p/b).x;//dot(b,p)/norm(b);
	Pt closest;
	if (sp>=0) {
		if (sp>1) closest=b;
		else closest=sp*b;
	}
	return abs(closest-p); // Note that actual closest Pt on line is closest + a
}

//Are lines perpendicular
bool areperp(cpt a,cpt b,cpt p,cpt q) {
	return deq(dot(b-a,q-p),0);
}

//Are lines parallel?
bool arepara(cpt a, cpt b, cpt p, cpt q) {
	return deq(det(b-a,q-p),0);
}

//Distance between 2 finite lines
double distfinitelineline(cpt a,cpt b,cpt p,cpt q) {
	if (!arepara(a,b,p,q) && !std::isnan(intersectline(a,b,p,q).x)) return 0;
	
	/*if (arepara(a,b,p,q)) {
		b-=a;p-=a;q-=a;
		double sp=dot(b,p)/norm(b);
		if (0<sp && sp<1) return det(b,p)/abs(b);
	}*/
	
	return min({
			distfinitelinept(a,b,p),
			distfinitelinept(a,b,q),
			distfinitelinept(p,q,a),
			distfinitelinept(p,q,b)
	});
}

//This is kind of unnecessary
double distpara(Pt a, Pt b, Pt p, Pt q) {
	return distlinept(a,b,p);
}

struct Circle {
	Pt c;double r;
	bool operator==(const Circle& b) const {return c==b.c && deq(r,b.r);}
};

// Number of intersections, pair containing intersections
// 3 means infinitely many intersections. This also happens with identical
// radius 0 circles.
pair<int,pair<Pt,Pt>> intersectcirclecircle(const Circle& a,const Circle& b) {
	Pt v=b.c-a.c;
	//  disjoint           || one inside other
	if (a.r+b.r+EPS<abs(v) || abs(a.r-b.r)>abs(v)+EPS) return {0,{}};
	if (abs(v)<EPS) return {3,{}};
	double X=(norm(a.r)-norm(b.r)+norm(v))/(2.0*abs(v));
	double Ysq=norm(a.r)-norm(X),Y;
	v/=abs(v);
	if (Ysq<0 || (Y=sqrt(Ysq))<EPS) return {1,{Pt{X,0}*v+a.c,{}}};
	return {2,{Pt{X,Y}*v+a.c,Pt{X,-Y}*v+a.c}};
}

// Compute the area of the intersection of two circles
double area_of_intersection(Circle c1, Circle c2) {
  if (c1.r < c2.r) swap(c1,c2);  double d = abs(c1.c - c2.c);
  if (d + c2.r <= c1.r + EPS) return pi*c2.r*c2.r;  // c2 is in c1
  if (d >= c1.r + c2.r - EPS) return 0.0;       // no overlap
  double alpha = acos((c1.r*c1.r+d*d-c2.r*c2.r)/(2*c1.r*d))*2;
  double beta = acos((c2.r*c2.r+d*d-c1.r*c1.r)/(2*c2.r*d))*2;
  double a1 = 0.5*beta*c2.r*c2.r-0.5*c2.r*c2.r*sin(beta);
  double a2 = 0.5*alpha*c1.r*c1.r-0.5*c1.r*c1.r*sin(alpha);
  return a1 + a2;
}

//untested
pair<int,pair<Pt,Pt>> intersectfinitelinecircle(cpt a,cpt b,Circle c) {
	Pt v=b-a;
	v/=abs(v);
	c.c=(c.c-a)/v;
	if (c.r+EPS<abs(c.c.y)) return {0,{}};
	double offsq=norm(c.r)-norm(c.c.y),off;
	if (offsq<0 || (off=sqrt(offsq))<EPS) {
		if (-EPS<c.c.x && c.c.x<abs(v)+EPS) return {1,{Pt{c.c.x,0}*v+a,{}}};
	}
	// Section much shorter without bounds check
	pair<int,pair<Pt,Pt>> ans;
	for (int sgn=-1;sgn<2;sgn+=2) {
		double X=c.c.x+sgn*off;
		if (-EPS<X && X<abs(v)+EPS) { // line bounds check
			if (ans.X==0) ans.Y.X=Pt{X,0}*v+a;
			else ans.Y.Y=Pt{X,0}*v+a;
			ans.X++;
		}
	}
	return ans;
}

//untested
Circle circlefrom3points(cpt a,cpt b,cpt c) {
	Pt v=b-a;
	// Circle not unique or infinite if points are colinear
	if (deq(det(v,c-a),0)) return {Pt(),-1};
	double X=abs(v);
	v/=abs(v);
	Pt p=(c-a)/v;
	Pt q(X/2,(norm(p.x)-norm(p.y)-p.x*X)/(2*p.y));
	return {q*v+a,abs(q)};
}

template<class T,int maxn>
struct Arr {
	int n=0;
	T a[maxn]{};
};

// Up to 4 common tangents for two circles (except when infinitely many). Each
// tangent is two points in Arr.
// 
// These points represent where the tangent touches each circle. If these points
// are the same then the second point is to the right of the first (when looking
// from the center of the first circle), and the distance between the two points
// is the distance between the centers of the circles.
// 
// Outer tangents are before inner tangents since they occur whenever inner
// tangents do. The first tangent in each group is the one which intersects the
// first circle to the left of the second circle (when looking from the center
// of the first circle).
//
// The radii should be positive. 0 radii should work but give multiple identical
// lines.
Arr<Pt,8> commontangents(const Circle& a,const Circle& b) {
	Arr<Pt,8> ans;
	if (a==b) {ans.n=9; return ans;} // infinitely many
	Pt v=b.c-a.c;
	double X=abs(v);
	v/=norm(X);
	int &n=ans.n;
	//Pt *aa=ans.a;
	for (int sgn=-1;sgn<2;sgn+=2) {
		Pt u=a.r+sgn*b.r;
		if (X+EPS<abs(u.x)) break;
		u.y=norm(X)-norm(u.x);
		if (u.y>=0) u.y=sqrt(u.y);

		ans.a[n++]=a.r*u;
		ans.a[n++]=(a.r+(u.y<EPS?X:u.y)*Pt(0,-1))*u;
		if (u.y>=EPS) {
			ans.a[n++]=a.r*conj(u);
			ans.a[n++]=(a.r-u.y*Pt(0,-1))*conj(u);
		}
	}

	for (int i=0;i<n;i++) ans.a[i]=ans.a[i]*v+a.c;

	return ans;
}

double yval(cpt a,cpt b,double X) {
	return deq(a.x,b.x)?a.y:(a.y*(b.x-X)+b.y*(X-a.x))/(b.x-a.x);
}

// Closest pair of points. Complexity: O(N log(N))
pair<Pt,Pt> closest_pair(vector<Pt> P) {
  sort(P.begin(), P.end(), [](auto p1, auto p2) { return p1.y < p2.y; });
  set<Pt> a; double d = INF; pair<Pt,Pt> cp{{0,0},{INF,0}};
  for (auto p = P.begin(), c = p; c != P.end(); c++) {
    while (p != c && p->y < c->y-d) a.erase(*p++);
    for (auto i=a.lower_bound(Pt{c->x-d,0}); i != a.end() && i->x < c->x+d; i++)
      if (abs(*c - *i) < d) d = abs(*c - *i), cp = {*c, *i};
    a.insert(*c);
  }
  return cp;
}

// Returns true if the given point is contained within the given circle
bool point_in_circle(const Pt& p, const Circle& c) {
  return abs(p - c.c) <= c.r + EPS;
}

// Construct a circle from two antipodal points on the boundary
Circle circle_from_diameter(cpt a, cpt b) {
  return {0.5*(a+b), abs(0.5*(a+b) - a)};
}

// Find the smallest circle that encloses all of the given points. Complexity: O(N)
Circle minimum_enclosing_circle(vector<Pt> P) {
  int N = (int)P.size(); random_shuffle(P.begin(), P.end());
  Circle c{P[0], 0};
  for (int i=1; i<N; i++) if (!point_in_circle(P[i], c)) {
    c = Circle{P[i],0};
    for (int j=0; j<i; j++) if (!point_in_circle(P[j], c)) {
      c = circle_from_diameter(P[i],P[j]);
      for (int k=0; k<j; k++) if (!point_in_circle(P[k], c))
        c = circlefrom3points(P[i],P[j],P[k]);
    }
  }
  return c;
}

// Find the area of the union of the given circles. Complexity: O(n^2 log(n))
double circle_union_area(const vector<Circle>& cir) {
  int n = (int)cir.size(); vector<bool> ok(n, 1); double ans = 0.0;
  for (int i=0; i<n; i++) for (int j=0; j<n; j++) if (i != j && ok[j]) 
    if (abs(cir[i].c - cir[j].c)+cir[i].r-cir[j].r < EPS) { ok[i] = false; break; }
  for (int i=0; i<n; i++) if (ok[i]) {
    bool flag = false;  vector<pair<double,double>> reg;
    for (int j=0; j<n; j++) if (i != j && ok[j]) {
      auto p = intersectcirclecircle(cir[i], cir[j]);
      if (p.X < 2) continue; else flag = true;
      auto ang1 = arg(p.Y.Y - cir[i].c), ang2 = arg(p.Y.X - cir[i].c);
      if (ang1 < 0) ang1 += 2*pi;
      if (ang2 < 0) ang2 += 2*pi;
      if (ang1 > ang2) reg.emplace_back(ang1, 2*pi), reg.emplace_back(0, ang2);
      else reg.emplace_back(ang1, ang2);
    }
    if (!flag) { ans += pi*cir[i].r*cir[i].r; continue; }
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
      if (ang < 0) ang += 2*pi;
      ans += 0.5 * cir[i].r*cir[i].r * (ang - sin(ang));
    }
  }
  return ans;
}

// TODO: Convex cut

// TODO: Area of polygon and circle intersection

// I wrote this and then realised it was less code if both lower and upper hulls
// are computed simultaneously. Also this code is wrong. do convex hull to get
// correct answer.
/*
vector<Pt> intersectupperhulls(vector<Pt>& U1,vector<Pt>& U2) {
	vector<Pt> r;
	double x=max(U1[0].x,U2[0].x),y;
	int i=upper_bound(U1.begin(),U1.end(),Pt(x,DBL_MAX))-U1.begin()-1;
	int j=upper_bound(U2.begin(),U2.end(),Pt(x,DBL_MAX))-U2.begin()-1;
	while (1) {
		y=min(i+1==U1.size()?U1[i].y:yval(U1[i],U1[i+1],x),
			  j+1==U2.size()?U2[j].y:yval(U2[j],U2[j+1],x));
		r.emplace_back(x,y);
		if (i+1==U1.size() || j+1==U2.size()) break;
		if (!arepara(U1[i],U1[i+1],U2[j],U2[j+1])) {
			Pt p=intersectline(U1[i],U1[i+1],U2[j],U2[j+1]);
			if (!std::isnan(p.x)) r.push_back(p);
		}
		if (U1[i+1].x<U2[j+1].x) x=U1[++i].x;
		else x=U2[++j].x;
	}
	r.resize(unique(r.begin(),r.end())-r.begin());
	return r;
}
*/
// Split convex hull into lower and uppper hull. Endpoints included
// Untested
pair<vector<Pt>,vector<Pt>> splithull(const vector<Pt>& p) {
	int i=min_element(p.begin(),p.end())-p.begin();
	int j=max_element(p.begin(),p.end())-p.begin();
	int n=p.size();
	vector<Pt> L,U;
	for (int k=i;1;k++) {
		if (k==n) k=0;
		L.push_back(p[k]);
		if (k==j) break;
	}
	for (int k=i;1;k--) {
		if (k==-1) k=n-1;
		U.push_back(p[k]);
		if (k==j) break;
	}
	
	return {L,U};
}

// TODO: Lat-long conversion and great circle distance

// TODO: Diameter of a convex polygon

// TODO: Minkowski sum of convex hull
/*
Wiki:

The Minkowski sum of two sets of position vectors A and B in Euclidean space is formed

by adding each vector in A to each vector in B, i.e. the set

A + B = {⃗a +⃗b |⃗a ∈ A,⃗b ∈ B}.

For all subsets S1 and S2 of a real vector-space, the convex hull of their Minkowski sum

is the Minkowski sum of their convex hulls Conv(S1 + S2) = Conv(S1) + Conv(S2).

Minkowski sums are used in motion planning of an object among obstacles. They are

used for the computation of the configuration space, which is the set of all admissi-
ble positions of the object. In the simple model of translational motion of an object

in the plane, where the position of an object may be uniquely specified by the posi-
tion of a fixed point of this object, the configuration space are the Minkowski sum of

the set of obstacles and the movable object placed at the origin and rotated 180 degrees.

int minkowski(point_t *h, point_t *h1, point_t *h2, int n, int m) {

2 point_t c = point_t(0, 0);

3 for (int i = 1; i <= m; ++i) c = c + h2[i];

4 c = c / m;

5 for (int i = 1; i <= m; ++i) h2[i] = h2[i] - c;

6 int cur = -1;

7 for (int i = 1; i <= m; ++i) {

8 if (dblcmp(cross(h2[i], h1[1] - h1[n])) >= 0) {

9 if (cur == -1 || cross(h2[i], h1[1] - h1[n]) > cross(h2[cur], h1[1] - h1[n]))

cur = i;

10 }

11 }

12 int cnt = 0;

13 h1[n + 1] = h1[1];

14 for (int i = 1; i <= n; ++i) {

15 while (true) {

16 h[++cnt] = h1[i] + h2[cur];

17 int next = (cur == m ? 1 : cur + 1);

18 if (dblcmp(cross(h2[cur], h1[i + 1] - h1[i])) < 0) cur = next;

19 else {

20 if (cross(h2[next], h1[i + 1] - h1[i]) > cross(h2[cur], h1[i + 1] - h1[i]))

cur = next;

21 else break;

22 }

23 }

24 }

25 for (int i = 1; i <= cnt; ++i) h[i] = h[i] + c;

26 for (int i = 1; i <= m; ++i) h2[i] = h2[i] + c;

27 return graham(h, cnt);

28 }
*/

// TODO: Test me and make me correct
// Untested
vector<Pt> intersecthulls(vector<Pt>& p,vector<Pt>& q) {
	auto P=splithull(p),Q=splithull(q);
	vector<Pt> hulls[4]{P.X,Q.X,P.Y,Q.Y},A[2];
	double X=DBL_MIN,xmax=DBL_MAX,nx,yl,yu,poss;
	int i[4];
	for (int j=0;j<4;j++) {
		X=max(X,hulls[j][0].x);
		xmax=min(xmax,hulls[j].back().x);
		//i[j]=upper_bound(hulls[j].begin(),hulls[j].end(),Pt(x,DBL_MAX))
		//	-hulls[j].begin()-1;
	}
	while (1) {
		yl=DBL_MIN, yu=DBL_MAX;

		for (int j=0;j<4;j++)
			for (;i[j]+1<hulls[j].size() && hulls[j][i[j]].x<X;i[j]++);

		for (int j=0;j<2;j++) {
			yl=max(yl,i[j]+1==hulls[j].size()?hulls[j][i[j]].y:
					yval(hulls[j][i[j]],hulls[j][i[j]+1],X));
			yu=min(yu,i[j+2]+1==hulls[j+2].size()?hulls[j+2][i[j+2]].y:
					yval(hulls[j+2][i[j+2]],hulls[j+2][i[j+2]+1],X));
		}

		if (yu+EPS>yl) for (int j=0;j<2;j++) {
			Pt pt(X,j==0?yl:yu);
			while ((A[j].size() && A[j].back()==pt) || (A[j].size()>1 &&
					deq(det(A[j].back()-A[j][A[j].size()-2],pt-A[j].back()),0)))
				A[j].pop_back();
			A[j].push_back(pt);
		}

		if (X+EPS>xmax) break;

		nx=DBL_MAX;
		for (int j=0;j<4;j++) 
			if (/*i[j]+1<hulls[j].size() &&*/ (poss=hulls[j][i[j]+1].x)>=X+EPS)
				nx=min(nx,poss);
		for (int j=0;j<4;j+=2)
			if (!arepara(hulls[j][i[j]],hulls[j][i[j]+1],hulls[j+1][i[j+1]],hulls[j+1][i[j+1]+1])) {
				Pt pt=intersectline(hulls[j][i[j]],hulls[j][i[j]+1],hulls[j+1][i[j+1]],hulls[j+1][i[j+1]+1]);
				if (!std::isnan(pt.x) && pt.x>=X+EPS) nx=min(nx,pt.x);
			}
	}

	int l=A[0].size(),u=A[1].size();
	A[0].resize(l+u-2);
	copy(A[1].rend()-u+1,A[1].rend()-1,A[0].begin()+l);
	return A[0];
}




	/*for (Pt &pt:P.x) pt.y*=-1;
	for (Pt &pt:Q.x) pt.y*=-1;
	vector<Pt> l=intersectupperhulls(P.x,Q.x),u=intersectupperhulls(P.y,Q.y);
	for (Pt &pt:L) pt.y*=-1;
	vector<Pt> L,U;
	double x=max(l[0].x,u[0].x),yl,yu;
	while (1)*/


namespace Test {
	void prtcommontangents(Circle a,Circle b) {
		Arr<Pt,8> ans=commontangents(a,b);
		for (int i=0;i<ans.n;i++) cout << ans.a[i] << endl;
	}

	void testcommontangents(Circle a,Circle b,vector<Pt> expected) {
		Arr<Pt,8> ans=commontangents(a,b);
		assert(ans.n==expected.size());
		if (ans.n<=8) for (int i=0;i<ans.n;i++) assert(ans.a[i]==expected[i]);
	}

	void solve() {
		Pt a{3,2};
		Pt b(2,-7);
		vector<Pt> aa{{1,2},{3,4},{5,6}};
		//for (Pt p:convexhull(aa)) cout << p << endl;
		//cout << endl;


		//cout << a + b << endl;
		//cout << a - b << endl;

		Pt c(0,0);
		Pt d(1e-10, 1e-10);

		//cout << (a == b) << endl;
		//cout << (c == d) << endl;

		cout << "dot det\n";
		//cout << dot(a,b) << endl;
		//cout << det(a,b) << endl;
		//cout << rotate({1,0},pi/2) << endl;


		a.x=0;
		//cout << a.x << endl;

		//cout << ((double*)a)[0] << endl;
		//a[0]=50;
		//cout << a[0] << endl;

		//cout << intersectline({10,1},{11,2},{11,1},{10,2}) << endl;
		cout << intersectline({10,1},{9,0},{11,1},{10,2}) << endl;

		cout << ptonseg({3,2},{5,4},{4,3}) << endl;
		cout << ptonseg({3,2},{5,4},{3,4}) << endl;
		cout << ptonseg({3,2},{5,4},{5,2}) << endl;

		vector<Pt> poly{{0,0},{2,2},{0,1}};
		cout << ptinpoly(poly,{1,1}) << endl;
		cout << ptinpoly(poly,{2,2}) << endl;
		cout << ptinpoly(poly,{-1,2}) << endl;
		cout << ptinpoly(poly,{4,3}) << endl;
		cout << ptinpoly(poly,{0.5,.75}) << endl;


		cout << fixed << setprecision(15);

		Arr<Pt,8> ans;
		testcommontangents(
				{{1,1},2},{{5,1},2},
				{{1,3},{5,3},{1,-1},{5,-1},{3,1},{3,-3}}
				);
		testcommontangents(
				{{-3,3},2},{{-3,8},2},
				{{-5,3},{-5,8},{-1,3},{-1,8},
				{-4.2,4.6},{-1.8,6.4},{-1.8,4.6},{-4.2,6.4}}
				);
		testcommontangents(
				{{-18,2},3.5},{{-13,-10},8.5},
				{{-15.535502958579881,4.485207100591716},
				{-7.014792899408285,-3.964497041420117},
				{-21.500000000000000,2.000000000000000},
				{-21.500000000000000,-10.000000000000000},
				{-15.514792899408285,-0.464497041420119},
				{-19.035502958579883,-4.014792899408285},
				{-18.000000000000000,-1.500000000000000},
				{-13.000000000000000,-1.500000000000000}}
				);
		testcommontangents(
				{{4,5},3},{{2,3},2},
				{{5.234313483298443,2.265686516701557},
				{2.822875655532295,1.177124344467705},
				{1.265686516701557,6.234313483298443},
				{0.177124344467705,3.822875655532295}}
				);
		testcommontangents(
				{{1,-5},1},{{1,-10},6},
				{{1.000000000000000,-4.000000000000000},
				{6.000000000000000,-4.000000000000000}}
				);
		testcommontangents(
				{{1,-5},0.5},{{1,-10},6},
				{}
				);
		testcommontangents(
				{{34.2,-57.88},39.9},{{34.2,-57.88},39.9},
				vector<Pt>(9)
				);
	}
}


namespace UVA_920 {
	void solve() {
		int T,n;
		cin >> T;
		cout << fixed << setprecision(2);
		while (T--) {
			cin >> n;
			vector<Pt> pts(n);
			for (Pt &p:pts) cin >> p;
			sort(pts.rbegin(),pts.rend());
			double ans=0,M=0;
			for (int i=1;i<n;i++) {
				Pt cept=intersectline(pts[i-1],pts[i],{-1,M},{30001,M});
				if (pts[i].y>M) ans+=abs(cept-pts[i]);
				M=max(M,pts[i].y);
			}

			cout << ans << '\n';
		}
	}
}

namespace UVA_10263 {
	Pt distfinitelinept(Pt a, Pt b, Pt p) {
		b-=a;p-=a;
		double sp=(p/b).x;//dot(b,p)/norm(b);
		Pt closest;
		if (sp>=0) {
			if (sp>1) closest=b;
			else closest=sp*b;
		}
		return closest+a;
	}
	void solve() {
		cout << fixed << setprecision(4);
		Pt M;
		int n;
		while (cin >> M) {
			Pt ans(1e100,1e100);
			cin >> n;
			vector<Pt> pts(n+1);
			for (Pt &p:pts) cin >> p;
			for (int i=1;i<=n;i++) {
				Pt cur=UVA_10263::distfinitelinept(pts[i-1],pts[i],M);
				if (abs(cur-M)<abs(ans-M)) ans=cur;
			}
			cout << ans.x << '\n' << ans.y << '\n';
		}
	}
}

namespace UVA_10927 {
	double safearg(cpt a) {
		double ag=arg(a);
		if (deq(ag,-pi)) return pi;
		return ag;
	}
	bool less(cpt a,cpt b) {
		double ag=safearg(a),bg=safearg(b);
		return ag<bg-EPS || (deq(ag,bg) && abs(a)<abs(b));
	}

	typedef pair<Pt,int> Pole;

	void solve() {
		int n;
		int cas=1;
		cout << fixed  << setprecision(0);
		while (cin >> n && n) {
			vector<Pole> poles(n);
			for (Pole &pol:poles) {
				cin >> pol.X >> pol.Y;
			}
			sort(poles.begin(),poles.end(), [](const Pole& a,const Pole& b) {return less(a.X,b.X);});

			vector<Pt> ans;

			for (int i=0,j=0;i<n;i=j) {
				for (;j<n && deq(safearg(poles[i].X),safearg(poles[j].X));j++);
				int M=-1;
				for (int k=i;k<j;k++) {
					if (poles[k].Y<=M) ans.push_back(poles[k].X);
					M=max(M,poles[k].Y);
				}
			}

			cout << "Data set " << cas++ << ":\n";
			if (ans.empty()) {
				cout << "All the lights are visible.\n";
			}
			else {
				sort(ans.begin(),ans.end());
				cout << "Some lights are not visible:\n";
				for (int i=0;i<ans.size();i++) {
					cout << "x = " << ans[i].x << ", y = " << ans[i].y;
					if (i+1==ans.size()) cout << ".\n";
					else cout << ";\n";
				}
			}
		}
	}

}

namespace UVA_378 {
	Pt intersectline(Pt a, Pt b, Pt p, Pt q) {
		Pt ab=b-a,qp=p-q,ap=p-a;
		double s=det(ap,qp)/det(ab,qp),t=det(ab,ap)/det(ab,qp);
		//double s,t;
		//tie(s,t)=solve(b-a,p-q,p-a);
		//a+t(b-a)=p+s(q-p)

		//Can also just use ptonseg.
		//if (-EPS<s && s<1+EPS //Answer is on ab
		//	&& -EPS<t && t<1+EPS) //Answer is on pq 
		return a+s*(b-a);
		return Pt(NAN,NAN);
	}

	void solve() {
		int N;
		cin >> N;
		cout << "INTERSECTING LINES OUTPUT\n" << fixed << setprecision(2);
		while (N--) {
			Pt a,b,c,d;
			cin >> a >> b >> c >> d;
			if (arepara(a,b,c,d)) {
				if (distlinept(a,b,c)>EPS) cout << "NONE\n";
				else cout << "LINE\n";
			}
			else {
				Pt ans=UVA_378::intersectline(a,b,c,d);
				cout << "POINT " << ans.x << ' ' << ans.y << '\n';
			}
		}
		cout << "END OF OUTPUT\n";
	}
}

namespace UVA_191 {
	void solve() {
		int n;
		cin >> n;
		while (n--) {
			Pt a,b,c,d;
			cin >> a >> b >> c >> d;
			if (c.x>d.x) swap(c.x,d.x);
			if (c.y>d.y) swap(c.y,d.y);
			if ((c.x-EPS<a.x && a.x<d.x+EPS && c.y-EPS<a.y && a.y<d.y+EPS)
					|| !std::isnan(intersectline(a,b,{c.x,c.y},{c.x,d.y}).x)
					|| !std::isnan(intersectline(a,b,{c.x,c.y},{d.x,c.y}).x)
					|| !std::isnan(intersectline(a,b,{c.x,d.y},{d.x,d.y}).x)
					|| !std::isnan(intersectline(a,b,{d.x,c.y},{d.x,d.y}).x))
				cout << "T\n";
			else cout << "F\n";
		}
	}
}

namespace AIZU_CGL_3_C {
	void solve() {
		int n;
		cin >> n;
		vector<Pt> poly(n);
		for (Pt &p:poly) cin >> p;
		int q;
		cin >> q;
		while (q--) {
			Pt p;
			cin >> p;
			cout << ptinpoly(poly,p) << '\n';
		}
	}
}


namespace UVA_634 {
	void solve() {
		int n;
		while (cin >> n && n) {
			vector<Pt> poly(n);
			for (Pt &p:poly) cin >> p;
			Pt q;
			cin >> q;
			if (ptinpoly(poly,q)) cout << "T\n";
			else cout << "F\n";
		}
	}
}

namespace UVA_453 {
	void prt(Pt a) {
		stringstream ss;
		ss << fixed << setprecision(3);
		ss << a.x << ' ' << a.y;
		string s1,s2;ss >> s1 >> s2;
		if (s1=="-0.000") s1="0.000";
		if (s2=="-0.000") s2="0.000";
		cout << "(" << s1 << "," << s2 << ")";
	}

	void solve() {
		cout << fixed << setprecision(3);
		Circle a,b;
		while (cin >> a.c >> a.r >> b.c >> b.r) {
			pair<int,pair<Pt,Pt>> ans=intersectcirclecircle(a,b);
			if (ans.X==0) cout << "NO INTERSECTION\n";
			else if (ans.X==3) {
				if (a.r>EPS) cout << "THE CIRCLES ARE THE SAME\n";
				else cout << a.c << '\n';
			}
			else {
				if (ans.X==1) prt(ans.Y.X);
				else {
					prt(min(ans.Y.Y,ans.Y.X));
					prt(max(ans.Y.Y,ans.Y.X));
					/*string s1,s2;
					ss >> s1 >> s2;
					cout << s1;
					if (s2!=s1) cout << s2;
					//cerr << s1 << s2 << ' ' << (s1!=s2) << ' ' << (s1==s2) << ' ' << (s2!=s1) << endl;*/
				}
				cout << '\n';
			}
		}
	}
}

int main() {
	
	Test::solve();
	//UVA_920::solve();
	//UVA_10263::solve();
	//UVA_10927::solve(); // Use EPS=1e-10
	//UVA_378::solve();
	//UVA_191::solve();
	//AIZU_CGL_3_C::solve();
	//UVA_634::solve();
	//UVA_453::solve(); // Use EPS=1e-4
}
