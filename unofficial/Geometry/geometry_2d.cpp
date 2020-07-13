#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef pair<int,int> pii;

//listings:geometry
// --------------------------- 2D Computational Geometry -------------------------------
#define x real()
#define y imag()
#define cpt const Pt&

const double EPS = 1e-9;
const double pi=acos(-1);
const double inf=1e100;
bool deq(double a,double b) {return abs(a-b)<EPS;}

typedef complex<double> cpx;
struct Pt : cpx {
  Pt() = default;	using cpx::cpx;
	Pt(cpx a) : cpx(a) {}
  double& x const { return (double&)*this; }
  double& y const { return ((double*)this)[1]; }
	bool operator ==(cpt b) const {return abs(*this-b) < EPS; }
	bool operator <(cpt b) const {return x<b.x || (x==b.x && y<b.y); }
};

bool epsless(cpt a,cpt b) {return a.x+EPS<b.x || (deq(a.x,b.x) && a.y<b.y);}

double dot(cpt a, cpt b) {return (conj(a) * b).x;} // Dot product
double det(cpt a, cpt b) {return (conj(a) * b).y;} // Determinant/"Cross Product"
double angle(cpt a, cpt b) {return arg(b - a);} // [-pi,pi] a to b with x axis
double angle (cpt a, cpt b, cpt c) {return arg((a-b)/(c-b));} //[-pi,pi]
double slope(cpt a, cpt b) {return (b.y-a.y)/(b.x-a.x);}

Pt rotate(cpt a, double theta) {return a * polar((double)1.0, theta);}
// rotate a around p by theta anticlockwise
Pt rotate(cpt a, cpt p, double theta) {return rotate(a - p,theta) + p;}
Pt project(cpt p, cpt v) {return v * dot(p, v) / norm(v);} // p onto v
Pt project(cpt p, cpt a, cpt b) {return a+project(p-a,b-a);} // p onto line (a,b)
// reflect p across the line (a,b)
Pt reflect(cpt p, cpt a, cpt b) {return a + conj((p - a) / (b - a)) * (b - a);}

bool collinear(Pt a, Pt b, Pt c) {return deq(det(b-a,c-b),0);}
bool areperp(cpt a,cpt b,cpt p,cpt q) {	return deq(dot(b-a,q-p),0); }
bool arepara(cpt a, cpt b, cpt p, cpt q) { return deq(det(b-a,q-p),0); }

// Orientation test (1 anticlockwise, -1 clockwise, 0 collinear)
int orient(cpt a, cpt b, cpt c) {
	double d=det(b-a,c-b);
	return d>EPS?1:d<-EPS?-1:0;
}

//Compare points by principal argument (-pi,pi] breaking ties by norm.
//0 is considered less than everything else.
bool argcomp(cpt a,cpt b) {
	if (b==0) return 0;
	if (a==0) return 1;
	double a1=arg(a),a2=arg(b);
	if (a1<-pi+EPS/2) a1+=2*pi;
	if (a2<-pi+EPS/2) a2+=2*pi;
	return a1+EPS<a2 || (deq(a1,a2) && norm(a)<norm(b));
}

// Point on line segment (including endpoints)
bool ptonseg(cpt a, cpt b, cpt p) {
	Pt u=b-a,v=p-a;
	return a==p || b==p || ((0<dot(u,v) && dot(u,v)<norm(u)) && deq(det(u,v),0));
}

// Signed area of polygon. Positive for anticlockwise orientation.
double polygonarea(const vector<Pt>& p) {
	double r=0;	int n=p.size();
	for (int j=0,i=n-1;j<n;i=j++) r+=det(p[i],p[j]);
	return r/2;
}

// Convex hull O(NlogN). Be careful of duplicate or very close points.
// if all points are colinear the middle points come up twice forwards and
// backwards e.g. a-b-c-d becomes a-b-c-d-c-b
// To remove colinear points change <-EPS and >EPS to <EPS and >-EPS.
vector<Pt> convexhull(vector<Pt> p) {
  sort(p.begin(),p.end(),epsless); p.resize(unique(p.begin(),p.end())-p.begin());
  int l=0,u=0; vector<Pt> L(p),U(p);
  if (p.size()<=2) return p;
  for (Pt& i:p) {
    while (l>1 && det(i-L[l-1],L[l-2]-i)<-EPS) l--;
    while (u>1 && det(i-U[u-1],U[u-2]-i)>EPS) u--;
    L[l++]=U[u++]=i;
  }
  L.resize(l+u-2); copy(U.rend()-u+1,U.rend()-1,L.begin()+l);
  return L;
}

// Point in polygon test O(N)
// Returns: 0 if not in polygon, 1 if on boundary, 2 if in interior
int ptinpoly(const vector<Pt>& p, cpt q) {
	int n=p.size(), i,j,r=0;
	for (j=0,i=n-1;j<n;i=j++) {
		if (ptonseg(p[i],p[j],q)) return 1;
		if (((p[i].y <= q.y && q.y < p[j].y) || (p[j].y <= q.y && q.y < p[i].y))
			&& q.x < (p[j].x-p[i].x) * (q.y-p[i].y)/(p[j].y-p[i].y) + p[i].x) r^=2;
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

Pt solve(cpt a, cpt b, cpt v) { // solves [a b]x==v with Cramer's rule.
	return Pt(det(v,b)/det(a,b),det(a,v)/det(a,b));
}

//Intersection of 2 line segments. Divides by 0 if they are parallel.
//Returns {nan,nan} if they don't intersect.
//Remove if statements below to get infinite lines.
Pt intersectline(Pt a, Pt b, Pt p, Pt q) {
	Pt ab=b-a,qp=p-q,ap=p-a;
	double s=det(ap,qp)/det(ab,qp),t=det(ab,ap)/det(ab,qp);
	if (-EPS<s && s<1+EPS //Answer is on ab
		&& -EPS<t && t<1+EPS) //Answer is on pq 
      return a+s*ab;
	return Pt(NAN,NAN);
}

Pt intersectlineexact(Pt a, Pt b, Pt p, Pt q) {
	Pt ab=b-a,qp=p-q,ap=p-a;
	double s=det(ap,qp)/det(ab,qp),t=det(ab,ap)/det(ab,qp);
	if (0<s && s<1 //Answer is on ab
		&& 0<t && t<1) //Answer is on pq 
      return a+s*ab;
	return Pt(NAN,NAN);
}

//Distance between infinite line and point.
double distlinept(cpt a, cpt b, cpt p) { return abs(det(b-a,p-a)/abs(b-a)); }

//Distance between finite line and point
double distfinitelinept(Pt a, Pt b, Pt p) {
	b-=a;p-=a;  Pt closest;  double sp=(p/b).x; //dot(b,p)/norm(b);
	if (sp>=0) {
		if (sp>1) closest=b;
		else closest=sp*b;
	}
	return abs(closest-p); // Note that actual closest Pt on line is closest + a
}

//Distance between 2 finite lines
double distfinitelineline(cpt a,cpt b,cpt p,cpt q) {
	if (!arepara(a,b,p,q) && !std::isnan(intersectlineexact(a,b,p,q).x)) return 0;
	return min({distfinitelinept(a,b,p),distfinitelinept(a,b,q),
	  distfinitelinept(p,q,a),distfinitelinept(p,q,b)});
}

struct Circle {
	Pt c;double r;
	bool operator==(const Circle& b) const {return c==b.c && deq(r,b.r);}
};

// Number of intersections, pair containing intersections
// 3 means infinitely many intersections. This also happens with identical
// radius 0 circles.
pair<int,pair<Pt,Pt>> intersect(const Circle& a,const Circle& b) {
	Pt v=b.c-a.c; // disjoint || one inside other
	if (a.r+b.r+EPS<abs(v)    || abs(a.r-b.r)>abs(v)+EPS) return {0,{}};
	if (abs(v)<EPS) return {3,{}};
	double X=(norm(a.r)-norm(b.r)+norm(v))/(2.0*abs(v)), Ysq=norm(a.r)-norm(X),Y;
	v/=abs(v);
	if (Ysq<0 || (Y=sqrt(Ysq))<EPS) return {1,{Pt{X,0}*v+a.c,{}}};
	return {2,{Pt{X,Y}*v+a.c,Pt{X,-Y}*v+a.c}};
}

pair<int,pair<Pt,Pt>> intersectfinitelinecircle(cpt a,cpt b,Circle c) {
	Pt v=b-a;	v/=abs(v);	c.c=(c.c-a)/v;
	if (c.r+EPS<abs(c.c.y)) return {0,{}};
	double offsq=norm(c.r)-norm(c.c.y),off;
	if (offsq<0 || (off=sqrt(offsq))<EPS)
  if (-EPS<c.c.x && c.c.x<abs(v)+EPS) return {1,{Pt{c.c.x,0}*v+a,{}}};
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

Circle circlefrom3points(cpt a,cpt b,cpt c) {
	Pt v=b-a; double X=abs(v); v/=X; Pt p=(c-a)/v;
	if (deq(det(v,c-a),0)) return {0,-1}; // Not unique or infinite if collinear
	Pt q(X/2,(norm(p.x)-norm(p.y)-p.x*X)/(2*p.y));
	return {q*v+a,abs(q)};
}

// Peter's custom array
template<class T,int maxn> struct Arr {	int n=0; T a[maxn]={}; };

// Each tangent is two points in Arr.
// These points represent where the tangent touches each circle. If these points
// are the same then the second point is to the right of the first (when looking
// from the center of the first circle), and the distance between the two points
// is the distance between the centers of the circles.
// Outer tangents are before inner tangents since they occur whenever inner
// tangents do. The first tangent in each group is the one which intersects the
// first circle to the left of the second circle (when looking from the center
// of the first circle).
// The radii should be positive. 0 radii should work but give multiple identical lines.
Arr<Pt,8> commontangents(const Circle& a,const Circle& b) {
	Arr<Pt,8> ans; Pt v=b.c-a.c; double X=abs(v); v/=norm(X); int &n=ans.n;
	if (a==b) {ans.n=9; return ans;} // infinitely many
	for (int sgn=-1;sgn<2;sgn+=2) {
		Pt u=a.r+sgn*b.r;
		if (X+EPS<abs(u.x)) break;
		u.y=norm(X)-norm(u.x), u.y=u.y>0?sqrt(u.y):0;
		ans.a[n++]=a.r*u, ans.a[n++]=(a.r+(u.y<EPS?X:u.y)*Pt(0,-1))*u;
		if (u.y>=EPS) ans.a[n++]=a.r*conj(u), ans.a[n++]=(a.r-u.y*Pt(0,-1))*conj(u);
	}
	for (int i=0;i<n;i++) ans.a[i]=ans.a[i]*v+a.c;
	return ans;
}

// Find the max dot product of a point in p with v. p must be a convex polygon
// where no three points are collinear and dot(p[l],v)<dot(p[l+1],v). O(log n)
int maxdot(int l,int r,const vector<Pt>& p,cpt v) {
	if (r-l<10) {
		int i=l;
		for (int j=l+1;j<r;j++) if (dot(p[i],v)<dot(p[j],v)) i=j;
		return i;
	}
	int m1=(2*l+r)/3,m2=(l+2*r)/3; double d1=dot(p[m1],v),d2=dot(p[m2],v);
	if (d1<dot(p[l],v)) return maxdot(l,m1,p,v);
	if (d1>d2) return maxdot(l,m2,p,v);
	return maxdot(m1+1,r,p,v);
}

// Min and max dot product of a point in p with v. p must be a convex polygon
// where no three points are collinear. Indices are returned. O(log n)
pii minmaxdot(const vector<Pt>& p,cpt v) {
	int i=deq(dot(p[0],v),dot(p[1],v)), n=p.size(),m,M;
	if (dot(p[i],v)<dot(p[i+1],v)) M=maxdot(i,n,p,v),m=maxdot(M,n,p,-v);
	else m=maxdot(i,n,p,-v),M=maxdot(m,n,p,v);
	for (int j=0;j<=i;j++) {
		if (dot(p[j],v)>dot(p[M],v)) M=j;
		if (dot(p[j],v)<dot(p[m],v)) m=j;
	}
	return {m,M};
}

//Returns convex hull of all points x within the convex polygon p, which satisfy
//det(b-a,x-a)>=0. Returned polygon may be degenerate if the cut runs across an
//edge. For p ordered counterclockwise, the cut polygon is on the left of a->b
vector<Pt> convexcut(cpt a,cpt b,const vector<Pt>& p) {
	int n=p.size();	vector<Pt> r;
	for (int i=n-1,j=0;j<n;i=j++) {
		double d1=det(b-a,p[i]-a),d2=det(b-a,p[j]-a);
		if (d1>-EPS) r.push_back(p[i]);
		if ((d1>EPS && d2<-EPS) || (d1<-EPS && d2>EPS))
			r.push_back(intersectline(a,b,p[i],p[j])); //infinite lines
	}
	return r;
}

// Facilitates queries for the pair of furthest visible points on a convex polygon from
// some external point. P must be a non-degenerate convex polygon. Returns an empty
// interval for p inside P. Complexity: O(N log(N)) pre-process, O(log(N)) per query.
struct PolygonTangents {
  vector<Pt> P; Pt c; int n; vi ids;
  PolygonTangents(vector<Pt> poly) : P(move(poly)), n(P.size()), ids(2*n) {
    for (auto p : P) c += 1.0/n * p;
    for (auto& p : P) p -= c;  sort(P.begin(), P.end(), argcomp);
    iota(ids.begin(),ids.begin()+n,0), iota(ids.begin()+n,ids.end(),0);
  }
  pii query(Pt p) { // Returns {i,j} such that P[i..j] are visible to p
    int a = lower_bound(P.begin(),P.end(),p-c,argcomp)-P.begin();
    int b = lower_bound(P.begin(),P.end(),-(p-c),argcomp)-P.begin();
    if (b < a) b += n;
    auto seen = [&](int i) { return orient(P[i?i-1:n-1],P[i],p-c) < 0; };
    int r = *partition_point(ids.begin()+a,ids.begin()+b,seen);
    int l = *partition_point(ids.rbegin()+n-a,ids.rbegin()+2*n-b,seen);
    return {l,r-1+(r?0:n)};
  }
  Pt operator[](int i) { return P[i]+c; }  // Get the (untranslated) i'th point
};

//Signed Area of polygon and circle intersection. Sign is determined by
//orientation of polygon. Divides by 0 if adjacent points are identical.
double areapolygoncircle(vector<Pt> p,Circle c) {
	int n=p.size();	double r=0;
	for (int i=n-1,j=0;j<n;i=j++) {
		Pt v=abs(p[j]-p[i])/(p[j]-p[i]), a=(p[i]-c.c)*v,b=(p[j]-c.c)*v;
		if (deq(a.y,0)) continue;
		double d=sqrt(max(0.0, norm(c.r)-norm(a.y)));
		r+=norm(c.r)*(atan2(b.y,min(b.x,-d))-atan2(a.y,min(a.x,-d))
				+atan2(b.y,max(b.x,d))-atan2(a.y,max(a.x,d)))
			+a.y*(min(max(a.x,-d),d)-min(max(b.x,-d),d));
	}
	return r/2;
}

// Closest pair of points. Complexity: O(N log(N))
pair<Pt,Pt> closest_pair(vector<Pt> P) {
  sort(P.begin(), P.end(), [](auto p1, auto p2) { return p1.y < p2.y; });
  set<Pt> a; double d = inf; pair<Pt,Pt> cp{{0,0},{inf,0}};
  for (auto p = P.begin(), c = p; c != P.end(); c++) {
    while (p != c && p->y < c->y-d) a.erase(*p++);
    for (auto i=a.lower_bound(Pt{c->x-d,0}); i != a.end() && i->x < c->x+d; i++)
      if (abs(*c - *i) < d) d = abs(*c - *i), cp = {*c, *i};
    a.insert(*c);
  }
  return cp;
}

//Diameter of convex polygon. Complexity: O(N)
double polygondiameter(const vector<Pt>& p) {
	int i=min_element(p.begin(),p.end())-p.begin(),ic=0,n=p.size(),ni=(i+1)%n;
	int j=max_element(p.begin(),p.end())-p.begin(),jc=0,nj=(j+1)%n;
	double r=0;
	while (ic<n || jc<n) {
		r=max(r,abs(p[j]-p[i]));
		if (det(p[ni]-p[i],p[j]-p[nj])>0) {
			i=ni++;ic++;
			if (ni==n) ni=0;
		}
		else {
			j=nj++;jc++;
			if (nj==n) nj=0;
		}
	}
	return r;
}

//Minimum width of a bounding rectangle of a convex polygon O(n)
//The polygon must have positive signed area.
double minboundingwidth(const vector<Pt>& p) {
	double r=DBL_MAX;	int n=p.size();
	for (int i=n-1,j=0,k=0,nk;j<n;i=j++) {
		Pt v=p[j]-p[i];v/=abs(v);
		for (;det(v,p[nk=k+1==n?0:k+1]-p[i])>det(v,p[k]-p[i]);k=nk);
		r=min(r,det(v,p[k]-p[i]));
	}
	return r;
}

//Minkowski sum of convex polygons O(n)
//Polygon is returned with the minimum number of points. i.e. No three points
//will be collinear. The input polygons must have positive signed area.
vector<Pt> minkowskisum(const vector<Pt>& p,const vector<Pt>& q) {
	vector<Pt> r;  int n=p.size(),m=q.size();
	int i=min_element(p.begin(),p.end(),epsless)-p.begin(),oi=i,ni=(i+1)%n;
	int j=min_element(q.begin(),q.end(),epsless)-q.begin(),oj=j,nj=(j+1)%m;
	do {
		r.push_back(p[i]+q[j]);
		Pt v=det(p[ni]-p[i],q[nj]-q[j])>0?p[ni]-p[i]:q[nj]-q[j];
		while (det(v,p[ni]-p[i])<EPS) {
			i=ni++;
			if (ni==n) ni=0;
		}
		while (det(v,q[nj]-q[j])<EPS) {
			j=nj++;
			if (nj==m) nj=0;
		}
	} while (i!=oi || j!=oj);
	return r;
}

// Returns true if the given point is contained within the given circle
bool point_in_circle(const Pt& p, const Circle& c) { return abs(p - c.c) <= c.r + EPS; }
// Construct a circle from two antipodal points on the boundary
Circle circle_from_diameter(cpt a, cpt b) { return {0.5*(a+b), abs(0.5*(a+b) - a)}; }

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
      auto p = intersect(cir[i], cir[j]);
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

// Find a pair of intersecting lines. Complexity: O(N log(N))
#define cl const Line&
struct Line {
	Pt u,v; //Endpoints
	double m() const {return (u.y-v.y)/(u.x-v.x);}
	double c() const {return yv(0);}
	double yv(double X) const {return det(u-X,v-X)/(u.x-v.x);}
	bool operator<(cl b) const {return u<b.u || (!(b.u<u) && v<b.v);}
};
namespace FindIntersection {
	typedef pair<double,int> pdi;
	const int maxn=300000;	Line segs[maxn];	pdi ord[2*maxn];
	int sgndiff(double a,double b) {return (a+EPS<b)-(b+EPS<a);}
	bool comp(int i,int j) {
		cl a=segs[i],b=segs[j];
		int by,bg;
		if (deq(a.u.x,b.u.x)) by=sgndiff(a.u.y,b.u.y);
		else if (a.u.x<b.u.x) by=sgndiff(0,det(a.v-a.u,b.u-a.u));
		else by=sgndiff(det(b.v-b.u,a.u-b.u),0);
		bg=sgndiff(0,det(a.v-a.u,b.v-b.u));
		return by==1 || (by==0 && (bg==1 || (bg==0 && i<j)));
	}
	set<int,bool(*)(int,int)> L(comp);
	pii checkpair(int i,int j) {
		cl a=segs[i],b=segs[j];	Pt ab=a.v-a.u,qp=b.u-b.v,ap=b.u-a.u;
		double d1=det(ap,qp),d2=det(ab,ap),d3=det(ab,qp);
		if (d3<0) d1*=-1,d2*=-1,d3*=-1;
		if (deq(d3,0)) {// Parallel
			Pt v=ab/abs(ab);  double c=dot(v,b.u),d=dot(v,b.v);
			if (d<c) swap(c,d);
			return {max(c,dot(v,a.u))+EPS<min(d,dot(v,a.v)) && deq(d1,0)?i:-1,j};
		}
		if (-EPS<d1 && d1<d3+EPS && -EPS<d2 && d2<d3+EPS) {
			if (EPS<d1 && d1+EPS<d3) return {i,j};
			if (EPS<d2 && d2+EPS<d3) return {j,i};
		}
		return {-1,0};
	}
	// Returns a pair of indices such that the first segment's interior
	// intersects with the other segment. First item is -1 if there are no such
	// segments.
	pii findintersection(const vector<Line>& lines) {
		int n=lines.size();copy(lines.begin(),lines.end(),segs);L.clear(); pii r;
		for (int i=0;i<n;i++) {
			if (epsless(segs[i].v,segs[i].u)) swap(segs[i].u,segs[i].v);
			ord[2*i]={segs[i].u.x,i};
			ord[2*i+1]={segs[i].v.x,i+n};
		}
		sort(ord,ord+2*n,[](const pdi& a,const pdi& b) {
      return a.X+EPS<b.X || (deq(a.X,b.X) && a.Y<b.Y); });
		for (int i=0;i<2*n;i++) {
			int j=ord[i].Y;
			if (j<n) {
				auto oit=L.insert(j).X,it=oit++;
				if (oit!=L.end() && (r=checkpair(*it,*oit)).X!=-1) return r;
				if (it!=L.begin() && (r=checkpair(*prev(it),*it)).X!=-1) return r;
			}
			else {
				auto it=L.erase(L.find(j-n));
				if (it!=L.begin() && it!=L.end()
					&& (r=checkpair(*prev(it),*it)).X!=-1) return r;
			}
		}
		return {-1,0};
	}
}

// Split convex hull into lower and uppper hull. Endpoints included
pair<vector<Pt>,vector<Pt>> splithull(vector<Pt> p) {
	rotate(p.begin(),min_element(p.begin(),p.end()),p.end());
	auto it=max_element(p.begin(),p.end());
	vector<Pt> L(p.begin(),it+1),U(it,p.end());
	U.push_back(p[0]), reverse(U.begin(),U.end());
	return {L,U};
}

//Intersect convex polygons O(n)
//Run convex hull to remove collinear points if required
//Beware of very steep but not vertical lines when polygon coordinates can
//differ by less than EPS without being equal. Undef defs if you want.
vector<Pt> intersectpolygons(const vector<Pt>&P,const vector<Pt>& Q) {
#define u(j) h[j][i[j]]
#define v(j) h[j][i[j]+1]
#define b(j) i[j]+1<h[j].size()
#define loop for (int j=0;j<4;j++)
#define yv(j) b(j)?det(u(j)-X,v(j)-X)/(u(j).x-v(j).x):u(j).x
	auto c=splithull(P),d=splithull(Q);	vector<Pt> h[]{c.X,d.X,c.Y,d.Y},L,U;
	int i[4]{};	double l=-inf,r=inf,X=-inf,nX;
	loop { r=min(r,h[j].back().x), l=max(l,u(j).x); }
	while (1) {
		nX=inf;
		loop if (b(j) && v(j).x>X+EPS) nX=min(nX,v(j).x);
		loop for (int k=j+1;k<4;k++) if (b(j) && b(k)) {
			double p=intersectline(u(j),v(j),u(k),v(k)).x;
			if (!std::isnan(p) && p>X+EPS) nX=min(nX,p);
		}
		if ((X=max(nX,l))>r+EPS) break;
		loop while (b(j) && v(j).x<X+(1-2*deq(X,r))*EPS) i[j]++;
		double m=max(yv(0),yv(1)),M=min(yv(2),yv(3));
		if (m<M+EPS) {
			L.emplace_back(X,m);
			if (!deq(m,M)) U.emplace_back(X,M);
		}
	}
	L.insert(L.end(),U.rbegin(),U.rend());
	return L;
}
//listings:/geometry

int main() {

}








