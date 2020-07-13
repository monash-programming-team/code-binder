#include<bits/stdc++.h>

using namespace std;
//using namespace rel_ops;

typedef long long ll;
typedef pair<int,int> ii;

// Can't use x and y for first and second
#define x real()
#define y imag()
#define X first
#define Y second
#define cpt const Pt&

//typedef long double ld;
//#define double ld
const double EPS = 1e-9;
const double pi=acos(-1);
const double inf=1e100;
bool deq(double a,double b) {return abs(a-b)<EPS;}

typedef complex<double> cpx;

struct Pt : cpx {
    Pt() = default;
	using cpx::cpx;
	Pt(cpx a) : cpx(a) {}
    double& x const {
        return (double&)*this;
    }
    double& y const {
        return ((double*)this)[1];
    }

	bool operator ==(cpt b) const {return abs(*this-b) < EPS; }
	bool operator <(cpt b) const {return x<b.x || (x==b.x && y<b.y); }
};

bool epsless(cpt a,cpt b) {return a.x+EPS<b.x || (deq(a.x,b.x) && a.y<b.y);}

//Allow points to be read in by input streams
istream& operator >>(istream& is, Pt& p) {
	return is >> p.x >> p.y;
}

double dot(cpt a, cpt b) {return (conj(a) * b).x;}// Dot product
double det(cpt a, cpt b) {return (conj(a) * b).y;}//Determinant/"Cross Product"
double angle(cpt a, cpt b) {return arg(b - a);}// [-pi,pi] a to b with x axis
double angle (cpt a, cpt b, cpt c) {return arg((a-b)/(c-b));}//[-pi,pi]
//double slope(cpt a, cpt b) {return tan(arg(b - a));}// m for line seg (a,b)
double slope(cpt a, cpt b) {return (b.y-a.y)/(b.x-a.x);}

Pt rotate(cpt a, double theta) {return a * polar((double)1.0, theta);}
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

/*
  vector<Pt> side[2];
  for (int i=0; i<2; i++, reverse(p.begin(), p.end())) {
  	auto it = find_if(p.begin(),p.end(),[&](auto q){ return !deq(q.x,p[0].x); });
  	side[i] = vector<Pt>(p.begin(), it), p.erase(p.begin(), it);
  	sort(side[i].begin(), side[i].end(), by_y);
  }

  for (int i=0; i<2; i++, reverse(p.begin(), p.end()))
    sort(p.begin(),find_if(p.begin(),p.end(),[&](auto q){ return !deq(q.x,p[0].x); }),
	  [](auto a, auto b) { return a.y < b.y; });
  
*/

void prt(vector<Pt> p) {
	for (Pt a:p) cerr << a << ' ';
	cerr << endl;
}

// Convex hull O(NlogN). Be careful of duplicate or very close points.
// if all points are colinear the middle points come up twice forwards and
// backwards e.g. a-b-c-d becomes a-b-c-d-c-b
// To remove colinear points change <-EPS and >EPS to <EPS and >-EPS.
vector<Pt> convexhull(vector<Pt> p) {
  sort(p.begin(),p.end(),epsless);p.resize(unique(p.begin(),p.end())-p.begin());
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

Pt intersectlineexact(Pt a, Pt b, Pt p, Pt q) {
	Pt ab=b-a,qp=p-q,ap=p-a;
	//if (deq(det(ab,qp),0)) return {1e100,1e100};
	double s=det(ap,qp)/det(ab,qp),t=det(ab,ap)/det(ab,qp);
	//double s,t;
	//tie(s,t)=solve(b-a,p-q,p-a);
	//a+t(b-a)=p+s(q-p)

	//Can also just use ptonseg.
	if (0<s && s<1 //Answer is on ab
			&& 0<t && t<1) //Answer is on pq 
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
	if (!arepara(a,b,p,q) && !std::isnan(intersectlineexact(a,b,p,q).x))
		return 0;
	
	return min({
			distfinitelinept(a,b,p),
			distfinitelinept(a,b,q),
			distfinitelinept(p,q,a),
			distfinitelinept(p,q,b)
	});
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

/*
double intersectionarea(Circle a,Circle b) {
	if (a.r<b.r) swap(
	Pt v=b.c-a.c;
	if (a.r+b.r<=abs(v)) return 0;
	if (abs(a.r-b.r)>=abs(v)) return pi*norm(min(a.r,b.r));
	double X=(norm(a.r)-norm(b.r)+norm(v))/(2.0*abs(v));
	double Ysq=norm(a.r)-norm(X),Y=Ysq<0?0:sqrt(Ysq);
*/

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

Circle circlefrom3points(cpt a,cpt b,cpt c) {
	Pt v=b-a;
	// Circle not unique or infinite if points are collinear
	if (deq(det(v,c-a),0)) return {0,-1};
	double X=abs(v);
	v/=abs(v);
	Pt p=(c-a)/v;
	Pt q(X/2,(norm(p.x)-norm(p.y)-p.x*X)/(2*p.y));
	return {q*v+a,abs(q)};
}

template<class T,int maxn>
struct Arr {
	int n=0;
	T a[maxn]={};
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
	for (int sgn=-1;sgn<2;sgn+=2) {
		Pt u=a.r+sgn*b.r;
		if (X+EPS<abs(u.x)) break;
		u.y=norm(X)-norm(u.x);
		u.y=u.y>0?sqrt(u.y):0;

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

Arr<Pt,4> tangentsthroughpt(Circle c,Pt p) {
	Arr<Pt,4> ans;
	Pt v=p-c.c;
	double X=abs(v),b=norm(X)-norm(c.r);
	v/=X;
	int &n=ans.n;
	Pt u=c.r/X*Pt(c.r,b>0?sqrt(b):0);
	if (X+EPS>c.r) {
		ans.a[n++]=u;
		ans.a[n++]={X,u.y<EPS?-X:0};
		if (u.y>=EPS) {
			ans.a[n++]=conj(u);
			ans.a[n++]={X,0};
		}
	}

	for (int i=0;i<n;i++) ans.a[i]=ans.a[i]*v+c.c;

	return ans;
}


// Find the max dot product of a point in p with v. p must be a convex polygon
// where no three points are collinear and dot(p[l],v)<dot(p[l+1],v). O(log n)
// untested
int maxdot(int l,int r,const vector<Pt>& p,cpt v) {
	if (r-l<3) {
		int i=l;
		if (l+1<r && dot(p[i],v)<dot(p[i+1],v)) i++;
		return i;
	}
	int m1=(2*l+r)/3,m2=(l+2*r)/3;
	double d1=dot(p[m1],v),d2=dot(p[m2],v);
	if (d1<dot(p[l],v)+EPS) return maxdot(l,m1,p,v);
	if (d1>d2) return maxdot(l,m2,p,v);
	return maxdot(m1+1,r,p,v);
}

// Min and max dot product of a point in p with v. p must be a convex polygon
// where no three points are collinear. Indices are returned. O(log n)
// untested
ii minmaxdot(const vector<Pt>& p,cpt v) {
	int i=deq(dot(p[0],v),dot(p[1],v));
	int n=p.size(),m,M;
	if (dot(p[i],v)<dot(p[i+1],v)) M=maxdot(i,n,p,v),m=maxdot(M,n,p,-v);
	else m=maxdot(i,n,p,-v),M=maxdot(m,n,p,v);
	for (int j=0;j<=i;j++) {
		if (dot(p[j],v)>dot(p[M],v)) M=j;
		if (dot(p[j],v)<dot(p[m],v)) m=j;
	}
	return {m,M};
}

// Find the max dot product of a point in p with v. p must be a convex polygon
// where no three points are collinear and dot(p[l],v)<dot(p[l+1],v). O(log n)
int maxdot(int l,int r,const vector<Pt>& p,cpt v) {
	if (r-l<10) {
		int i=l;
		for (int j=l+1;j<r;j++) if (dot(p[i],v)<dot(p[j],v)) i=j;
		return i;
	}
	int m1=(2*l+r)/3,m2=(l+2*r)/3;
	double d1=dot(p[m1],v),d2=dot(p[m2],v);
	if (d1<dot(p[l],v)) return maxdot(l,m1,p,v);
	if (d1>d2) return maxdot(l,m2,p,v);
	return maxdot(m1+1,r,p,v);
}

// Min and max dot product of a point in p with v. p must be a convex polygon
// where no three points are collinear. Indices are returned. O(log n)
ii minmaxdot(const vector<Pt>& p,cpt v) {
	int i=deq(dot(p[0],v),dot(p[1],v));
	int n=p.size(),m,M;
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
	int n=p.size();
	vector<Pt> r;
	for (int i=n-1,j=0;j<n;i=j++) {
		double d1=det(b-a,p[i]-a),d2=det(b-a,p[j]-a);
		if (d1>-EPS) r.push_back(p[i]);
		if ((d1>EPS && d2<-EPS) || (d1<-EPS && d2>EPS))
			r.push_back(intersectline(a,b,p[i],p[j])); //infinite lines
	}
	return r;
}

//Signed Area of polygon and circle intersection. Sign is determined by
//orientation of polygon. Divides by 0 if adjacent points are identical.
double areapolygoncircle(vector<Pt> p,Circle c) {
	int n=p.size();
	double r=0;
	for (int i=n-1,j=0;j<n;i=j++) {
		Pt v=abs(p[j]-p[i])/(p[j]-p[i]);
		Pt a=(p[i]-c.c)*v,b=(p[j]-c.c)*v;
		if (deq(a.y,0)) continue;
		double d=norm(c.r)-norm(a.y);
		if (d<0) d=0;
		d=sqrt(d);
		r+=norm(c.r)*(atan2(b.y,min(b.x,-d))-atan2(a.y,min(a.x,-d))
				+atan2(b.y,max(b.x,d))-atan2(a.y,max(a.x,d)))
			+a.y*(min(max(a.x,-d),d)-min(max(b.x,-d),d));
	}
	return r/2;
}

//Diameter of convex polygon
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

//Minimum width of a bounding rectange of a convex polygon O(n)
//The polygon must have positive signed area.
double minboundingwidth(const vector<Pt>& p) {
	double r=DBL_MAX;
	int n=p.size();
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
	vector<Pt> r;
	int n=p.size(),m=q.size();
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

// Split convex hull into lower and uppper hull. Endpoints included
pair<vector<Pt>,vector<Pt>> splithull(vector<Pt> p) {
	rotate(p.begin(),min_element(p.begin(),p.end()),p.end());
	auto it=max_element(p.begin(),p.end());
	vector<Pt> L(p.begin(),it+1),U(it,p.end());
	U.push_back(p[0]);
	reverse(U.begin(),U.end());
	return {L,U};
}

//Intersect convex polygons O(n)
//Run convex hull to remove collinear points if required
//Beware of very steep but not vertical lines when polygon coordinates can
//differ by less than EPS without being equal.
vector<Pt> intersectpolygons(const vector<Pt>&P,const vector<Pt>& Q) {
#define u(j) h[j][i[j]]
#define v(j) h[j][i[j]+1]
#define b(j) i[j]+1<h[j].size()
#define loop for (int j=0;j<4;j++)
//#define yv(j) b(j)?(u(j).y*(v(j).x-X)+v(j).y*(X-u(j).x))/(v(j).x-u(j).x):u(j).x
#define yv(j) b(j)?det(u(j)-X,v(j)-X)/(u(j).x-v(j).x):u(j).x //Test me!
	auto c=splithull(P),d=splithull(Q);
	vector<Pt> h[]{c.X,d.X,c.Y,d.Y},L,U;

	int i[4]{};
	double l=-inf,r=inf,X=-inf,nX;
	loop {
		r=min(r,h[j].back().x);
		l=max(l,u(j).x);
	}

	while (1) {
		nX=inf;
		loop if (b(j) && v(j).x>X+EPS)
			nX=min(nX,v(j).x);
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
#undef u
#undef v
#undef b
#undef loop
#undef yv
}

#define cl const Line&
struct Line {
	Pt u,v; //Endpoints
	double m() const {return (u.y-v.y)/(u.x-v.x);}
	double c() const {return yv(0);}
	double yv(double X) const {return det(u-X,v-X)/(u.x-v.x);}
	bool operator<(cl b) const {return u<b.u || (!(b.u<u) && v<b.v);}
};

namespace FindIntersection {
	const int maxn=300000;
	typedef pair<double,int> pdi;
	Line segs[maxn];
	pdi ord[2*maxn];

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

	ii checkpair(int i,int j) {
		cl a=segs[i],b=segs[j];
		Pt ab=a.v-a.u,qp=b.u-b.v,ap=b.u-a.u;
		double d1=det(ap,qp),d2=det(ab,ap),d3=det(ab,qp);
		if (d3<0) d1*=-1,d2*=-1,d3*=-1;
		if (deq(d3,0)) {// Parallel
			Pt v=ab/abs(ab);
			double c=dot(v,b.u),d=dot(v,b.v);
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
	ii findintersection(const vector<Line>& lines) {
		int n=lines.size();copy(lines.begin(),lines.end(),segs);L.clear();
		ii r;
		for (int i=0;i<n;i++) {
			if (epsless(segs[i].v,segs[i].u)) swap(segs[i].u,segs[i].v);
			ord[2*i]={segs[i].u.x,i};
			ord[2*i+1]={segs[i].v.x,i+n};
		}
		sort(ord,ord+2*n,[](const pdi& a,const pdi& b) {
				return a.X+EPS<b.X || (deq(a.X,b.X) && a.Y<b.Y);
				});
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


//For integer coordinates fitting in a square of size s the approximate minimum
//distances between distinct points are:
//	- 1/s^3 for intersection points (1/s^4 may also be possible).
//	- 1/s^2 for intersection points with an integer point
namespace BentleyOttmann {
	const int maxn=300000; // Max number of lines
	double X;
	bool post;
	bool comp(cl a,cl b) {
		double y1=a.yv(X),y2=b.yv(X),m1=a.m(),m2=b.m();
		return y1+EPS<y2 ||
			(deq(y1,y2) && ((m1!=m2 && post==(m1<m2)) || (m1==m2 && a<b)));
	}
	set<Line,bool(*)(cl,cl)> L(comp);
	typedef pair<double,int> pdi;
	priority_queue<pdi,vector<pdi>,greater<pdi>> adds,removes;
	priority_queue<Pt> inters;
	//pair<double,int> ord[2*maxn];
	
	void intersect(cl a,cl b) {
		// Can check for infinitely many points of overlap here.
		Pt poss[]{a.u,a.v,b.u,b.v};
		for (int j=0;j<4;j++) if ((j>1 &&ptonseg(a.u,a.v,poss[j]))
				|| (j<2 && ptonseg(b.u,b.v,poss[j])))
			inters.push(-poss[j]);

		if (!arepara(a.u,a.v,b.u,b.v)) {
			Pt p=intersectline(a.u,a.v,b.u,b.v);
			if (!std::isnan(p.x)) inters.push(-p);
		}
	}
	
	void addline(Line l) {
		L.insert(l);
		auto it=L.upper_bound(l);
		if (it!=L.end()) intersect(l,*it);
		if (--it!=L.begin()) intersect(l,*(--it));
	}


	vector<Pt> intersectionpoints(vector<Line> lines) {
		int n=lines.size();
		// To be continued...
	}

}

/*struct BentleyOttmann {
	const int maxn=300000;
	double X;
	bool post;
	bool comp(cl a,cl b) {
		double y1=a.yv(X),y2=b.yv(X),m1=a.m(),m2=b.m();
		return y1+EPS<y2 ||
			(deq(y1,y2) && ((m1!=m2 && post==(m1<m2)) || (m1==m2 && a<b)));
	}
	set
*/

namespace Test {
	void testargcomp() {
		vector<Pt> inp{{3,6},{-1,8},{3,-5},{6,-10},{0,-11},{-3,5},{-2,-7},{0,25},{-3,0},{0,0},{-2,0},{1,0},{0,0}};
		vector<Pt> ans{{0,0},{0,0},{-2,-7},{0,-11},{3,-5},{6,-10},{1,0},{3,6},{0,25},{-1,8},{-3,5},{-2,0},{-3,0}};
		sort(inp.begin(),inp.end(),argcomp);
		assert(inp==ans);
	}

	void testintersectfinitelinecircle() {
		Pt a(-1,sqrt(3)+5),b(1,sqrt(3)+5);
		auto ans=intersectfinitelinecircle(a,b,{{0,5},2});
		assert(ans.X==2);
		assert(ans.Y.X==a);
		assert(ans.Y.Y==b);

		a={0,sqrt(3)+5};
		ans=intersectfinitelinecircle(a,b,{{0,5},2});
		assert(ans.X==1);
		assert(ans.Y.X==b);

		a={-1,7};b={1,7};
		ans=intersectfinitelinecircle(a,b,{{0,5},2});
		assert(ans.X==1);
		a={0,7};
		assert(ans.Y.X==a);

		a={-1,5},b={-4,-6};
		ans=intersectfinitelinecircle(a,b,{{3,4},20});
		assert(ans.X==0);
	}

	void testcirclefrom3points() {
		Pt a(sqrt(3),-1),b(sqrt(2),sqrt(2)),c(-1,0),cent(-2,4);
		Circle circ=circlefrom3points(a+cent,b+cent,c+cent);
		assert(circ.c==cent);
		assert(circ.r==2);

		a=(-1,0);b=(1,0);c=(5,0);
		circ=circlefrom3points(a,b,c);
		assert(deq(circ.r,-1));
	}

	void prtcommontangents(Circle a,Circle b) {
		Arr<Pt,8> ans=commontangents(a,b);
		for (int i=0;i<ans.n;i++) cout << ans.a[i] << endl;
	}

	void testcommontangents(Circle a,Circle b,vector<Pt> expected) {
		Arr<Pt,8> ans=commontangents(a,b);
		assert(ans.n==expected.size());
		if (ans.n<=8) for (int i=0;i<ans.n;i++) assert(ans.a[i]==expected[i]);
	}

	void testtangentsthroughpt(Circle c,Pt p,vector<Pt> expected) {
		Arr<Pt,4> ans=tangentsthroughpt(c,p);
		for (int i=0;i<ans.n;i++) cerr << ans.a[i] << ' ';
		cerr << endl;
		assert(ans.n==expected.size());
		if (ans.n<=4) for (int i=0;i<ans.n;i++) assert(ans.a[i]==expected[i]);
	}

	void testareapolygoncircle() {
		vector<Pt> p{{5,0},{-8,9},{10,-6}};
		Circle c{0,3};
		assert(deq(areapolygoncircle(p,c),3.865865554267046));
	}

	
	void prtintersect(vector<Pt> p,vector<Pt> q) {
		auto ans=intersectpolygons(p,q);
		prt(ans);
	}

	vector<Pt> naivepolygonintersect(vector<Pt> p,vector<Pt> q) {
		vector<Pt> pts;
		int n=p.size(),m=q.size();
		for (int i=n-1,j=0;j<n;i=j++) for (int k=m-1,l=0;l<m;k=l++)
			if (!arepara(p[i],p[j],q[k],q[l])) {
				Pt a=intersectline(p[i],p[j],q[k],q[l]);
				if (!std::isnan(a.x) && ptinpoly(p,a) && ptinpoly(q,a)) pts.push_back(a);
			}
		for (Pt a:p) if (ptinpoly(q,a)) pts.push_back(a);
		for (Pt a:q) if (ptinpoly(p,a)) pts.push_back(a);
		return convexhull(pts);
	}

	vector<Pt> randconvex(int n, int b=30) {
		vector<Pt> pts(n);
		for (Pt &a:pts) a={rand()%b,rand()%b};
		return convexhull(pts);
	}

	Pt randpt(int b=30) {return {rand()%(2*b+1)-b,rand()%(2*b+1)-b};}

	// Convex hull O(NlogN)
	// if all points are colinear the middle points come up twice forwards and
	// backwards e.g. a-b-c-d becomes a-b-c-d-c-b
	// To remove colinear points change <-EPS and >EPS to <EPS and >-EPS.
	vector<Pt> convexhull_nocolinear(vector<Pt> p) {
	  sort(p.begin(),p.end()); p.resize(unique(p.begin(),p.end())-p.begin());
	  int l=0,u=0;
	  vector<Pt> L(p),U(p);
	  if (p.size()<=2) return p;
	  for (Pt& i:p) {
		while (l>1 && det(i-L[l-1],L[l-2]-i)<EPS) l--;
		while (u>1 && det(i-U[u-1],U[u-2]-i)>-EPS) u--;
		L[l++]=U[u++]=i;
	  }
	  L.resize(l+u-2);
	  copy(U.rend()-u+1,U.rend()-1,L.begin()+l);

	  return L;
	}

	void testminmaxdot() {
		vector<Pt> poly;
		int M=0;
		for (int n=110;n<INT_MAX;n++) for (int ct=0;ct<1000000;ct++) {
			if (ct==0) cerr << n << ' ' << M << endl;
			poly=randconvex(n, 10000);
			poly=convexhull_nocolinear(poly);
			if (poly.size() < 3) continue;
			M=max(M,(int)poly.size());
			int numdo=poly.size();
			numdo=numdo>20?20:numdo*numdo*numdo/100;
			for (int ct2=0;ct2<numdo;ct2++) {
				Pt p=randpt(10000);
				if (p==0) continue;
				ii ans=minmaxdot(poly,p); Pt minpt, maxpt;
				double biggest = -DBL_MAX/2, smallest = DBL_MAX/2;
				for (int i=0;i<(int)poly.size();i++) {
					if (dot(poly[i], p) > biggest) biggest = dot(poly[i], p), maxpt = poly[i];
					if (dot(poly[i], p) < smallest) smallest = dot(poly[i], p), minpt = poly[i];
				}
				if (!deq(smallest, dot(poly[ans.X], p))) {
					prt(poly);
					cerr << "Test point = " << p << endl;
					cerr << "Expected minimum = " << smallest << ", found " << dot(poly[ans.X], p) << endl;
					cerr << "Expected point = " << minpt << ", gave ans = " << poly[ans.X] << endl;
					return;
				}
				if (!deq(biggest, dot(poly[ans.Y], p))) {
					prt(poly);
					cerr << "Test point = " << p << endl;
					cerr << "Expected maximum = " << biggest << ", found " << dot(poly[ans.Y], p) << endl;
					cerr << "Expected point = " << maxpt << ", gave ans = " << poly[ans.Y] << endl;
					return;
				}
			}
		}
	}

	bool pequal(const vector<Pt>& p1, const vector<Pt>& p2) {
		bool good = true;
		for (int i=0; i<(int)p1.size(); i++)
			good = good && (p1[i] == p2[i]);
		return good;
	}

	Pt realintersectline(Pt a, Pt b, Pt p, Pt q) {
		Pt ab=b-a,qp=p-q,ap=p-a;
		double s=det(ap,qp)/det(ab,qp),t=det(ab,ap)/det(ab,qp);
		return a+s*ab;
	}

	vector<Pt> realconvexcut(cpt a,cpt b,const vector<Pt>& p) {
		int n=p.size();
		vector<Pt> r;
		for (int i=n-1,j=0;j<n;i=j++) {
			double d1=det(b-a,p[i]-a),d2=det(b-a,p[j]-a);
			if (d1>-EPS) r.push_back(p[i]);
			if ((d1>EPS && d2<-EPS) || (d1<-EPS && d2>EPS))
				r.push_back(realintersectline(a,b,p[i],p[j])); //infinite lines
		}
		return r;
	}

	void test_convex_cut() {
		for (int t=1; ; t++) {
			if (t % 1000 == 0) cout << "Test " << t << "           \r" << flush;
			auto poly=randconvex(5000, 1000);
			auto p1 = randpt(1000), p2 = randpt(1000);
			if (p1 == p2) continue;
			auto res1 = realconvexcut(p1,p2,poly);
			auto res2 = realconvexcut(p2,p1,poly);
			if (!deq(polygonarea(res1)+polygonarea(res2),polygonarea(poly))) {
				cerr << "poly = "; prt(poly);
				cerr << "res1 = "; prt(res1);
				cerr << "res2 = "; prt(res2);
				cerr << "cut line = " << p1 << ", " << p2 << endl;
				return;
			}
		}
	}

	void manual_test_convex_cut() {
		vector<Pt> poly = {{0,0}, {2,0}, {2,2}, {0,2}};
		vector<pair<Pt,Pt>> tests = {
			{{0,1},{1,1}},
			{{1,1},{2,0}},
			{{2,1},{1,2}},
			{{2,-1},{-1,2}},
			{{0,2},{2,2}},
			{{2,2},{0,2}},
			{{200,100}, {200,200}},
			{{-100,-100}, {-200,-100}}
		};
		vector<vector<Pt>> answers = {
			{{2,1},{2,2},{0,2},{0,1}},
			{{0,2},{2,0},{2,2}},
			{{2,1},{2,2},{1,2}},
			{{0,0},{1,0},{0,1}}
		};
	
		for (int i=0; i<(int)tests.size(); i++) {
			auto cut = realconvexcut(tests[i].X, tests[i].Y, poly);
			cout << "TEST " << (i+1) << ": "; prt(cut);
		}
	
	}
	
	void testhull() {
		vector<Pt> poly = {{0,0}, {0,1}, {0,2}, {2,2}, {2,1}, {2,0}, {1,0}};
		prt(convexhull(poly));
	}

	void wtf() {
		vector<Pt> p1 = {{0,4}, {29,2}, {29,17}, {25,21}, {0,15}},
			p2 = {{11,10}, {15,3}, {27,23}, {17,29}, {11,18}, {11,15}};
		
		auto p = naivepolygonintersect(p1, p2);
		cerr << "After convex hull: "; prt(p);
	}

	void testintersectpolygons() {
		vector<Pt> p{{2,1},{0,1},{0,0},{2,0}},q{{1,1},{1,0},{3,0},{3,1}};
		auto ans=intersectpolygons(p,q);
		for (Pt a:ans) cerr << a << ' ';
		cerr << endl;

		p={{2,-2},{3,0},{2,2},{1,0}};q={{0,0},{1,-2},{2,0},{1,2}};
		prtintersect(p,q);

		cerr << "begin" << endl;
		while (1) {
			p={},q={};
			while (deq(polygonarea(p),0)) p=randconvex(8);
			while (deq(polygonarea(q),0)) q=randconvex(8);
			vector<Pt> r=intersectpolygons(p,q),s=naivepolygonintersect(p,q);
			if (!deq(polygonarea(r),polygonarea(s))) {
				prt(p);
				prt(q);
				prt(r);
				prt(s);
				
				r=convexhull_nocolinear(r);
				s=convexhull_nocolinear(s);
				prt(r);
				prt(s);
				
				cerr << "Area(r) = " << polygonarea(r) << endl;
				cerr << "Area(s) = " << polygonarea(s) << endl;
				
				return;
			}
		}
	}

	Line randline(int m) {
		return {{rand()%m,rand()%m},{rand()%m,rand()%m}};
	}

	void prt(Line l) {
		cerr << l.u.x << ' ' << l.u.y << ' ' << l.v.x << ' ' << l.v.y << endl;
	}

	void testfindintersection() {
		while (1) {
			vector<Line> lines{randline(10),randline(10)};
			int r=FindIntersection::findintersection(lines).X;
			if (r==-1) {
				for (int i=0;i<lines.size();i++) for (int j=i+1;j<lines.size();j++)
					if (distfinitelineline(lines[i].u,lines[i].v,lines[j].u,lines[j].v)<EPS) {
						for (Line &l:lines) prt(l);
						cerr << r << endl;
						assert(0);
					}
			}
			else {
				int c=-1;
				for (int i=0;i<lines.size();i++) if (i!=r)
					if (distfinitelineline(lines[i].u,lines[i].v,lines[r].u,lines[r].v)<EPS)
						c=i;
				if (c==-1) {
					for (Line &l:lines) prt(l);
					cerr << r << endl;
					assert(0);
				}
			}
		}
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

		testargcomp();

		testintersectfinitelinecircle();

		testcirclefrom3points();

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

		testtangentsthroughpt(
				{{0,0},7},{14,0},
				{{7.0/2,7.0*sqrt(3)/2},
				{14,0},
				{7.0/2,-7.0*sqrt(3)/2},
				{14,0}}
				);
		testtangentsthroughpt(
				{{-8,23},8},{-8,15},
				{{-8,15},{-16,15}}
				);
		testtangentsthroughpt(
				{{83,29},20},{90,25},{}
				);

		testareapolygoncircle();

		//testintersectpolygons();

		testfindintersection();
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

namespace AIZU_CGL_7_H {
	void solve() {
		int n,r;
		cin >> n >> r;
		vector<Pt> p(n);
		for (Pt &a:p) cin >> a;
		cout << fixed << setprecision(10) << areapolygoncircle(p,{0,r}) << endl;
	}
}

namespace Timus_1894 {
	void solve() {
		int n,m;
		cin >> n >> m;
		vector<Pt> p(n),q(m);
		for (Pt &a:p) cin >> a;
		for (Pt &a:q) {
			cin >> a;
			a=-a;
		}
		vector<Pt> sum=minkowskisum(p,q);
		double d=DBL_MAX;
		for (int i=sum.size()-1,j=0;j<sum.size();i=j++)
			d=min(d,distfinitelinept(sum[i],sum[j],0));
		d=max((double)0.0,d-60);
		if (d<=0) d=0;
		cout << fixed << setprecision(10) << d << endl;
	}
}

namespace UVA_1111 {
	void solve() {
		int n,cas=1;
		cout << fixed << setprecision(2);
		while (cin >> n&&n) {
			vector<Pt> p(n);
			for (Pt &a:p) cin >> a;
			p=convexhull(p);
			double r=minboundingwidth(p);
			if (r<=0) r=0;
			cout << "Case " << cas++ << ": " << r+0.005 << '\n';
		}
	}
}

namespace AIZU_CGL_4_B {
	void solve() {
		int n;
		cin >> n;
		vector<Pt> p(n);
		for (Pt &a:p) cin >> a;
		cout << fixed << setprecision(10) << polygondiameter(p) << endl;
	}
}

namespace Kattis_asteroids {
	Pt a[10],b[10],v1,v2;int n,m;

	double dist(double t) {
		double r=1e200;
		for (int i=n-1,j=0;j<n;i=j++) for (int k=m-1,l=0;l<m;k=l++)
			r=min(r,distfinitelineline(a[i]+t*v1,a[j]+t*v1,b[k]+t*v2,b[l]+t*v2));
		return r;
	}

	double areaoverlap(double t) {
		vector<Pt> p(n),q(m);
		for (int i=0;i<n;i++) p[i]=a[i]+t*v1;
		for (int i=0;i<m;i++) q[i]=b[i]+t*v2;

		return polygonarea(intersectpolygons(p,q));
	}

	void solve() {
		cin >> n;
		for (int i=0;i<n;i++) cin >> a[n-1-i];
		cin >> v1;
		cin >> m;
		for (int i=0;i<m;i++) cin >> b[m-1-i];
		cin >> v2;

		if (v1==v2) {
			cout << "never\n";
			return;
		}

		double l=0,r=1e6;
		for (int ct=0;ct<200;ct++) {
			double m1=(2*l+r)/3,m2=(l+2*r)/3;
			double d1=dist(m1),d2=dist(m2);
			double a1=areaoverlap(m1),a2=areaoverlap(m2);
			if (deq(a1,0) && deq(a2,0)) {
				if (d1<d2+EPS/3) r=m2;
				else l=m1;
			}
			else {
				if (a1+EPS>a2) r=m2;
				else l=m1;
			}
		}

		if ((deq(dist(l),0) || !deq(areaoverlap(l),0)) && !(l<0))
			cout << fixed << setprecision(10) << abs(l) << endl;
		else cout << "never\n";
	}
}

namespace Timus_1469 {
	void solve() {
		int n;
		cin >> n;
		vector<Line> lines(n);
		for (Line &l:lines) cin >> l.u >> l.v;
		ii r=FindIntersection::findintersection(lines);
		if (r.X!=-1) {
			cout << "YES\n";
			cout << r.X+1 << ' ' << r.Y+1 << endl;
			return;
		}
		vector<pair<Pt,int>> ends(2*n);
		for (int i=0;i<n;i++) {
			ends[2*i]={lines[i].u,i};
			ends[2*i+1]={lines[i].v,i};
		}
		sort(ends.begin(),ends.end(),[](pair<Pt,int> a,pair<Pt,int> b) {
				return epsless(a.X,b.X) || (!epsless(b.X,a.X) && a.Y<b.Y);
				});
		for (int i=0;i+1<2*n;i++) if (ends[i].X==ends[i+1].X) {
			cout << "YES\n";
			cout << ends[i].Y+1 << ' ' << ends[i+1].Y+1 << endl;
			return;
		}
		cout << "NO\n";
	}
}

int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);

	//Test::solve();
	//UVA_920::solve();
	//UVA_10263::solve();
	//UVA_10927::solve(); // Use EPS=1e-10
	//UVA_378::solve();
	//UVA_191::solve();
	//AIZU_CGL_3_C::solve();
	//UVA_634::solve();
	//UVA_453::solve(); // Use EPS=1e-4
	//AIZU_CGL_7_H::solve();// Trim file to fit in file size limit
	//Timus_1894::solve();
	//UVA_1111::solve();
	//AIZU_CGL_4_B::solve();
	//Kattis_asteroids::solve();
	//Timus_1469::solve();
	//Test::testminmaxdot();
	//Test::test_convex_cut();
	Test::testintersectpolygons();
	//Test::wtf();
}
