#include<bits/stdc++.h>

using namespace std;

// Can't use x and y for first and second
#define x real()
#define y imag()
//#define Pt Pt//const Pt&

const double EPS = 1e-9;
const double pi=acos(-1);
bool dequal(double a,double b) {return abs(a-b)<EPS;}
enum Orientation {CCW, CW, CNEITHER};


//using namespace std;
typedef complex<double> cpx;

struct Pt : public cpx {
    Pt() = default;
	using cpx::cpx;
	Pt(cpx a) : cpx(a) {}
    double& real() const {
        return (double&)*this;
    }
    double& imag() const {
        return ((double*)this)[1];
    }

	bool operator ==(const Pt& b) const {return abs(*this-b) < EPS; }
	bool operator <(const Pt& b) const {return y<b.y || (y==b.y && x<b.x); }
	//by y first then break ties with x
};

//Allow points to be read in by input streams
istream& operator >>(istream& is, Pt& p) {
	is >> p.x >> p.y;
	return is;
}

double dot(Pt a, Pt b) {return (conj(a) * b).x;}// Dot product
double det(Pt a, Pt b) {return (conj(a) * b).y;}//Determinant/"Cross Product"
//double dist(Pt a, Pt b) {return abs(a - b);}// Euclidean distance from a to b
double dist2(Pt a, Pt b) {return norm(a - b);}// Euclidean distance squared
double angle(Pt a, Pt b) {return arg(b - a);}// [-pi,pi] a to b with x axis
double angle (Pt a, Pt b, Pt c) {return arg((a-b)/(c-b));}//[-pi,pi]
double slope(Pt a, Pt b) {return tan(arg(b - a));}// m for line segment (a,b)

Pt rotate(Pt a, double theta) {return a * polar(1.0, theta);}//anticlockwise
//around p by theta anticlockwise
Pt rotate(Pt a, Pt p, double theta) {return rotate(a - p,theta) + p;}
Pt project(Pt p, Pt v) {return v * dot(p, v) / norm(v);}// p onto v
Pt project(Pt p, Pt a, Pt b) {return a+project(p-a,b-a);}//p onto line (a,b)
//reflect p across the line (a,b)
Pt reflect(Pt p, Pt a, Pt b) {return a + conj((p - a) / (b - a)) * (b - a);}

//bool colinear(Pt a, Pt b, Pt c) {return dequal(det(b-a,c-b),0);}

// Orientation test (1 anticlockwise, -1 clockwise, 0 colinear)
int orient(Pt a, Pt b, Pt c) {
	double d=det(b-a,c-b);
	return d>EPS?1:d<-EPS?-1:0;
}

// Point on line segment (including endpoints)
bool ptonseg(Pt a, Pt b, Pt p) {
	Pt u=b-a,v=p-a;
	return a==p || b==p ||
		((0 < dot(u,v) && dot(u,v) < norm(u)) && dequal(det(u,v),0));
}

// Signed area of polygon
// Positive for anticlockwise orientation
double polygonarea(const vector<Pt>& p) {
	double r=0;
	int n=p.size();
	for (int i=0;i<n;i++) r+=det(p[i],p[(i+1)%n]);
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
const bool boundary=1; //is boundary in polygon?
bool ptinpoly(vector<Pt> p, Pt q) {
	int n=p.size();
	int i,j,r=0;
	for (j=0,i=n-1;j<n;i=j++) { //that trick to avoid modding
		if (ptonseg(p[i],p[j],q)) return boundary;
		if (((p[i].y <= q.y && q.y < p[j].y)
			|| (p[j].y <= q.y && q.y < p[i].y))
			&& q.x < (p[j].x-p[i].x) * (q.y-p[i].y)/(p[j].y-p[i].y) + p[i].x)
			r=!r;
	}
	return r;
}

Pt solve(Pt a, Pt b, Pt v) {// solves [a b]x==v with Cramer's rule.
	return Pt(det(v,b)/det(a,b),det(a,v)/det(a,b));
}

//Distance between infinite line and point.
double linept(Pt a, Pt b, Pt p) {
	return abs(det(b-a,p-a)/abs(b-a));
}

//Distance between finite line and point (probably works in higher dimensions)
double flinept(Pt a, Pt b, Pt p) {
	b-=a;p-=a;
	double sp=dot(b,p)/norm(b);
	Pt closest;
	if (sp>=0) {
		if (sp>1) closest=b;
		else closest=sp*b;
	}
	return abs(closest-p);
}



//Are lines parallel?
bool arepara(Pt a, Pt b, Pt p, Pt q) {
	return dequal(det(b-a,q-p),0);
}

//Distance between 2 finite lines
double flineline(Pt a,Pt b,Pt p,Pt q) {
	if (arepara(a,b,p,q)) {
		b-=a;p-=a;q-=a;
		double sp=dot(b,p)/norm(b);
		if (0<sp && sp<1) return det(b,p)/abs(b);
	}
	return min({flinept(a,b,p),flinept(a,b,q),flinept(p,q,a),flinept(p,q,b)});
}

//This is kind of unnecessary
double distpara(Pt a, Pt b, Pt p, Pt q) {
	return linept(a,b,p);
}

//Intersection of 2 line segments. Divides by 0 if they are parallel.
//Returns {nan,nan} if they don't intersect.
//Uncomment if statements below to get infinite lines.
Pt intersectline(Pt a, Pt b, Pt p, Pt q) {
	Pt ab=b-a,qp=p-q,ap=p-a;
	double t=det(ab,ap)/det(ab,qp),s=det(ap,qp)/det(ab,qp);
	//a+t(b-a)=p+s(q-p)
	
	//Can also just use ptonseg.
	if (-EPS<t && t<1+EPS //Answer is on ab
		&& -EPS<s && s<1+EPS) //Answer is on pq 
		return a+t*ab;
	return Pt(NAN,NAN);
}


/*
// Intersection of lines (a,b) and (p,q)
Pt intersection(Pt a, Pt b, Pt p, Pt q) {		
  double c1 = cross(p - a, b - a), c2 = cross(q - a, b - a);
  return (c1 * q - c2 * p) / (c1 - c2); // undefined if parallel
}
*/




int main() {
	int N,M,cas=1;
	while (cin >> N && N) {
		cin >> M;
		vector<Pt> poly(N);
		for (Pt &p:poly) cin >> p;
		vector<Pt> quer(M);
		for (Pt &p:quer) cin >> p;
		if (cas!=1) printf("\n");
		printf("Problem %d:\n",cas++);
		//poly=convexhull(poly);
		for (Pt p:quer) {
			if (ptinpoly(poly,p)) printf("Within\n");
			else printf("Outside\n");
		}
	}
}
