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
	bool operator <(const Pt& b) const {return x<b.x || (x==b.x && y<b.y); }
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
    while (l>1 && det(i-L[l-1],L[l-2]-i)<EPS) l--;
    while (u>1 && det(i-U[u-1],U[u-2]-i)>-EPS) u--;
    L[l++]=U[u++]=i;
  }
  L.resize(l+u-2);
  copy(U.rend()-u+1,U.rend()-1,L.begin()+l);

  return L;
}

int rnd(double k) {
	return k>=0?k+0.5:k-0.5;
}

int main() {
	int T;
	cin >> T;
	cout << T << endl;
	while (T--) {
		int n;
		cin >> n;
		vector<Pt> pts(n);
		for (Pt &p:pts) cin >> p;
		vector<Pt> hull = convexhull(pts);
		int i=min_element(hull.begin(),hull.end(),[](Pt a,Pt b){
				return a.y<b.y-EPS || (dequal(a.y,b.y) && a.x < b.x);
				})-hull.begin();

		cout << hull.size()+1<< endl;
		for (int c=0;c<=hull.size();c++)
			cout << rnd(hull[(i+c)%hull.size()].x) << ' ' <<
				rnd(hull[(i+c)%hull.size()].y) << '\n';
		if (T) cout << "-1\n";
		cin >> i;
	}
}
