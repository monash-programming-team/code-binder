// Extreme visible points on a polygon query
//
// Author: Daniel Anderson
// Date: 30-04-2017
// Reliability: 2
// Tested on: Brute-force, PNW2016:E
//
// Given a convex polygon, allows queries of the form, what range of vertices
// are visible to a given point that is strictly outside the polygon.

#include<bits/stdc++.h>
using namespace std;

#include "../code-template/debug.h"

#define X first
#define Y second

typedef vector<int> vi;
typedef pair<int,int> pii;

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
double dot(cpt a, cpt b) {return (conj(a) * b).x;}
double det(cpt a, cpt b) {return (conj(a) * b).y;}

int orient(cpt a, cpt b, cpt c) {
	double d=det(b-a,c-b);
	return d>EPS?1:d<-EPS?-1:0;
}

bool argcomp(cpt a,cpt b) {
	if (b==0) return 0;
	if (a==0) return 1;
	double a1=arg(a),a2=arg(b);
	if (a1<-pi+EPS/2) a1+=2*pi;
	if (a2<-pi+EPS/2) a2+=2*pi;
	return a1+EPS<a2 || (deq(a1,a2) && norm(a)<norm(b));
}

// Facilitates queries for the pair of furthest visible points on a convex polygon from
// some external point. P must be a non-degenerate convex polygon. Returns an empty
// interval for p inside P. Complexity: O(N log(N)) pre-process, O(log(N)) per query.
struct ExtremePoints {
  vector<Pt> P; Pt c; int n; vi ids;
  ExtremePoints(vector<Pt> poly) : P(move(poly)), n(P.size()), ids(2*n) {
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

namespace Test {
  
  bool epsless(cpt a,cpt b) {return a.x+EPS<b.x || (deq(a.x,b.x) && a.y<b.y);}
  
  vector<Pt> convexhull(vector<Pt> p) {
    sort(p.begin(),p.end(),epsless); p.resize(unique(p.begin(),p.end())-p.begin());
    int l=0,u=0; vector<Pt> L(p),U(p);
    if (p.size()<=2) return p;
    for (Pt& i:p) {
      while (l>1 && det(i-L[l-1],L[l-2]-i)<EPS) l--;
      while (u>1 && det(i-U[u-1],U[u-2]-i)>-EPS) u--;
      L[l++]=U[u++]=i;
    }
    L.resize(l+u-2); copy(U.rend()-u+1,U.rend()-1,L.begin()+l);
    return L;
  }
  
  vector<Pt> randconvex(int n, int b=30) {
    vector<Pt> pts(n);
    for (Pt &a:pts) a={double(rand()%b),double(rand()%b)};
    return convexhull(pts);
  }

  Pt randpt(int b=30) {return {double(rand()%(2*b+1)-b),double(rand()%(2*b+1)-b)};}

  pii naive_extreme(vector<Pt> P, const Pt& p) {
    int n = (int)P.size();
    Pt c; for (auto p : P) c += 1.0/n * p;
    for (auto& p : P) p -= c; sort(P.begin(), P.end(), argcomp);
    auto seen = [&](int i) { return orient(P[i?i-1:n-1],P[i],p-c) < 0; };
    if (seen(0)) {  // circular
      int a = n-1; while(seen(a)) a--;
      int b = 0; while(seen(b+1)) b++;
      return {a,b};
    } else {  // contiguous
      int a = 0; while (!seen(a+1)) a++;
      int b = n-1; while (!seen(b)) b--;
      return {a,b};
    }
  }
  
  void prt(vector<Pt> p) {
    for (Pt a:p) cerr << a << ' ';
    cerr << endl;
  }
  
  double polygonarea(const vector<Pt>& p) {
    double r=0;	int n=p.size();
    for (int j=0,i=n-1;j<n;i=j++) r+=det(p[i],p[j]);
    return r/2;
  }
  
  bool ptonseg(cpt a, cpt b, cpt p) {
    Pt u=b-a,v=p-a;
    return a==p || b==p || ((0<dot(u,v) && dot(u,v)<norm(u)) && deq(det(u,v),0));
  }
    
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
  
  void test() {
    int tc=0;
    for (int n=5; ; n++) {
      for (int tt=0;tt<10000;tt++) {
        auto P = randconvex(n, 50);
        auto p = randpt(100);
        if (abs(polygonarea(P)) < EPS) continue;  // Degenerate polygon
        if (tc++ % 100 == 0) cout << "Test " << tc++ << ", n="<<n << "              \r" << flush;
        ExtremePoints xtreme(P);
        auto ans2 = xtreme.query(p);
        if (point_in_convex_polygon(P, p, true)) {
          if (ans2.X != ans2.Y) {
            cerr << "Internal point not returning empty interval" << endl;
            DEBUG(P);
            DEBUG(p);
            DEBUG(ans2);
            assert(ans2.X == ans2.Y);
          }
        } else {
          auto ans1 = naive_extreme(P,p);
          if (ans1 != ans2) {
            DEBUG(P);
            DEBUG(p);
            DEBUG(ans1, ans2);
            assert(ans1 == ans2);
          }
        }
      }
    }
  }
}

int main() {
  vector<Pt> P = { {0,0}, {1,0}, {1,1}, {0,1} };
  ExtremePoints ext(P);
  Pt p = {0.25, 0.5};
  DEBUG(ext.query(p));
  
  Test::test();
}
