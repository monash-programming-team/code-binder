#include <bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef pair<int,int> ii;

//listings:geometry
// -------------------------- Integer Computational Geometry ---------------------------
#define x real()
#define y imag()
#define cpt const Pt&

typedef complex<ll> cpx;
struct Pt : cpx {
	Pt() = default; using cpx::cpx;
	Pt(cpx a) : cpx(a) {}
	ll& x const { return (ll&)*this; }
	ll& y const { return ((ll*)this)[1]; }
	bool operator<(cpt b) {return x<b.x || (x==b.x && y<b.y);}
};

ll dot(cpt a,cpt b) { return (conj(a)*b).x; }
ll det(cpt a,cpt b) { return (conj(a)*b).y; }

//Compare points by principal argument (-pi,pi] breaking ties by norm.
//0 is considered less than everything else.
bool argcomp(cpt a, cpt b) {
	if (b==(ll)0) return 0;
	if (a==(ll)0) return 1;
	bool r1=a.y>0 || (a.y==0 && a.x<0), r2=b.y>0 || (b.y==0 && b.x<0); ll d=det(a,b);
	return r1<r2 || (r1==r2 && (d>0 || (d==0 && norm(a)<norm(b))));
}
//listings:/geometry

int main() {

}
