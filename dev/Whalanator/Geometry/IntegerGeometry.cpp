#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef pair<int,int> ii;

#define x real()
#define y imag()
#define X first
#define Y second
#define cpt const Pt&

typedef complex<ll> cpx;

struct Pt : cpx {
	Pt() = default;
	using cpx::cpx;
	Pt(cpx a) : cpx(a) {}
	ll& x const {
		return (ll&)*this;
	}
	ll& y const {
		return ((ll*)this)[1];
	}

	bool operator<(cpt b) {return x<b.x || (x==b.x && y<b.y);}
};

istream& operator>>(istream& is,Pt& p) {
	return is >> p.x >> p.y;
}

ll dot(cpt a,cpt b) {return (conj(a)*b).x;}
ll det(cpt a,cpt b) {return (conj(a)*b).y;}

bool argcomp(cpt a,cpt b) {
	if (b==(ll)0) return 0;
	if (a==(ll)0) return 1;
	bool r1=a.y>0 || (a.y==0 && a.x<0);
	bool r2=b.y>0 || (b.y==0 && b.x<0);
	ll d=det(a,b);
	return r1<r2 || (r1==r2 && (d>0 || (d==0 && norm(a)<norm(b))));
}

namespace Test {
	void testargcomp() {
		vector<Pt> inp{{3,6},{-1,8},{3,-5},{6,-10},{0,-11},{-3,5},{-2,-7},{0,25},{-3,0},{0,0},{-2,0},{1,0},{0,0}};
		vector<Pt> ans{{0,0},{0,0},{-2,-7},{0,-11},{3,-5},{6,-10},{1,0},{3,6},{0,25},{-1,8},{-3,5},{-2,0},{-3,0}};
		sort(inp.begin(),inp.end(),argcomp);
		assert(inp==ans);
	}

	void solve() {
		testargcomp();
	}
}

int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);

	Test::solve();
}
