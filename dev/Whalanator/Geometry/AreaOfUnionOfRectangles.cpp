#include <bits/stdc++.h>

#define x first
#define y second

using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef pair<int,int> ii;

template <class T, class U>
struct SegmentTree {
	T I,t[4]; int N,h;
	vector<T> A;

	SegmentTree(const vector<T>& data, T I=T()): I(I), N(data.size()),
		h(sizeof(int)*8-__builtin_clz(N)), A(2*N,I) {
			copy(data.begin(),data.end(),A.begin()+N);
			for (int i=N-1;i;i--) op(i);
	}

	void op(int i) {A[i].op(A[2*i],A[2*i+1]);}

	void prop(int i) {
		A[2*i].us(A[i].U);
		A[2*i+1].us(A[i].U);
		A[i].NU();
	}

	void push(int i) {for (int j=h;j;j--) prop(i>>j);}

	void update(int l, int r, U v) {
		push(l+=N);
		push((r+=N)-1);
		bool cl=0,cr=0;
		for (;l<r;l/=2,r/=2) {
			if (cl) op(l-1);
			if (cr) op(r);
			if (l&1) A[l++].us(v), cl=1;
			if (r&1) A[--r].us(v), cr=1;
		}

		if (l==1 && cr) op(1);
		else for (l--;r>0;l/=2,r/=2) {
			if (cl && l) op(l);
			if (cr && (!cl || (l!=r && r!=1))) op(r);
		}
	}

	T query(int l, int r) {
		push(l+=N);
		push((r+=N)-1);
		t[0]=t[2]=I;
		int i=0,j=2;
		for (;l<r;l/=2,r/=2) {
			if (l&1) t[i^1].op(t[i],A[l++]), i^=1;
			if (r&1) t[j^1].op(A[--r],t[j]), j^=1;
		}

		t[i^1].op(t[i],t[j]);
		return t[i^1];
	}
};

template<class T>
struct UnionOfRect {
	int m=0,U=0;
	T nm=1,l=1;

	void op(UnionOfRect& b,UnionOfRect& c) {
		if (b.m<c.m) {
			m=b.m;
			nm=b.nm;
		}
		else {
			m=c.m;
			nm=c.nm;
			if (b.m==c.m) nm+=b.nm;
		}
		l=b.l+c.l;
	}

	void us(int v) {
		m+=v;
		U+=v;
	}

	void NU() {U=0;}

	T nonzerolen() {
		return l-(m?0:nm);
	}
};

template<class T>
struct Rect {
	pair<T,T> l,u; // lower left corner and upper right corner
};

// You will get runtime error if all y values are the same
// Segment tree becomes size 0.
template<class T>
T areaofunionofrect(const vector<Rect<T>>& rect) {
	T r=0;
	int n=rect.size(),m;
	vector<T> ycoord(2*n);
	vector<pair<pair<T,int>,pair<T,T>>> sides(2*n);
	for (int i=0;i<n;i++) {
		ycoord[2*i]=rect[i].l.y;
		ycoord[2*i+1]=rect[i].u.y;
		sides[2*i]={{rect[i].l.x,1},{rect[i].l.y,rect[i].u.y}};
		sides[2*i+1]={{rect[i].u.x,-1},{rect[i].l.y,rect[i].u.y}};
	}
	sort(ycoord.begin(),ycoord.end());
	ycoord.resize(unique(ycoord.begin(),ycoord.end())-ycoord.begin());

	sort(sides.begin(),sides.end());
	
	vector<UnionOfRect<T>> stinit(m=ycoord.size()-1);
	for (int i=0;i<m;i++) stinit[i].l=stinit[i].nm=ycoord[i+1]-ycoord[i];
	SegmentTree<UnionOfRect<T>,int> st(stinit);

	T x=INT_MIN; // Less than or equal to min x value
	for (auto &i:sides) {
		r+=(i.x.x-x)*st.query(0,m).nonzerolen();
		x=i.x.x;
		int a=lower_bound(ycoord.begin(),ycoord.end(),i.y.x)-ycoord.begin();
		int b=lower_bound(ycoord.begin(),ycoord.end(),i.y.y)-ycoord.begin();
		if (a!=b) st.update(a,b,i.x.y);
	}

	return r;
}

namespace SPOJ_AREAU {
	void solve() {
		int n;
		cin >> n;
		vector<Rect<int>> r(n);
		for (auto &re:r) cin >> re.l.x >> re.l.y >> re.u.x >> re.u.y;
		cout << areaofunionofrect(r) << endl;
	}
}

namespace LIVE_ARCHIVE_2184 {
	void solve() {
		int n;
		cout << fixed << setprecision(2);
		int cas=1;
		while (cin >> n && n) {
			vector<Rect<double>> r(n);
			for (auto &a:r) cin >> a.l.x >> a.l.y >> a.u.x >> a.u.y;
			cout << "Test case #" << cas++ << '\n';
			cout << "Total explored area: " << areaofunionofrect(r) << "\n\n";
		}
	}
}

namespace UVA_688 {
	void solve() {
		int n;
		cout << fixed << setprecision(2);
		int cas=1;
		while (cin >> n && n) {
			vector<Rect<double>> r;
			for (int i=0;i<n;i++) {
				double x,y,ra;
				cin >> x >> y >> ra;
				if (ra>1e-12) r.push_back({{x-ra,y-ra},{x+ra,y+ra}});
			}
			n=r.size();
			cout << cas++ << ' ' << (!n?0.0:areaofunionofrect(r)) << '\n';
		}
	}
}

int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);

	//SPOJ_AREAU::solve(); //This problem won't let you submit to it.
	//LIVE_ARCHIVE_2184::solve();
	UVA_688::solve();
}
