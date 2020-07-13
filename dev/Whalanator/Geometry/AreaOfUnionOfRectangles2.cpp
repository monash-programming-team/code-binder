#include <bits/stdc++.h>

#define X first
#define Y second

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

//listings:union_of_rect
// Area of union of rectangles. Include SegmentTree code. Complexity: O(N log(N))
template<class T> struct UnionOfRect {
	int m=0,U=0;  T nm=1,l=1;
	void op(UnionOfRect& b,UnionOfRect& c) {
		if (b.m<c.m) m=b.m,	nm=b.nm;
		else { m=c.m, nm=c.nm; if (b.m==c.m) nm+=b.nm; }
		l=b.l+c.l;
	}
	void us(int v) { m+=v, U+=v; }
	void NU() {U=0;}
	T nonzerolen() { return l-(m?0:nm);	}
};

template<class T> struct Rect {	pair<T,T> l,u; }; // lower left and upper right corners

// You will get runtime error if all y values are the same. 
template<class T> T areaofunionofrect(const vector<Rect<T>>& rect) {
	int n=rect.size(),m;  T r=0; 
	vector<T> ys(2*n);	vector<pair<pair<T,int>,pair<T,T>>> sides(2*n);
	for (int i=0;i<n;i++) {
		ys[2*i]=rect[i].l.Y, ys[2*i+1]=rect[i].u.Y;
		sides[2*i]={{rect[i].l.X,1},{rect[i].l.Y,rect[i].u.Y}};
		sides[2*i+1]={{rect[i].u.X,-1},{rect[i].l.Y,rect[i].u.Y}};
	}
	sort(ys.begin(),ys.end());  sort(sides.begin(),sides.end());
	ys.resize(unique(ys.begin(),ys.end())-ys.begin());
	vector<UnionOfRect<T>> stinit(m=ys.size()-1);
	for (int i=0;i<m;i++) stinit[i].l=stinit[i].nm=ys[i+1]-ys[i];
	SegmentTree<UnionOfRect<T>,int> st(stinit);  // Include SegmentTree code
	T x = sides[0].X.X;
	for (auto &i:sides) {
		r+=(i.X.X-x)*st.query(0,m).nonzerolen();  x=i.X.X;
		int a=lower_bound(ys.begin(),ys.end(),i.Y.X)-ys.begin();
		int b=lower_bound(ys.begin(),ys.end(),i.Y.Y)-ys.begin();
		if (a!=b) st.update(a,b,i.X.Y);
	}
	return r;
}
//listings:/union_of_rect

namespace SPOJ_AREAU {
	void solve() {
		int n;
		cin >> n;
		vector<Rect<int>> r(n);
		for (auto &re:r) cin >> re.l.X >> re.l.Y >> re.u.X >> re.u.Y;
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
			for (auto &a:r) cin >> a.l.X >> a.l.Y >> a.u.X >> a.u.Y;
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
