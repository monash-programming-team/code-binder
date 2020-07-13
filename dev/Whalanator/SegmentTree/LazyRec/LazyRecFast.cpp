#include <bits/stdc++.h>

using namespace std;

typedef vector<int> vi;

/* 
 * O(N) construction, O(log N) queries and updates
 * Memory = O(N) = 4N([AT]+[UT]+8B)
 *
 * 0-based and ranges are inclusive.
 *
 * For N<=0 segfaults on construction
 *
 * Treats r<l as an empty range (I returned)
 * Only considers the part of the range that is in [0,N-1]
 * e.g. l=-10 and r=N+100 result is for range [0,N-1]
*/

//A Type, U Type
template <class AT, class UT>
struct SegmentTree {
	AT I,res[2],*rc,*ro; UT NU,v;//I is identity, NU is no update constant
	int N,i,j;
	vector<AT> A; vector<UT> U; vi L,R;

	//Implement these in specialisation or can inline if you only want one type
	//of segment tree
	//Then just uncomment build call
	virtual void op(AT& r, AT& i, AT& j)=0;
	virtual void us(int p, UT& v)=0;

	//Usually don't need this.
	//Use other constructor with identity array
	SegmentTree(int N, AT I, UT NU):
		I(I), NU(NU), N(N), A(4*N,I), U(4*N,NU), L(4*N), R(4*N) {
			res[0]=res[1]=I;
			build(1,0,N-1);
		}
	void build(int p, int i, int j) {
		L[p]=i;R[p]=j;
		if (i!=j) {
			build(2*p,i,(i+j)/2);
			build(2*p+1,(i+j)/2+1,j);
		}
	}
	//--------------------------------

	SegmentTree(const vector<AT>& data, AT I, UT NU):
		I(I), NU(NU), N(data.size()), A(4*N,I), U(4*N,NU), L(4*N), R(4*N) {
			/*build(1,0,N-1,data);*/
			res[0]=res[1]=I;
		}

	void build(int p, int i, int j,const vector<AT>& data) {
		L[p]=i;R[p]=j;
		if (i==j) {A[p]=data[i];return;}
		build(2*p,i,(i+j)/2,data);
		build(2*p+1,(i+j)/2+1,j,data);
		op(A[p],A[2*p],A[2*p+1]);
	}

	//Next two methods are very similar.
	void update(int p) {
		if (i<=L[p] && R[p]<=j) {return us(p,v);}
		if (R[p]<i || j<L[p]) return;

		us(2*p,U[p]);//propogate pending updates to children
		update(2*p);
		us(2*p+1,U[p]);//this happens in update and query
		update(2*p+1);
		U[p]=NU;

		op(A[p],A[2*p],A[2*p+1]);
	}

	void query(int p) {
		if (i<=L[p] && R[p]<=j) {
			op(*ro,*rc,A[p]);
			swap(ro,rc);
			return;
		}
		if (R[p]<i || j<L[p]) return;

		us(2*p,U[p]);//propogate pending updates to children
		query(2*p);
		us(2*p+1,U[p]);//this happens in update and query
		query(2*p+1);
		U[p]=NU;
	}

	void update(int l, int r, UT k) {i=l;j=r;v=k;update(1);}
	AT query(int l, int r) {
		i=l;j=r;
		res[0]=I;
		rc=res;ro=res+1;
		query(1);
		return *rc;
	}
};

/*
 * Remember to not do anything if v==NU
 * us should set both A[p] and U[p] otherwise.
 * Recognise identity
 * Call build (from specialisation) if using Array constructor
 *
 * If update can change as you propagate, just remember the original
 * range for the update and it should be possible to work out what
 * to do for any particular segment.
*/

namespace Samples {
	struct RangeMin : SegmentTree<int,int> {
		RangeMin(int N): SegmentTree(N,INT_MAX,INT_MIN) {}
		RangeMin(const vi& data): SegmentTree(data,INT_MAX,INT_MIN) {build(1,0,N-1,data);}
		void op(int& r, int& i, int& j) {r=min(i,j);}
		void us(int p, int& v) {if (v!=NU) A[p]=U[p]=v;}

	};

	struct RangeSum : SegmentTree<int,int> {
		RangeSum(int N): SegmentTree(N,0,INT_MIN) {}
		RangeSum(const vi& data): SegmentTree(data,0,INT_MIN) {build(1,0,N-1,data);}
		void op(int& r, int& i, int& j) {r=i+j;}
		void us(int p, int& v) {
			if (v!=NU) {
				U[p]=v;
				A[p]=(R[p]-L[p]+1)*v;
			}
		}
	};
}

namespace SPOJ_LITE {
	typedef pair<int,int> ii;
#define x first
#define y second

	struct STT : SegmentTree<ii,int> {
		STT(int N): SegmentTree(N,ii(0,0),0) {}
		STT(const vector<ii> data): SegmentTree(data,ii(0,0),0) {build(1,0,N-1,data);}
		void op(ii& r, ii& i, ii& j) {
			r.x=i.x+j.x;
			r.y=i.y+j.y;
		}
		void us(int p, int& v) {
			if (v!=0) {
				U[p]^=1;
				A[p].x=A[p].y-A[p].x;
			}
		}
	};

	int N,M,t,i,j;
	void solve() {
		cin >> N >> M;
		STT st(vector<ii>(N,ii(0,1)));
		while (M--) {
			cin >> t >> i >> j;i--;j--;
			if (t) cout << st.query(i,j).x << '\n';
			else st.update(i,j,1);
		}
	}
#undef x
#undef y
}

int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);
	//RangeMin rm(20);
	//cout << rm.query(0,19) << endl;

	SPOJ_LITE::solve();
}
