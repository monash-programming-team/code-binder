#include <bits/stdc++.h>

using namespace std;

typedef vector<int> vi;

/* 
 * O(N) construction, O(log N) queries and updates
 * Memory = O(N) = 4N([AT]+[UT]+8B)
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
	AT I,r; UT NU,v;//I is identity, NU is no update constant
	int N,i,j;
	vector<AT> A; vector<UT> U; vi L,R;

	//Implement these in specialisation or can inline if you only want one type of segment tree
	//Then just uncomment build call
	virtual AT op(const AT& i, const AT& j)=0;
	virtual AT& us(int p, const UT& v)=0;

	//Usually don't need this.
	//Use other constructor with identity array
	SegmentTree(int N, AT I=0, UT NU=INT_MIN):
		I(I), NU(NU), N(N), A(4*N,I), U(4*N,NU), L(4*N), R(4*N) {build(1,0,N-1);}
	void build(int p, int i, int j) {
		L[p]=i;R[p]=j;
		if (i!=j) {
			build(2*p,i,(i+j)/2);
			build(2*p+1,(i+j)/2+1,j);
		}
	}
	//--------------------------------

	SegmentTree(const vector<AT>& data, AT I=0, UT NU=INT_MIN):
		I(I), NU(NU), N(data.size()), A(4*N), U(4*N,NU), L(4*N), R(4*N) {/*build(1,0,N-1,data);*/}

	AT build(int p, int i, int j,const vector<AT>& data) {
		L[p]=i;R[p]=j;
		if (i==j) return A[p]=data[i];
		return A[p]=op(build(2*p,i,(i+j)/2,data),build(2*p+1,(i+j)/2+1,j,data));
	}

	AT update(int p) {
		if (i<=L[p] && R[p]<=j) return us(p,v);
		if (R[p]<i || j<L[p]) return v==NU?I:A[p];
		
		us(2*p,U[p]);//propagate pending updates to children
		us(2*p+1,U[p]);//only place this happens
		U[p]=NU;

		r=op(update(2*p),update(2*p+1));
		if (v!=NU) A[p]=r;
		return r;
	}

	AT update(int l, int r, UT k) {i=l;j=r;v=k;return update(1);}
	AT query(int l, int r) {return update(l,r,NU);}
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

struct RangeMin : SegmentTree<int,int> {
	RangeMin(int N): SegmentTree(N,INT_MAX) {}
	RangeMin(const vi& data): SegmentTree(data,INT_MAX) {build(1,0,N-1,data);}
	int op(const int& i,const int& j) {return min(i,j);}
	int& us(int p, const int& v) {return v==NU?A[p]:A[p]=U[p]=v;}
};

struct RangeSum : SegmentTree<int,int> {
	RangeSum(int N): SegmentTree(N) {}
	RangeSum(const vi& data): SegmentTree(data) {build(1,0,N-1,data);}
	int op(const int& i, const int& j) {return i+j;}
	int& us(int p, const int& v) {
		if (v==NU) return A[p];
		U[p]=v;
		return A[p]=(R[p]-L[p]+1)*v;
	}
};

int main() {
	RangeMin rm(20);
	cout << rm.query(0,19) << endl;
	return 0;
}
