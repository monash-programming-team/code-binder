#include <bits/stdc++.h>

using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;

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

//Type of segment tree (e.g. RangeMin), A type, U type
template<class T,class AT,class UT>
struct HeavyLight {
	int N;
	vi O,P,H,D; //Left start in segment tree for vertex. Right..., Parent, Head, Depth
	T st;

	//HeavyLight(int N): N(N), L(N), R(N), H(N), D(N), st(N)
	HeavyLight( vvi& al,const vector<AT>& data):
		N(al.size()), O(N),H(N),P(N),D(N),st(N) {build(al,data);}

	/*int ma(int p,int par,vvi & al) {//modify adjacency list placing heavy edge first.
	  P[p]=par;
	  int i,Mc,M=0,S=1,t;
	  for (int c=0;c<al[p].size();c++) if ((i=al[p][c])!=par) {
	  t=ma(i,p,al);
	  S+=t;
	  if (t>M) {
	  Mc=c;
	  M=t;
	  }
	  }
	  if (M) swap(al[p][0],al[p][Mc]);
	  return S;
	  }*/

	int setPH(int p,int par,const vvi& al) {
		P[p]=par;
		H[p]=p;
		int i,M=0,S=1,t;
		for (int c:al[p]) if (c!=par) {
			t=setPH(c,p,al);
			S+=t;
			if (t>M) {
				i=c;
				M=t;
			}
		}
		if (M) H[i]=p;
		return S;
	}

	void build(const vvi& al,const vector<AT>& data) {
		vector<AT> rord;//Reordered data
		setPH(0,0,al);
		build(0,0,0,al,data,rord);
		st.build(1,0,N-1,rord);
	}

	void build(int p,int h,int d,const vvi& al,const vector<AT>& data,vector<AT>& rord) {
		O[p]=rord.size();
		rord.push_back(data[p]);
		//L[p]=p==h?rord.size():L[h]; //To get actual value for non heads
		H[p]=h;
		D[p]=d;
		for (int c:al[p]) if (c!=P[p] && H[c]==p) build(c,h,d+1,al,data,rord);
		for (int c:al[p]) if (c!=P[p] && H[c]==c) build(c,c,d+1,al,data,rord);
	}

	int lca(int i,int j) {
		while (H[i]!=H[j]) if (D[H[i]]<D[H[j]]) i=P[H[i]];
		else j=P[H[j]];
		return D[i]<D[j]?i:j;
	}

	AT update(int i,int j,UT v) {
		AT r=st.I;
		while (H[i]!=H[j]) if (D[H[i]]<D[H[j]]) {
			r=st.op(r,st.update(O[H[i]],O[i],v));
			i=P[H[i]];
		}
		else {
			r=st.op(r,st.update(O[H[j]],O[j],v));
			j=P[H[j]];
		}
		if (D[i]>D[j]) swap(i,j);
		return st.op(r,st.update(O[i],O[j],v));
	}

	AT query(int i,int j) {return update(i,j,st.NU);}
};

int main() {
	return 0;
}

