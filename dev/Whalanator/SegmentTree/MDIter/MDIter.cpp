#include <bits/stdc++.h>

using namespace std;

/*
 * Multidimensional Segment Tree
 *
 * O(2^D * N1N2...ND) construction
 * O(log(N1)log(N2)...log(ND)) queries and updates
 * Memory = 2^D * N1N2...ND * [AT]
 *
 * For N<=0 segfaults on build call
 *
 * Treats ri<li as an empty range (I returned)
 * Generally segfaults for other invalid ranges. There are special
 * cases where wrong answer is received instead.
 *
 * Non commutative problems make no sense in more than one dimension
 * and are not supported.
 *
 * The using feature may not be supported on compilers that use c++0x.
*/

//Multidimensional vector
template<int D,class T>
struct Vec {typedef vector<typename Vec<D-1,T>::type> type;};
template<class T>
struct Vec<0,T> {typedef T type;};

template<int D,class T>
using MDV = typename Vec<D,T>::type;

//Segment tree uses static members of problem it solves

struct RangeSum {
	static const int I=0;
	static int op(const int& x,const int& y) {return x+y;}
};

template<class T,class AT,int D>
struct SegmentTree {
	template<int _D>
	using ST = SegmentTree<T,AT,_D>;

	int N;
	vector<ST<D-1> > A;

	template<class...Ts>
	SegmentTree(int N,Ts...Ns): N(N), A(2*N,ST<D-1>(Ns...)) {}

	void build(const MDV<D,AT>& data) {
		for (int c=0;c<N;c++) A[c+N].build(data[c]);
		for (int c=N-1;c;c--) A[c].merge(A[2*c],A[2*c+1]);
	}

	void merge(const ST<D>& L, const ST<D>& R) {
		for (int c=1;c<2*N;c++) A[c].merge(L.A[c],R.A[c]);
	}

	template<class...Ts>
	void merge(const ST<D>& L, const ST<D>& R, int i, Ts...is) {
		for (i+=N;i;i/=2) A[i].merge(L.A[i],R.A[i],is...);
	}

	template<class...Ts>
	void update(const AT& v, int i, Ts...is) {
		for (A[i+=N].update(v,is...);i/=2;) A[i].merge(A[2*i],A[2*i+1],is...);
	}

	template<class...Ts>
	AT query(int i,int j,Ts...limits) {
		AT r=T::I;
		for (i+=N,j+=N;i<=j;i/=2,j/=2) {
			if (i&1) r=T::op(r,A[i++].query(limits...));
			if (~j&1) r=T::op(r,A[j--].query(limits...));
		}
		return r;
	}
};

template<class T,class AT>
struct SegmentTree<T,AT,0> {
	typedef SegmentTree<T,AT,0> ST;
	AT a;
	SegmentTree() {a=T::I;}
	void build(const AT& data) {a=data;}
	void merge(const ST& L,const ST& R) {a=T::op(L.a,R.a);}
	void update(const AT& v) {a=v;}
	AT query() {return a;}
};

const int numtests=100000,numqueries=10000;

int main() {
	srand(time(0));
	for (int t=0;t<numtests;t++) {
		//printf("Test: %d\n",t);
		//cerr << t << endl;
		vector<int> Ns(5);
		for (int c=0;c<5;c++) /*printf("%d\n",*/Ns[c]=rand()%25+1/*)*/;
		MDV<5,int> A{Ns[0],{Ns[1],{Ns[2],{Ns[3],vector<int>(Ns[4])}}}};

		for (int c=0;c<Ns[0];c++)
			for (int d=0;d<Ns[1];d++)
				for (int e=0;e<Ns[2];e++)
					for (int f=0;f<Ns[3];f++)
						for (int g=0;g<Ns[4];g++)
							/*printf("%d\n",*/A[c][d][e][f][g]=rand()%401-200/*)*/;

		SegmentTree<RangeSum,int,5> st(Ns[0],Ns[1],Ns[2],Ns[3],Ns[4]);
		st.build(A);

		vector<int> pt(5),bounds(10);
		int v;

		for (int q=0;q<numqueries;q++) {
			int type=rand()%2;
			//printf("%d",type);
			if (type==0) {
				for (int c=0;c<5;c++) /*printf(" %d",*/pt[c]=rand()%Ns[c]/*)*/;
				/*printf(" %d",*/v=rand()%401-200/*)*/;
				//cout << endl;
				st.update(v,pt[0],pt[1],pt[2],pt[3],pt[4]);
				A[pt[0]][pt[1]][pt[2]][pt[3]][pt[4]]=v;
			}
			else {
				for (int c=0;c<5;c++) {
					bounds[2*c]=rand()%Ns[c];
					bounds[2*c+1]=rand()%(Ns[c]-bounds[2*c])+bounds[2*c];
					//printf(" %d %d",bounds[2*c],bounds[2*c+1]);
				}
				//cout << endl;


				int r=0;
				for (int c=bounds[0];c<=bounds[1];c++)
					for (int d=bounds[2];d<=bounds[3];d++)
						for (int e=bounds[4];e<=bounds[5];e++)
							for (int f=bounds[6];f<=bounds[7];f++)
								for (int g=bounds[8];g<=bounds[9];g++)
									r+=A[c][d][e][f][g];
				if (r!=st.query(bounds[0],bounds[1],bounds[2],bounds[3],bounds[4],bounds[5],
							bounds[6],bounds[7],bounds[8],bounds[9]))
				{fprintf(stderr,"Failed on test %d, query %d\n",t,q);assert(0);}
			}
		}
	}

	return 0;
}
