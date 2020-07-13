/* 
 * Modified by Daniel for the world finals binder.
 *
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
#include <bits/stdc++.h>
using namespace std;

//listings:mdst
// Multidimensional vector required for multidimensional segment tree
template<class T, int D> struct Vec {typedef vector<typename Vec<T,D-1>::type> type;};
template<typename T> struct Vec<T,0> { typedef T type; };
template<typename T, int D> using MDV = typename Vec<T,D>::type;

// Multidimensional segment tree supporting point updates and ranged queries.
// The operation and an identity element must be provided by the template traits Op.
// Build initial data with st.build(A) where A is a D-dimensional vector of type T.
// Query ranges are closed hyperrectangles, updates are on single D-dimensional points.
// Complexity: O(2^D N1 N2..ND) to build, O(log(N1)log(N2)..log(ND)) to query.
template<typename T, typename Op, int D> struct SegmentTree {
	template<int _D> using ST = SegmentTree<T,Op,_D>;
	int N; vector<ST<D-1>> A;
	void build(const MDV<T,D>& data) {
		for (int c=0;c<N;c++) A[c+N].build(data[c]);
		for (int c=N-1;c;c--) A[c].merge(A[2*c],A[2*c+1]);
	}
	void merge(const ST<D>& L, const ST<D>& R) {
		for (int c=1;c<2*N;c++) A[c].merge(L.A[c],R.A[c]);
	}
	template<typename...Ts> void merge(const ST<D>& L, const ST<D>& R, int i, Ts...is) {
		for (i+=N;i;i/=2) A[i].merge(L.A[i],R.A[i],is...);
	}  // Create a segment tree with dimensions N1 * N2 * N3 * ...
  template<typename...Ts> SegmentTree(int N,Ts...Ns): N(N), A(2*N,ST<D-1>(Ns...)) {}
  // Set the value at (i1, i2, i3, ...) to v 
	template<typename...Ts> void update(const T& v, int i, Ts...is) {
		for (A[i+=N].update(v,is...);i/=2;) A[i].merge(A[2*i],A[2*i+1],is...);
	}  // Perform a ranged query on the range ([i1,j1] * [i2,j2] * ...)
	template<typename...Ts> T query(int i,int j,Ts...limits) {
		T r = Op::I;
		for (i+=N,j+=N;i<=j;i/=2,j/=2) {
			if (i&1) r = Op::op(r,A[i++].query(limits...));
			if (~j&1) r = Op::op(r,A[j--].query(limits...));
		}
		return r;
	}
};
template<typename T, typename Op> struct SegmentTree<T,Op,0> {
	typedef SegmentTree<T,Op,0> ST;  T a;
	SegmentTree() { a = Op::I; }
	void build(const T& data) {a=data;}
	void merge(const ST& L,const ST& R) { a = Op::op(L.a,R.a); }
	void update(const T& v) { a=v; }
	T query() { return a; }
};

// Example: Op for a multidimensional segment tree for ranged sums
struct RangeSumOp {
	static const int I = 0;
	static int op(const int& x, const int& y) { return x+y; }
};
//listings:/mdst

const int numtests=10,numqueries=10000;

int main() {
	srand(time(0));
	for (int t=0;t<numtests;t++) {
		//printf("Test: %d\n",t);
		//cerr << t << endl;
		vector<int> Ns(5);
		for (int c=0;c<5;c++) /*printf("%d\n",*/Ns[c]=rand()%25+1/*)*/;
		MDV<int,5> A{Ns[0],{Ns[1],{Ns[2],{Ns[3],vector<int>(Ns[4])}}}};

		for (int c=0;c<Ns[0];c++)
			for (int d=0;d<Ns[1];d++)
				for (int e=0;e<Ns[2];e++)
					for (int f=0;f<Ns[3];f++)
						for (int g=0;g<Ns[4];g++)
							/*printf("%d\n",*/A[c][d][e][f][g]=rand()%401-200/*)*/;

		SegmentTree<int,RangeSumOp,5> st(Ns[0],Ns[1],Ns[2],Ns[3],Ns[4]);
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
