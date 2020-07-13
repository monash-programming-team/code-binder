#include <bits/stdc++.h>

using namespace std;

#define x first
#define y second

typedef vector<int> vi;

// Iterative Segment Tree
//
// Author      : Peter Whalan
// Date        : July 6, 2016
// Reliability : 3
// Tested On   : CF_750E, HackerCup2017_Fighting, CF_759C
//
// Performs point updates and range queries on array. Accepts a custom segment
// class, T, for custom operations. e.g. RangeMin. T stores the value of the
// segment and has the methods required by the segment tree. op is used to find
// the answer for a segment given the answer on 2 subsegments that combine to
// make the segment. The identity can be passed in as an argument to the
// constructor. Otherwise it is default contructed.
//
// queryRight finds the largest x such that b returns true for the range [l,x).
// x is never greater than N. queryLeft is similar but the range is [x,r). This
// method is analogous to binary search and so b must be decreasing i.e. if b is
// true for [0,x) and y<x then b is true for [0,y). The return value is a pair
// with x and the segment itself. Any extra arguments to this function are
// passed to b.
//
// Second constructor segfaults for N==0. Change loop condition to c>0 to not
// segfault. Sometimes a segment representing a prefix of the array is merged to the
// right of one representing a suffix due to wrap around. Ensure that this
// doesn't cause a runtime error. The resulting segment never gets queried so
// its value doesn't matter.
//
// Complexity: O(N) to build, O(log N) to update and query. These complexities
// are for how many times op is called.

template<class T>
struct SegmentTree {
	T I,t[4]; int N;
	vector<T> A;

	SegmentTree(int N, T I=T()): I(I), N(N), A(2*N,I) {}
	SegmentTree(const vector<T>& data, T I=T()): I(I), N(data.size()), A(2*N,I)
	{
		copy(data.begin(),data.end(),A.begin()+N);
		for (int i=N-1;i;i--) A[i].op(A[2*i],A[2*i+1]);
	}

	void update(int i,T v) {
		for (A[i+=N]=v;i/=2;) A[i].op(A[2*i],A[2*i+1]);
	}

	T query(int l,int r) {
		t[0]=t[2]=I;
		int i=0,j=2;
		for (l+=N,r+=N;l<r;l/=2,r/=2) {
			if (l&1) t[i^1].op(t[i],A[l++]), i^=1;
			if (r&1) t[j^1].op(A[--r],t[j]), j^=1;
		}
		t[i^1].op(t[i],t[j]);
		return t[i^1];
	}

	template<class...Ts>
	pair<int,T> queryRight(int l,Ts...args) {
		int r=l,w=1,p=0;
		t[0]=I;

		for (l+=N;r+2*w<=N && (t[1-p].op(t[p],A[l]),t[1-p].b(args...));)
			if (l&1) l++,r+=w,p^=1;
			else l/=2,w*=2;
		for (;w;l*=2,w/=2) if (r+w<=N && (t[1-p].op(t[p],A[l]),t[1-p].b(args...)))
			l++,r+=w,p^=1;

		return {r,t[p]};
	}

	template<class...Ts>
	pair<int,T> queryLeft(int r,Ts...args) {
		int l=r,w=1,p=0;
		t[0]=I;

		for (r+=N;l>=2*w && (t[1-p].op(A[r-1],t[p]),t[1-p].b(args...));)
			if (r&1) r--,l-=w,p^=1;
			else r/=2,w*=2;
		for (;w;r*=2,w/=2) if (l>=w && (t[1-p].op(A[r-1],t[p]),t[1-p].b(args...)))
			r--,l-=w,p^=1;

		return {l,t[p]};
	}
};

namespace Testing {
	template<class T>
		struct Slow {
			T I,t[4]; int N;
			vector<T> A;

			Slow(int N, T I=T()): I(I), N(N), A(2*N,I) {}
			Slow(const vector<T>& data, T I=T()): I(I), N(data.size()), A(data) {}

			void update(int i,T v) {
				A[i]=v;
			}

			T query(int l,int r) {
				t[0]=I;
				int p=0;
				for (int i=l;i<r;i++) t[p^1].op(t[p],A[i]), p^=1;
				return t[p];
			}

			template<class...Ts>
				pair<int,T> queryRight(int l,Ts...args) {
					int r=l,p=0;
					t[0]=I;
					while (r<N && (t[1-p].op(t[p],A[r]),t[1-p].b(args...))) p^=1,r++;
					return {r,t[p]};
				}

			template<class...Ts>
				pair<int,T> queryLeft(int r,Ts...args) {
					int l=r,p=0;
					t[0]=I;
					while (l && (t[1-p].op(A[l-1],t[p]),t[1-p].b(args...))) p^=1,l--;
					return {l,t[p]};
				}
		};

	struct RangeMin {
		int a=INT_MAX;
		void op(RangeMin& b,RangeMin& c) {a=min(b.a,c.a);}
		bool b(int v) {return a>=v;}
	};

	int rnd() {return rand()%20;}

	void testall() {
		vector<RangeMin> data;
		for (int n=1;n<200;n++) {
			data.emplace_back();
			for (int ct=0;ct<1000;ct++) {
				stringstream ss;
				for (int j=0;j<n;j++) data[j].a=rnd();
				for (auto d:data) ss << d.a << ' ';
				ss << '\n';
				SegmentTree<RangeMin> st(data);
				Slow<RangeMin> slow(data);
				for (int ct2=0;ct2<1000;ct2++) {
					int t=rand()%4;
					if (t==0) {
						int i=rand()%n,v=rnd();
						ss << i << ' ' << v << '\n';
						slow.update(i,{v});
						st.update(i,{v});
					}
					else if (t==1) {
						int r=rand()%n+1,l=rand()%r;
						int slowa=slow.query(l,r).a;
						int sta=st.query(l,r).a;
						ss << "t=1 " << l << ' ' << r << '\n';
						if (slowa!=sta) {
							cerr << slowa << ' ' << sta << '\n';
							cerr << "It is wrong\n";
							string s;
							while (getline(ss,s)) cerr << s << '\n';
							assert(0);
						}
					}
					else if (t==2) {
						int l=rand()%(n+1),v=rnd();
						ss << "t==2 " << l << ' ' << v << '\n';
						auto slowa=slow.queryRight(l,v);
						auto sta=st.queryRight(l,v);
						if (slowa.x!=sta.x) {
							cerr << slowa.x << ' ' << sta.x << '\n';
							cerr << "It is wrong\n";
							string s;
							while (getline(ss,s)) cerr << s << '\n';
							assert(0);
						}
					}
					else {
						int r=rand()%(n+1),v=rnd();
						ss << "t==3 " << r << ' ' << v << '\n';
						auto slowa=slow.queryLeft(r,v);
						auto sta=st.queryLeft(r,v);
						if (slowa.x!=sta.x) {
							cerr << slowa.x << ' ' << sta.x << '\n';
							cerr << "It is wrong\n";
							string s;
							while (getline(ss,s)) cerr << s << '\n';
							assert(0);
						}
					}
				}
			}
		}
	}

	void testqueryRight() {
		vector<RangeMin> data;
		for (int n=1;n<200;n++) {
			data.emplace_back();
			data.back().a=-n;
			SegmentTree<RangeMin> st(data);
			Slow<RangeMin> slow(data);
			for (int i=0;i<=n;i++) for (int j=0;j<=n;j++)
				assert(st.queryRight(i,-j).x==slow.queryRight(i,-j).x);
		}
	}

	void testqueryLeft() {
		vector<RangeMin> data;
		for (int n=1;n<200;n++) {
			data.emplace_back();
			data.back().a=n;
			SegmentTree<RangeMin> st(data);
			Slow<RangeMin> slow(data);
			for (int i=0;i<=n;i++) for (int j=0;j<=n;j++)
				assert(st.queryLeft(i,j+1).x==slow.queryLeft(i,j+1).x);
		}
	}

	void solve() {
		testall();
		testqueryRight();
		testqueryLeft();

		SegmentTree<RangeMin> st(10);
		cout << st.query(0,6).a << endl;
		RangeMin a;
		a.a=8;
		st.update(7,a);
		cout << st.query(0,10).a << endl;
	}
}

namespace CF_750E {
	template<class T,class U> bool smin(T& a,U b) {return a>b?(a=b,1):0;}

	string bad="2016";
	string good="2017";

	struct TT {
		int a[5][5];
		void op(TT& i, TT& j) {
			for (int c=0;c<5;c++) for (int d=0;d<5;d++) a[c][d]=INT_MAX/2;
			for (int c=0;c<5;c++) for (int d=c;d<5;d++) {
				for (int e=c;e<=d;e++) smin(a[c][d],i.a[c][e]+j.a[e][d]);
			}
		}
		TT() {
			for (int c=0;c<5;c++) for (int d=0;d<5;d++) {
				a[c][d]=(c==d?0:INT_MAX/2);
			}
		}
	};

	int n,q;
	string s;

	int sum[200001];
	int l6[200001];
	void solve() {
		cin >> n >> q >> s;

		vector<TT> data(n);
		for (int i=0;i<n;i++) {
			for (int c=0;c<4;c++) {
				if (bad[c]==s[i]) {
					data[i].a[c][c]=1;
					data[i].a[c][c+1]=0;
				}
				else data[i].a[c][c]=0;
			}
			data[i].a[4][4]=0;
		}

		SegmentTree<TT> st(data);

		for (int i=0;i<n;i++) sum[i+1]=sum[i]+(s[i]=='6');

		l6[0]=-1;
		for (int i=0;i<n;i++) {
			l6[i+1]=l6[i];
			if (s[i]=='7') l6[i+1]=i;
		}

		for (int qc=0;qc<q;qc++) {
			int a,b;
			cin >> a >> b;
			a--;b--;
			if (l6[b+1]<a) {
				cout << -1 << '\n';
				continue;
			}
			int c1=st.query(a,l6[b+1]+1).a[0][3];
			if (c1==INT_MAX/2) {
				cout << -1 << '\n';
				continue;
			}
			cout << c1+sum[b+1]-sum[l6[b+1]] << '\n';
		}

	}
}

namespace HackerCup2017_Fighting {
	typedef long long ll;
	const int mod=1000000000+7;

	int addm(int& a,int b) {return (a+=b)<mod?a:a-=mod;}
	struct ST {
		int a[2][2]={};

		void op(ST& b,ST& c) {
			for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
				a[i][j]=((ll)b.a[i][0]*c.a[0][j]+(ll)b.a[i][1]*c.a[1][j])%mod;
			}
		}
	};

	int N,M;
	int w1,aw,bw,d1,ad,bd,s1,as,bs;
	int W[800000],Z[800000],S[800000];

	void solve() {
		int T;
		cin >> T;
		for (int cas=1;cas<=T;cas++) {
			cin >> N >> M >> w1 >> aw >> bw >> d1 >> ad >> bd >> s1 >> as >> bs;
			W[0]=w1;
			Z[0]=max(1,min(N,w1+d1-1));
			S[0]=s1;
			int d2;
			for (int i=1;i<M;i++) {
				W[i]=(((ll)aw*W[i-1]+bw)%N)+1;
				d2=((ll)ad*d1+bd)%3;
				swap(d1,d2);
				Z[i]=max(1,min(N,W[i]+d1-1));
				S[i]=(((ll)as*S[i-1]+bs)%1000000000)+1;
			}

			for (int i=0;i<M;i++) {W[i]--;Z[i]--;}

			vector<ST> data(N);
			for (ST &i: data) i.a[0][0]=1;
			ST I;
			I.a[0][0]=I.a[1][1]=1;
			SegmentTree<ST> st(data,I);

			int ans=0;
			for (int i=0;i<M;i++) {
				ST up=st.A[st.N+Z[i]];
				if (W[i]<Z[i]) addm(up.a[1][0],S[i]);
				else if (W[i]==Z[i]) addm(up.a[0][0],S[i]);
				else addm(up.a[0][1],S[i]);
				st.update(Z[i],up);

				addm(ans,st.query(0,N).a[0][0]);

			}

			cout << "Case #" << cas << ": " << ans << '\n';
		}
	}
}

namespace CF_759C {
#define x first
#define y second
	
	typedef pair<int,int> ii;

	set<ii> st;
	int m;
	int psh[100000];

	struct ST {
		int s=0,M=0;

		ST(int a,int b):s(a),M(b) {}
		ST(): s(0), M(0) {}

		void op(ST& b,ST& c) {
			s=b.s+c.s;
			M=max(c.M,b.M+c.s);
		}

		bool b() {return M<1;}
	};

	void solve() {
		cin >> m;

		SegmentTree<ST> st(m);

		for (int ct=0;ct<m;ct++) {
			int i,j,k;
			cin >> i >> j;
			if (j) cin >> k;
			i--;
			if (j) psh[i]=k;

			if (j) st.update(i,ST(1,1));
			else st.update(i,ST(-1,0));

			int pos=st.queryLeft(m).x-1;

			if (pos<0) cout << "-1\n";
			else cout << psh[pos] << '\n';
		}
	}

#undef x
#undef y
}
int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);

	Testing::solve();
	//CF_750E::solve();
	//HackerCup2017_Fighting::solve();
	//CF_759C::solve();
}
