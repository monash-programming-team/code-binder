#include <bits/stdc++.h>

using namespace std;

#define x first
#define y second

typedef vector<int> vi;

// Iterative Segment Tree with Lazy Updates
//
// Author      : Peter Whalan
// Date        : January 19, 2017
// Reliability : 5
// Tested On   : Custom Tests, SPOJ_LITE, Xtreme10_Safety, HackerCup2017_BigTop,
//	             CF_806E
//
// Performs range updates and queries on an array. Accepts a custom segment
// class, T, which contains both the value of the segment and any updates which
// need to be propagated to children. See below for examples of the T class.
// Intervals are 0-based and half-open.
//
// I is the identity. By default it default constructs the class T.
//
// op is used to find the answer for a segment given the answer on 2 subsegments
// that combine to make the segment.
//
// U is the type of an update. T must have a member named U which stores the
// update for a particular segment.
//
// queryRight finds the largest x such that b returns true for the range [l,x).
// x is never greater than N. queryLeft is similar but the range is [x,r). This
// method is analogous to binary search and so b must be decreasing i.e. if b is
// true for [0,x) and y<x then b is true for [0,y). The return value is a pair
// with x and the segment itself. Any extra arguments to this function are
// passed to b.
//
// Segfaults for N==0. Avoid using empty or invalid ranges although some may
// work. Sometimes a segment representing a prefix of the array is merged to the
// right of one representing a suffix due to wrap around. Ensure that this
// doesn't cause a runtime error. The resulting segment never gets queried so
// its value doesn't matter.
//
// It is a good idea to make identity updates fast since these will often occur.
//
// Complexity: O(N) to build, O(log N) to update and query. These complexities
// are for how many times op is called.

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

	void branchr(int i) {for (int j=__builtin_ctz(i);j;j--) prop(i>>j);}

	template<class...Ts>
	pair<int,T> queryRight(int l,Ts...args) {
		int r=l,w=1,p=0;
		t[0]=I;

		if (r<N) for (push(l+=N);r+2*w<=N && (t[1-p].op(t[p],A[l]),t[1-p].b(args...));)
			if (l&1) branchr(++l),r+=w,p^=1;
			else l/=2,w*=2;
		for (;w;l*=2,w/=2) if (r+w<=N && (prop(l/2),t[1-p].op(t[p],A[l]),t[1-p].b(args...)))
			l++,r+=w,p^=1;

		return {r,t[p]};
	}

	void branchl(int i) {for (int j=__builtin_ctz(i--);j;j--) prop(i>>j);}

	template<class...Ts>
	pair<int,T> queryLeft(int r,Ts...args) {
		int l=r,w=1,p=0;
		t[0]=I;

		for (push(r+=N-1);l>=2*w && (t[1-p].op(A[r],t[p]),t[1-p].b(args...));)
			if (~r&1) branchl(r--),l-=w,p^=1;
			else r/=2,w*=2;
		for (;w;r=2*r+1,w/=2) 
			if (l>=w && (prop(r/2),t[1-p].op(A[r],t[p]),t[1-p].b(args...))) r--,l-=w,p^=1;

		return {l,t[p]};
	}
};

namespace Testing {
	template <class T, class U>
		struct Slow {
			T I,t[4]; int N;
			vector<T> A;

			Slow(const vector<T>& data, T I=T()): I(I), N(data.size()), A(data) {}

			void update(int l, int r, U v) {
				for (int i=l;i<r;i++) A[i].us(v);
			}

			T query(int l, int r) {
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
		int a=INT_MAX,U=INT_MIN;

		void op(RangeMin& b,RangeMin& c) {a=min(b.a,c.a);}
		void us(int v) {
			if (v!=INT_MIN) a=U=v;
		}
		void NU() {U=INT_MIN;}
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
				SegmentTree<RangeMin,int> st(data);
				Slow<RangeMin,int> slow(data);
				for (int ct2=0;ct2<1000;ct2++) {
					int t=rand()%4;
					if (t==0) {
						int r=rand()%n+1,l=rand()%r,v=rnd();
						ss << l << ' ' << r << ' ' << v << '\n';
						slow.update(l,r,v);
						st.update(l,r,v);
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
			SegmentTree<RangeMin,int> st(data);
			Slow<RangeMin,int> slow(data);
			for (int i=0;i<=n;i++) for (int j=0;j<=n;j++)
				assert(st.queryRight(i,-j).x==slow.queryRight(i,-j).x);
		}
	}

	void testqueryLeft() {
		vector<RangeMin> data;
		for (int n=1;n<200;n++) {
			data.emplace_back();
			data.back().a=n;
			SegmentTree<RangeMin,int> st(data);
			Slow<RangeMin,int> slow(data);
			for (int i=0;i<=n;i++) for (int j=0;j<=n;j++)
				assert(st.queryLeft(i,j+1).x==slow.queryLeft(i,j+1).x);
		}
	}

	void solve() {
		cout << (sizeof(int)*8 - __builtin_clz(0)) << endl;
		testall();
		testqueryRight();
		testqueryLeft();

		/*SegmentTree<RangeMin, int> st(vector<RangeMin>(20));
		cout << st.query(0,20).a << endl;
		st.update(3,7,22);
		cout << st.query(0,3).a << ' ' << st.query(1,4).a << ' ' << st.query(3,7).a
			<< ' ' << st.query(6,8).a << ' ' << st.query(7,20).a << endl;
		st.update(4,5,0);
		cout << st.query(0,20).a << ' ' << st.query(4,5).a << ' ' << st.query(3,7).a << endl;
		*/
		//Fail
		//cout << st.query(0,0).a << endl;
	}
}

namespace SPOJ_LITE {
	struct STT {
		int a=0,w=0,U=0;
		void op(STT& b,STT& c) {
			a=b.a+c.a;
			w=b.w+c.w;
		}
		void us(int v) {
			if (v) {
				a=w-a;
				U^=1;
			}
		}
		void NU() {U=0;}
	};

	int N,M,t,i,j;
	void solve() {
		cin >> N >> M;
		STT base;
		base.w=1;
		SegmentTree<STT,int> st(vector<STT>(N,base));
		while (M--) {
			cin >> t >> i >> j;i--;
			if (t) cout << st.query(i,j).a << '\n';
			else st.update(i,j,1);
		}
	}
}

namespace Xtreme10_Safety {
	typedef long long ll;
	int N,M; string s;
	const int mod=1000000000+7,b=13;

	int chash[26][300001],p[300001];

	struct Update {
		int t=1,i,j,k,l=0;

		Update() {}
		Update(int a,int b,int c,int d,int e):
			t(a),i(b),j(c),k(d),l(e) {}
	};

	struct STL {
		int h[26],i=0,l=-1,r;
		Update U;

		void op(STL& b,STL& c) {
			if (b.l==-1) {*this=c; return;}
			if (c.l==-1) {*this=b; return;}
			i=0;
			l=b.l;
			r=c.r;
			if (c.r<c.l) return;
			for (int i=0,u=b.i,v=c.i;i<26;i++,u++,v++) {
				if (u==26) u=0;
				if (v==26) v=0;
				h[i]=((ll)b.h[u]*p[c.r-c.l] + c.h[v])%mod;
			}
		}

		void us(Update& u) {
			if (u.t==2) {
				i=u.l;
				for (int i=0;i<26;i++)
					h[i]=(chash[i][u.k+r-u.i]-(ll)p[r-l]*chash[i][u.k+l-u.i])%mod;

				U=u;
			}
			else if (u.t==3) {
				i=(i+u.l)%26;
				U.l=(U.l+u.l)%26;
				if (U.t==1) U.t=3;
			}
		}

		void NU() {
			U=Update();
		}
	};
	
	void solve() {
		cin >> s >> M; N=s.size();
		for (int i=0;i<26;i++) {
			chash[i][0]=0;
			for (int j=1;j<=N;j++) chash[i][j]=((ll)b*chash[i][j-1]+(s[j-1]-'a'+i)%26)%mod;
		}
		p[0]=1;
		for (int i=1;i<=N;i++) p[i]=(ll)p[i-1]*b%mod;

		vector<STL> data(N);
		for (int i=0;i<N;i++) {
			data[i].l=i;
			data[i].r=i+1;
			for (int j=0;j<26;j++) data[i].h[j]=(s[i]-'a'+j)%26;
		}

		SegmentTree<STL,Update> st(data);

		while (M--) {
			int t,i,j,k;
			cin >> t >> i >> j; i--;
			if (t!=3) cin >> k;

			if (t==1) {
				STL a,b;
				k--;
				a=st.query(i,j);
				b=st.query(k,k+j-i);
				if ((a.h[a.i]-b.h[b.i])%mod==0) cout << "Y\n";
				else cout << "N\n";
			}
			else if (t==2) {
				k--;
				st.update(i,j,Update(t,i,j,k,0));
			}
			else {
				st.update(i,j,Update(t,i,j,0,1));
			}
		}

	}
}

namespace HackerCup2017_BigTop {
#define x first
#define y second

	typedef long long ll;
	typedef pair<int,int> ii;

	struct ST {
		ll a=0;
		int l,r,U=INT_MIN;

		void op(ST& b,ST& c) {
			a=b.a+c.a;
			l=b.l;r=c.r;

			//w=b.w+c.w;
			//m=b.m;
			//M=c.M;
		}

		void us(ll v) {
			if (v==INT_MIN) return;
			a=(ll)(l+r-2*v)*(r-l);
			//a=w*v;
			U=v;
		}

		void NU() {U=INT_MIN;}
	};

	ll N,x1,ax,bx,cx,h1,ah,bh,ch;
	ll M;

	int X[800000],H[800000];

	int Mch=100000;

	void solve() {
		cout << fixed << setprecision(8);

		int T;
		cin >> T;
		for (int cas=1;cas<=T;cas++) {
			cin >> N >> x1 >> ax >> bx >> cx >> h1 >> ah >> bh >> ch;
			Mch=max(h1,ch);
			M=max(x1,cx)+4*Mch;
			vector<ST> data(M);
			for (int i=0;i<M;i++) {
				data[i].l=i;
				data[i].r=i+1;
			}
			SegmentTree<ST,int> st(data);

			set<ii> mins;

			X[0]=x1;
			H[0]=h1;
			for (int i=1;i<N;i++) {
				X[i]=((ax*X[i-1]+bx)%cx)+1;
				H[i]=((ah*H[i-1]+bh)%ch)+1;
			}

			double sum=0;

			for (int i=0;i<N;i++) {
				int v=X[i]-H[i]+2*Mch;
				int s=v,e=X[i]+H[i]+2*Mch;
				auto it=mins.upper_bound(ii(v,INT_MAX));
				if (it!=mins.begin()) {
					it--;
					s=max(s,it->y);
				}

				if (e>s) {

					mins.insert(ii(v,e));
					it=mins.lower_bound(ii(v,e));
					it++;
					while (it!=mins.end()) {
						if (it->y <= e) mins.erase(it++);
						else break;
					}

					st.update(s,e,v);
				}

				double x=st.query(0,M).a;
				sum+=x;

			}

			cout << "Case #" << cas << ": " << sum/4 << '\n';
		}
	}
#undef x
#undef y
}

namespace CF_806E {
	const int inf=INT_MAX/10;
	typedef pair<int,int> ii;
#define x first
#define y second
	struct ST {
		int sp=inf,sn=-inf,U=0;

		void op(ST &b,ST &c) {
			sp=min(b.sp,c.sp);
			sn=max(b.sn,c.sn);
		}
		void us(int v) {
			sp+=v;
			sn+=v;
			U+=v;
		}
		void NU() {U=0;}
		bool b(bool plus) {
			if (plus) return sn<-1;
			return sp>1;
		}

		int val() {
			if (sp<-sn) return sp;
			return sn;
		}
	};

	int n,diff[500000],ord[500000],ans[500000];
	ii a[500000];

	void solve() {
		ios::sync_with_stdio(0);
		cin.tie(0);

		int v=0;
		cin >> n;
		for (int i=0;i<n;i++) {cin >> a[i].x;a[i].y=i;}
		sort(a,a+n);

		for (int i=0;i<n;i++) {
			diff[i]=a[i].x-v;
			v+=(diff[i]>0)-(diff[i]<0);
		}

		vector<ST> data(n);
		for (int i=0;i<n;i++) {
			if (diff[i]>=0) data[i].sp=diff[i];
			if (diff[i]<=0) data[i].sn=diff[i];
		}
		SegmentTree<ST,int> st(data);

		for (int i=0;i<n;i++) ord[a[i].y]=i;

		for (int i=n-1;i>=0;i--) {
			ans[i]=v;
			ST seg=st.query(ord[i],ord[i]+1);
			int d;
			if (seg.sp<-seg.sn) d=seg.sp;
			else d=seg.sn;
			int sgn=(d>0)-(d<0);
			st.push(ord[i]+n);
			st.A[ord[i]+n].sp=inf;
			st.A[ord[i]+n].sn=-inf;
			for (int j=ord[i]+n;j/=2;) st.A[j].op(st.A[2*j],st.A[2*j+1]);
			if (sgn) {
				int u=st.queryRight(ord[i]+1,sgn==1).x;
				if (ord[i]+1<u) st.update(ord[i]+1,u,sgn);
				if (u<n) {
					st.push(u+n);
					st.A[u+n].us(sgn);
					st.A[u+n].sp=st.A[u+n].sn=st.A[u+n].val();
					if (st.A[u+n].sp<0) st.A[u+n].sp=inf;
					if (st.A[u+n].sn>0) st.A[u+n].sn=-inf;
					for (int j=u+n;j/=2;) st.A[j].op(st.A[2*j],st.A[2*j+1]);
				}
				else v-=sgn;
			}
		}

		for (int i=0;i<n;i++) cout << ans[i] << '\n';
	}
#undef x
#undef y
}

int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);

	//Testing::solve();
	//SPOJ_LITE::solve();
	//Xtreme10_Safety::solve();
	//HackerCup2017_BigTop::solve();
	CF_806E::solve();
}

