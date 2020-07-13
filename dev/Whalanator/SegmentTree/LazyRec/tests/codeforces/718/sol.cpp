#include <bits/stdc++.h>

#define x first
#define y second

using namespace std;

typedef long long ll;
typedef long double ld;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> ii;
typedef pair<ll,ll> pll;

//template<class T,class U> bool smin(T &a,U b) { return a > b ? (a = b, true) : false; }
//template<class T,class U> bool smax(T &a,U B) { return a < b ? (a = b, true) : false; }
//
const ll mod=1000000000+7;


//A Type, U Type
/*template <class AT, class UT>
struct SegmentTree {
	AT I,res[2],*rc,*ro; UT NU,v;//I is identity, NU is no update constant
	int N,i,j;
	vector<AT> A; vector<UT> U; vi L,R;

	//Implement these in specialisation or can inline if you only want one type of segment tree
	//Then just uncomment build call
	virtual void op(AT& r, AT& i, AT& j)=0;
	virtual void us(int p, UT& v)=0;
	//virtual AT op(const AT& i, const AT& j)=0;
	//virtual void us(int p, const UT& v)=0;

	//Usually don't need this.
	//Use other constructor with identity array
	SegmentTree(int N, AT I, UT NU):
		I(I), NU(NU), N(N), A(4*N,I), U(4*N,NU), L(4*N), R(4*N) {build(1,0,N-1);}
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
			/*build(1,0,N-1,data);*//*
			res[0]=res[1]=I;
		}

	void build(int p, int i, int j,const vector<AT>& data) {
		L[p]=i;R[p]=j;
		if (i==j) {A[p]=data[i];return;}
		build(2*p,i,(i+j)/2,data);
		build(2*p+1,(i+j)/2+1,j,data);
		op(A[p],A[2*p],A[2*p+1]);
	}

	void update(int p) {
		if (i<=L[p] && R[p]<=j) {return us(p,v);}// return;}
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
		//res[0]=I;
		for (int i=0;i<2;i++) for (int j=0;j<2;j++) res[0][i][j]=0;
		rc=res;ro=res+1;
		query(1);
		return *rc;
	}
};
*/
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

typedef vector<vector<ll>> Mat;

Mat operator*(const Mat& a,const Mat& b) {
	Mat r(2,{0,0});
	for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
			for (int k=0;k<2;k++)
				r[i][j]=(r[i][j]+a[i][k]*b[k][j])%mod;
		}
	return r;
}

Mat operator+(const Mat& a,const Mat& b) {
	Mat r(2,{0,0});
	for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
		r[i][j]=(a[i][j]+b[i][j]);
		if (r[i][j]>=mod) r[i][j]-=mod;
	}
	return r;
}

/*struct Mat {
	ll m[2][2];

	Mat operator* (Mat b) {
		Mat r;
		for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
			r.m[i][j]=0;
			for (int k=0;k<2;k++)
				r.m[i][j]=(r.m[i][j]+m[i][k]*b.m[k][j]%mod)%mod;
		}
		return r;
	}

	Mat operator+ (Mat b) {
		Mat r;
		for (int i=0;i<2;i++) for (int j=0;j<2;j++)
			r.m[i][j]=(m[i][j]+b.m[i][j])%mod;
		return r;
	}

	Mat operator== (Mat b) {
		for (int i=0;i<2;i++) for (int j=0;j<2;j++)
			if (m[i][j]!=b.m[i][j]) return 0;
		return 1;
	}
};*/

Mat addid() {
	return Mat(2,{0,0});
}

Mat id() {
	Mat r(2,{0,0});
	for (int i=0;i<2;i++) for (int j=0;j<2;j++) r[i][j]=(i==j);
	return r;
}

const Mat I=id();
Mat temp=id();

struct Pup : SegmentTree<Mat, Mat> {
	Pup(const vector<Mat>& data): SegmentTree(data,addid(),id()) {build(1,0,N-1,data);}
	void op(Mat &r,Mat& i,Mat& j) {
		for (int c=0;c<2;c++) for (int d=0;d<2;d++) {
			r[c][d]=(i[c][d]+j[c][d]);
			if (r[c][d]>=mod) r[c][d]-=mod;
		}
	}
	void us(int p,Mat& v) {
		for (int c=0;c<2;c++) for (int d=0;d<2;d++) {
			temp[c][d]=0;
			for (int e=0;e<2;e++) temp[c][d]=(temp[c][d]+A[p][c][e]*v[e][d])%mod;
		}
		for (int c=0;c<2;c++) for (int d=0;d<2;d++) A[p][c][d]=temp[c][d];
		for (int c=0;c<2;c++) for (int d=0;d<2;d++) {
			temp[c][d]=0;
			for (int e=0;e<2;e++) temp[c][d]=(temp[c][d]+U[p][c][e]*v[e][d])%mod;
		}
		for (int c=0;c<2;c++) for (int d=0;d<2;d++) U[p][c][d]=temp[c][d];
	}
};

Mat fi() {
	Mat r=I;
	for (int i=0;i<2;i++) for (int j=0;j<2;j++) r[i][j]=(i!=j || (i==j && j==0));
	return r;
}

Mat p1;

Mat *r1,*r2,r[2];

void fpo2(int e) {
	if (e==0) return;
	fpo2(e/2);
	for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
		(*r2)[i][j]=0;
		for (int k=0;k<2;k++) (*r2)[i][j]=((*r2)[i][j]+(*r1)[i][k]*(*r1)[k][j])%mod;
	}
	if (e%2) for (int i=0;i<2;i++) for (int j=0;j<2;j++) {
		(*r1)[i][j]=0;
		for (int k=0;k<2;k++) (*r1)[i][j]=((*r1)[i][j]+(*r2)[i][k]*p1[k][j])%mod;
	}
	else swap(r1,r2);
}

Mat fpo(int e) {
	r[0]=r[1]=id();
	r1=r;r2=r+1;
	fpo2(e);
	return *r1;
}



int n,m;
vector<Mat> init;

int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);
	p1=fi();

	cin >> n >> m;
	init.resize(n);

	for (int i=0;i<n;i++) {
		int a;cin >> a;a--;
		init[i]=fpo(a);
	}

	Pup st(init);
	for (int q=0;q<m;q++) {
		int t,l,r,x;
		cin >> t >> l >> r;
		l--;r--;
		if (t==1) {
			cin >> x;
			st.update(l,r,fpo(x));
		}
		else cout << st.query(l,r)[0][0] << '\n';
	}

}
