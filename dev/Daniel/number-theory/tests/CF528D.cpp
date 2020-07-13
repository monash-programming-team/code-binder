// Codeforces 528D - Fuzzy Search
#include <bits/stdc++.h>
using namespace std;;

#define debug(x) cerr << #x << " = " << x << endl;

typedef long long int ll;
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef ll big;

ll gcd(ll a, ll b, ll& x, ll& y) {
	if (b == 0) { y = 0; x = (a < 0) ? -1 : 1; return (a < 0) ? -a : a; }
	else { ll g = gcd(b, a%b, y, x); y -= a/b*x; return g; }
}

ll inv(ll a, ll m) { ll x, y; gcd(m,a,x,y); return ((y % m) + m) % m; }

template<typename T> struct convolution {
	const T m, r, ord;
	T mult(T x, T y) { return big(x) * y % m; }
	
	void ntt(vector<T> & a, int invert = 0) {
		int n = (int) a.size();
    T ninv = inv(n, m), rinv = inv(r, m);
		for (int i=1, j=0; i<n; ++i) {
			int bit = n >> 1;
			for (; j>=bit; bit>>=1)	j -= bit;
			j += bit;
			if (i < j) swap (a[i], a[j]);
		}
		for (int len=2; len<=n; len<<=1) {
			T wlen = invert ? rinv : r;
			for (int i=len; i<ord; i<<=1) wlen = mult(wlen, wlen);
			for (int i=0; i<n; i+=len) {
				T w = 1;
				for (int j=0; j<len/2; ++j) {
					T u = a[i+j],  v = mult(a[i+j+len/2], w);
					a[i+j] = u + v < m ? u + v : u + v - m;
					a[i+j+len/2] = u - v >= 0 ? u - v : u - v + m;
					w = mult(w, wlen);
				}
			}
		}
		if (invert) for (int i=0; i<n; ++i)	a[i] = mult(a[i], ninv);
	}
  
	vector<T> multiply(vector<T>& a, vector<T>& b) {
		vector<T> fa(a.begin(), a.end()), fb(b.begin(), b.end());
		int n = 1;  while (n < 2 * (int)max(a.size(), b.size())) n*=2;
		fa.resize(n), fb.resize(n);
		ntt(fa), ntt(fb);
		for(int i=0;i<n;i++) fa[i] = mult(fa[i], fb[i]);
		ntt(fa, 1);
		fa.resize(n);
		return fa;
	}
};

struct FenwickTree {
	int N;  vi A;
	FenwickTree(int n): N(n+1), A(N) {}
	void adjust(int b,int v) { for (;b;b-=b&-b) A[b]+=v; }
	void adjust(int a,int b,int v) { adjust(b,v); adjust(a,-v); }
	int query(int i) { int r=0; for (i++;i<N;i+=i&-i) r+=A[i]; return r;	}
};

const string alph = "ACGT";
int n, m, k;
string S, T;

int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);
	cin >> n >> m >> k >> S >> T;
  
  vvi S_char(4, vi(n)), T_char(4, vi(m));
  
  // Build indicator sequences for S
  vector<FenwickTree> indicator(4, FenwickTree(n));
  for (int i=0; i<n; i++)
    indicator[alph.find(S[i])].adjust(max(i-k, 0), min(i+k+1,n), 1);
  
  for (int c=0; c<4; c++) for (int i=0; i<n; i++) 
    S_char[c][i] = (indicator[c].query(i) > 0);

  // Build indicator sequences for T
  for (int j=0; j<m; j++) T_char[alph.find(T[j])][m-1-j] = 1;
  
  // Convolutions
  convolution<int> conv{7340033, 5, 1 << 20};
  vvi matches(4);
  for (int c=0; c<4; c++) matches[c] = conv.multiply(S_char[c], T_char[c]);
  
  // Find matches
  int ans = 0;
  for (int i=m-1; i<n; i++) {
    int correct = 0;
    for (int c=0; c<4; c++) correct += matches[c][i];
    if (correct == m) ans++;
  }
  cout << ans << endl;
}
