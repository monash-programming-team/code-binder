// ACM ICPC 2013-2014. NEERC. Eastern Subregional Contest
// Problem G. Cipher Message 3
//
// Verdict: AC
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

// Fast convolution using NTT
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

// KMP Prefix function
vi prefix(const vi &pat){
	int m = (int)pat.size(), j = -1;
	vi P(m + 1, -1);
	for (int i = 0; i < m; i++){
		while(j >= 0 && pat[i] != pat[j]) j = P[j];
		j++;
		P[i+1] = j;
	}
	return P;
}
// KMP Pattern finding
vi match(const vi &txt, const vi &pat, const vi &P){
	int n = (int)txt.size(), m = (int)pat.size(), j = 0;
	vi idx;
	for (int i = 0; i < n; i++){
		while (j >= 0 && txt[i] != pat[j]) j = P[j];
		j++;
		if (j == m){
			idx.push_back(i - j + 1);
			j = P[j];
		}
	}
	return idx;
}
int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);
	
	int n, m; 
	cin >> n >> m;
	vi t_msb(n), t_lsb(n), p_msb(m), p_lsb(m);
  
  // Read the string
	for (int i = 0; i < n; i++){
		string str; cin >> str;
		t_lsb[i] = (int)str.back() - '0';
		t_msb[i] = 0;
		for (int j = 0; j < 7; j++){
			t_msb[i] *= 2;
			t_msb[i] += (int)str[j] - '0';
		}
	}
  // Read the pattern
	for (int i = 0; i < m; i++){
		string str; cin >> str;
		p_lsb[i] = (int)str.back() - '0';
		p_msb[i] = 0;
		for (int j = 0; j < 7; j++){
			p_msb[i] *= 2;
			p_msb[i] += (int)str[j] - '0';
		}
	}
  // Find matches in the most significant bits
	vi P = prefix(p_msb);
	vi indx = match(t_msb, p_msb, P);
  
  // Reverse the pattern for the convolution
  reverse(p_lsb.begin(), p_lsb.end());
  // Indicator functions for zero and one bits
  vi p_ones = p_lsb, p_zeros = p_lsb; 
  for (auto& x : p_zeros) x ^= 1;
  vi t_ones = t_lsb, t_zeros = t_lsb; 
  for (auto& x : t_zeros) x ^= 1;
  
  // Compute the convolution with the string
  convolution<int> conv{7340033LL, 5, 1 << 20};
  vi ones_match = conv.multiply(t_ones, p_ones);
  vi zeros_match = conv.multiply(t_zeros, p_zeros);
  
  // Find the best place to put the pattern  
	int res = INT_MAX, id = -1;
	for (int x : indx){
    int correct = ones_match[x + m - 1] + zeros_match[x + m - 1];
    int diff = m - correct;
    if (diff < res) res = diff, id = x;
	}
  
  // Output the answer
	if (res == INT_MAX) cout << "No" << endl;
	else {
		cout << "Yes" << endl;
		cout << res << ' ' << id + 1 << endl;
	}
}
