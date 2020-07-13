// Number theoretic transform for fast convolution
//
// Author: Daniel Anderson, based on e-maxx.ru
// Date: 10-12-2016
// Reliability: 5
// Tested on: UVA12879, CF632E, CF300D, NEERC-EAST_13G, CF528D
//
// Note on big mods: Using any of the mods bigger than 1e9 will require
// __int128 for the intermediate calculations. Just follow the typedef
// below (alias __int128 as "big"). This will however slow the code down
// by roughly a factor of five, so only use __int128 if its actually
// necessary. An alternative way to recover large numbers is to solve
// the problem for multiple small mods then use the Chinese remainder
// algorithm, as this will be much faster than solving NTT with a big mod.
//
#include<bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef ll big;  // use __int128 when required and available

ll gcd(ll a, ll b, ll& x, ll& y) {
	if (b == 0) { y = 0; x = (a < 0) ? -1 : 1; return (a < 0) ? -a : a; }
	else { ll g = gcd(b, a%b, y, x); y -= a/b*x; return g; }
}

ll inv(ll a, ll m) { ll x, y; gcd(m,a,x,y); return ((y % m) + m) % m; }

//listings:ntt
// Integer convolution mod m using number theoretic transform.
// m = modulo, r = a primitive root, ord = order of the root
// (Must be a power of two). The length of the given input
// vectors must not exceed n = ord.   Complexity: O(n log(n))
//
// Usable coefficients::
//   m                    | r           | ord       | __int128 required
//------------------------|-------------|-----------|--------------------
//  7340033               | 5           | 1 << 20   | No
//  469762049             | 13          | 1 << 25   | No
//  998244353             | 31          | 1 << 23   | No
//  1107296257            | 8           | 1 << 24   | No
//  10000093151233        | 366508      | 1 << 26   | Yes
//  1000000523862017      | 2127080     | 1 << 26   | Yes
//  1000000000949747713   | 465958852   | 1 << 26   | Yes
//
// In general, you may use mod = c * 2^k + 1 which has a primitive
// root of order 2^k, then use number theory to find a generator.
template<typename T> struct convolution {
	const T m, r, ord;
	T mult(T x, T y) { return big(x) * y % m; }
	void ntt(vector<T> & a, int invert = 0) {
		int n = (int)a.size(); T ninv = inv(n, m), rinv = inv(r, m);  // Modular inverses
		for (int i=1, j=0; i<n; ++i) {
			int bit = n >> 1;	for (; j>=bit; bit>>=1)	j -= bit;
			j += bit;	if (i < j) swap (a[i], a[j]);
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
	// Compute the convolution a * b -- Complexity: O(n log(n))
	vector<T> multiply(vector<T>& a, vector<T>& b) {
		vector<T> fa(a.begin(), a.end()), fb(b.begin(), b.end());
		int n = 1;  while (n < 2 * (int)max(a.size(), b.size())) n*=2;
		fa.resize(n), fb.resize(n);	ntt(fa), ntt(fb);
		for(int i=0;i<n;i++) fa[i] = mult(fa[i], fb[i]);
		ntt(fa, 1);	fa.resize(n);
		return fa;
	}
};
//listings:/ntt

// Compute the convolution a^k -- Complexity: O(n log(n) log(k))
template<typename T> vector<T> power(convolution<T> c, vector<T> a, int k) {
  vector<T> res(1,1);
  for (; k; k /= 2) {
    if (k & 1) res = c.multiply(res, a);
    a = c.multiply(a, a);
  }
  return res;
}

// ------------------------------------------------------------------
//						               TEST PROBLEMS
// ------------------------------------------------------------------

// Verdict: AC
void solve_UVA12879() {
  convolution<ll> conv{1000000000949747713LL, 465958852, 1 << 26};
  int N, k, M, d;
  while (cin >> N) {
    vector<ll> dists;
    for (int i=0; i<N; i++) {
      cin >> k;
      if ((int)dists.size() < k + 1) dists.resize(k+1);
      dists[k] = 1;
    }
    auto bot = conv.multiply(dists, dists);
    cin >> M;  int ans = 0;
    for (int i=0; i<M; i++) {
      cin >> d;
      ans += (((int)bot.size() > d && bot[d] != 0)
        || ((int)dists.size() > d && dists[d] != 0));
    }
    cout << ans << '\n';
  }
}

// Verdict: AC
void solve_CF632E() {
  convolution<ll> conv{1107296257LL, 8, 1 << 24};
  int n, k; cin >> n >> k;
  vi a(n); for (int i=0; i<n; i++) cin >> a[i];
  vector<ll> ans(1001); for (auto x: a) ans[x] = 1;
  auto res = power(conv, ans, k);
  bool space = false;
  for (int i=0; i<(int)res.size(); i++) {
    if (res[i]) {
      if (space) cout << ' ';
      cout << i;  space = true;
    }
  }
  cout << endl;
}

// Verdict: AC
void solve_CF300D() {
  // Pre-compute DP values
  vector<vector<ll>> dp(31, vector<ll>(1001));
  vector<ll> d;
  convolution<ll> conv{7340033LL, 5, 1 << 20};
  for (int i=0; i<=29; i++) {
    dp[i][0] = 1;
    d = conv.multiply(dp[i], dp[i]);
    d.resize(1001);
    dp[i+1] = conv.multiply(d, d);
    dp[i+1].insert(dp[i+1].begin(),0);
    dp[i+1].resize(1001);
  }
  // Answer queries
  int tc; cin >> tc;
  for (int i = 0; i < tc; i++) {
    ll n, k; cin >> n >> k;
    int it = 0;
    while (n > 1 && (n&1)) it++, n /= 2;
    cout << dp[it][k] << '\n';
  }
}

// ------------------------------------------------------------------
//						TEST PARAMETER DATA AND CONVOLUTION
// ------------------------------------------------------------------

// Parameters for NTT
vector<vector<ll>> params = {
  {7340033LL, 5, 1 << 20},                     // good
  {469762049LL, 13, 1 << 25},                  // good
  {998244353LL, 31, 1 << 23},                  // good
  {1107296257LL, 8, 1 << 24},                  // good
  {10000093151233LL, 366508, 1 << 26},         // good (requires int128)
  {1000000523862017LL, 2127080, 1 << 26},      // good (requires int128)
  {1000000000949747713LL, 465958852, 1 << 26}  // good (required int128)
};

// Expmod
ll expmod(big a, big b, big m) {
	big res=1;  a %= m;
	for(; b; b /= 2) { if (b&1) res=res*a%m;  a=a*a%m; }
	return res;
}

// Prime test
vi val = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}; // n <= 2^64
bool is_prime(ll n) {
	if (n < 2) return false;
	ll s = __builtin_ctzll(n-1), d = (n-1) >> s;
	for (int v : val) {
    if (v >= n) break;
		ll x = expmod(v, d, n);
		if (x == 1 || x == n - 1) continue;
		for (ll r=1; r<s; r++) if ((x = big(x)*x % n) == n - 1) goto nextPr;
		return false;
		nextPr:;
	}
	return true;
}

// Brute-force test for NTT parameters and convolutions  -- ONLY WORKS IF YOU
// HAVE __int128, otherwise the biggest mods overflow.
void test_ntt() {
  // Check that the parameters are good
  for (auto& pset : params) {
    ll mod = pset[0], g = pset[1], ord = pset[2];
    cerr << "Parameters:: (" << mod << ", " << g << ", " << ord << ")" << endl << '\t';
    assert(is_prime(mod));            // check that mod is prime
    assert(expmod(g,ord,mod) == 1);   // check that the order is correct
    cerr << "\tParameters are valid... testing convolution..." << endl << '\t';
    // If the parameters are good, do a convolution to check
    vector<ll> a = {1,1}, b = {1,1}, ans = {1,2,1};
    convolution<ll> conv{mod, g, ord};
    auto res = conv.multiply(a,b);
    for (int i = 0; i < 3; i++) assert(res[i] == ans[i]);
    cerr << "\tConvolution was good." << endl;
  }
}

int main() {
  //test_ntt();
	//solve_UVA12879();
  solve_CF632E();
  //solve_CF300D();
}
