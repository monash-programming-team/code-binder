// All of the number theory algorithms
//
// Author: Daniel Anderson, Peter Whalan, Xin Wei Chow, e-maxx.ru, Darcy's binder
// Date: 10-12-2016
// Reliability: See individual functions
// Tested on: See individual functions
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef pair<int,int> pii;
//typedef ll big;		// Use __int128 when available

//listings:int128
typedef __int128 big;  // Use this if necessary. Mainly needed for huge prime testing.
//listings:/int128

// Tested on: UVA1230, UVA374
//listings:expmod
// Binary exponentiation - compute a^b mod m. Complexity O(log(n))
ll expmod(big a, big b, big m) {
	big res=1%m;  a %= m;
	for(; b; b /= 2) { if (b&1) res=res*a%m;  a=a*a%m; }
	return res;
}
//listings:/expmod

// Tested on: UVA10104, UVA10090
//listings:euclidean
// Extended Euclidean Algorithm. Finds x,y such that
// ax + by = gcd(a,b). Returns gcd(a,b). Compexity: O(log(min(a,b)))
ll gcd(ll a, ll b, ll& x, ll& y) {
	if (b == 0) { y = 0; x = (a < 0) ? -1 : 1; return (a < 0) ? -a : a; }
	else { ll g = gcd(b, a%b, y, x); y -= a/b*x; return g; }
}
//listings:/euclidean

// Tested on: Brute force, UVA11904
//listings:inverse
// Multiplicative inverse of a mod m, for a,m coprime. Complexity: O(log(a))
ll inv(ll a, ll m) { ll x, y; gcd(m,a,x,y); return ((y % m) + m) % m; }
//listings:/inverse

// Tested on: ECNA07G, UVA756
//listings:cra
// Chinese Remainder Algorithm. Solves x = a[i] mod m[i] for x mod lcm(m)
// for m[i] pairwise coprime. In general x = x0 + t*lcm(m) for all t.
ll cra(vi& a, vi& m) {
	int n = (int)a.size();	big u = a[0], v = m[0]; ll p, q, r, t;
	for (int i = 1; i < n; ++i) {
		r = gcd(v, m[i], p, q); t = v;
		if ((a[i] - u) % r != 0) { return -1; }  // no solution!
		v = v/r * m[i];  u = ((a[i] - u)/r * p * t + u) % v;
	}
	if (u < 0) u += v;
	return u;
}
//listings:/cra

// Tested on: SPOJ-ETF, UVA10179, CF284A
//listings:phi
// Euler Phi Function - Count the integers coprime to n. Facts:
// (1) If p is prime, phi(p) = p - 1. (2) If p is prime, then
// phi(p^k) = p^k - p^(k-1). (3) If a and b are relatively
// prime, then phi(ab) = phi(a)phi(b). (4) If a and b are relatively
// prime, then a^phi(m) = 1 mod m.  Complexity: O(sqrt(n))
ll phi(ll n) {
	ll res = n;
	for (ll i = 2; i*i <= n; ++i)	if (n % i == 0) {
    while (n % i == 0) n /= i;
    res -= res / i;
	}
	if (n > 1) res -= res / n;
	return res;
}
//listings:/phi

// Tested on: Brute force, UVA10394, UVA543
//listings:sieve
// Sieve for primality testing up to 10^8. Complexity: O(n log(log(n)))
vector<bool> isprime;
void sieve(int n) {
	isprime.assign(n + 1, 1);	isprime[0] = isprime[1] = 0;
	for (ll i = 2; i * i <= n; ++i) if (isprime[i])
		for (ll j = i*i; j <= n; j += i) isprime[j] = 0;
}

// Sieve for factoring up to 10^7. Complexity: O(n)
// fac contains a prime factor, pr is a list of primes.
vi fac, pr;
void fast_sieve(int n) {
	fac.assign(n + 1, 0);
	for (ll i = 2; i <= n; ++i) {
		if (fac[i] == 0) fac[i] = i, pr.push_back(i);
		for (int p : pr) if (p > fac[i] || i * p > n) break; else fac[i * p] = p;
	}
}
//listings:/sieve

// Primality check via trial division in O(sqrt(n))
bool trial_division(ll n) {
	if (n <= 1) return false;  if (n <= 3) return true;
  if (n % 2 == 0 || n % 3 == 0) return false;
  for (ll i = 5; i * i <= n; i += 6)
  if (n % i == 0 || n % (i + 2) == 0) return false;
  return true;
}

/*
//listings:mr_seeds
// Deterministic Miller-Rabin primality test. Complexity: O(log(n))
vi val = {2, 7, 61};                                    // n <= 2^32
vi val = {2, 13, 23, 1662803};                          // n <= 10^12
vi val = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};  // n <= 2^64 (Needs __int128)
//listings:/mr_seeds
*/

// Tested on: Brute force, UVA10394, UVA543
// vi val = {2, 7, 61}; // n <= 2^32
// vi val = {2, 13, 23, 1662803}; // n <= 10^12
vi val = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}; // n <= 2^64 (Requires __int128)
//listings:prime_test
bool is_prime(ll n) {
	if (n < 2) return false;
	ll s = __builtin_ctzll(n-1), d = (n-1) >> s;
	for (int v : val) {
    if (v >= n) break;
		ll x = expmod(v, d, n);
		if (x == 1 || x == n - 1) continue;
		for (ll r=1; r<s; r++) if ((x = ((big(x)*x) % n)) == n - 1) goto nextPr;
		return false;
		nextPr:;
	}
	return true;
}
//listings:/prime_test

// Tested on: Brute-force, UVA583, UVA10392
//listings:factorise
// Prime factors in O(log(n)) using precomputed fast_sieve(N >= n).
vi fast_factors(int n) {
	vi res;
	while (n > 1) {
		int f = fac[n];
		while (n % f == 0) n /= f;  // remove while to include duplicates
		res.push_back(f);
	}
	return res;
}

// Prime factors in O(sqrt(n)) with no precomputation.
vector<ll> slow_factors(ll n) {
	vector<ll> res;
	for (ll i = 2; i*i <= n; ++i)	if (n % i == 0) {  // change if to while for duplicates
    res.push_back(i);
    while (n % i == 0) n /= i;  // remove while to include duplicates
  }
	if (n > 1) res.push_back(n);
	return res;
}

// Finds one (not necessarily prime) factor of n.
// Works best on semi-primes (n = pq for p, q distinct primes)
// Does not work well on perfect powers -- check those separately.
// Expected complexity: O(n^(1/4)) (only a heuristic)
ll F(ll x,ll n,ll c){ x=big(x)*x%n-c; return (x < 0 ? x + n : x); }
ll pollardRho(ll n) {
	ll i,c,b,x,y,z,g;
	for(g=0, c=3; g%n == 0; c++)
		for (g=b=x=y=z=1; g == 1; b *= 2, g = __gcd(z,n), z = 1, y = x)
			for (i=0; i<b; i++) { x = F(x,n,c); z = (big)z * abs(x-y) % n; }
	return g;
}

// Factorise a huge number (n <= 10^18). Expected Complexity: O(n^(1/3))
vector<ll> factor_huge(ll n) {
	vector<ll> res;
	for (ll i = 2; i*i*i <= n; ++i)	if (n % i == 0) { // change if to while for duplicates
    res.push_back(i);
    while (n % i == 0) n /= i;  // remove while to include duplicates
  }
  if (n == 1) return res; // Below, push_back(sqrt(n)) twice for duplicates
	if (ll(sqrt(n))*ll(sqrt(n)) == n) return res.push_back(sqrt(n)), res;
	if (is_prime(n)) return res.push_back(n), res; 
	ll q = pollardRho(n);  res.push_back(q);  res.push_back(n/q);
	return res;
}
//listings:/factorise

// Fast factors with duplicates
vi fast_factors_wdupe(int n) {
	vi res;
	while (n > 1) {
		int f = fac[n];
		n /= f;
		res.push_back(f);
	}
	return res;
}

// Slow factors with duplicates
vector<ll> slow_factors_wdupe(ll n) {
	vector<ll> res;
	for (ll i = 2; i*i <= n; ++i)
		while (n % i == 0) {
			res.push_back(i);
			n /= i;
		}
	if (n > 1) res.push_back(n);
	return res;
}

// Factoring huge numbers with duplicates
vector<ll> factor_huge_wdupe(ll n) {
	vector<ll> res;
	for (ll i = 2; i*i*i <= n; ++i)
		while (n % i == 0) {
			res.push_back(i);
			n /= i;
		}
	if (is_prime(n)) return res.push_back(n), res;
	if (ll(sqrt(n))*ll(sqrt(n)) == n)
    return res.push_back(sqrt(n)), res.push_back(sqrt(n)), res;
	ll q = pollardRho(n);  res.push_back(q);  res.push_back(n/q);
	return res;
}

// Tested on: SPOJ - PROOT, CF284A, Brute-force
//listings:primitive_root
// Find a primitive root modulo n. g is a primitive root modulo n if
// all coprimes to n are congruent to a power of g (mod n), ie. for any a
// such that gcd(a,n) = 1, there is k such that g^k = a (mod n) where k
// is the index or discrete logarithm of a to g (mod n). A primitive root
// exists only if n = 1,2,4 or n is a power of an odd prime or twice
// the power of an odd prime. The number of primitive roots is phi(phi(n)).
// Complexity: O(g log(phi(n)) log(n)). Returns -1 if no root exists.
ll primitive_root(ll n) {
	ll tot = phi(n);				          // if n is prime, can use tot = n - 1
  auto fact = slow_factors(tot);		// use fast_factors if you need
	for (ll res=2; res<n; ++res) {
		bool ok = __gcd(res, n) == 1;
		for (int i = 0; i < (int)fact.size() && ok; ++i)
			ok &= expmod(res, tot / fact[i], n) != 1;
		if (ok) return res;  // Can add to a vector and find all of them if needed
	}
	return -1;
}
//listings:/primitive_root

vector<ll> all_primitive_root(ll n) {
	ll tot = phi(n);				          // if n is prime, can use tot = n - 1
  auto fact = slow_factors(tot);		// use fast_factors if you need
  vector<ll> ans;
	for (ll res=2; res<n; ++res) {
		bool ok = __gcd(res,n) == 1;
		for (int i = 0; i < (int)fact.size() && ok; ++i)
			ok &= expmod(res, tot / fact[i], n) != 1;
		if (ok) ans.push_back(res);
	}
	return ans;
}

// Tested on: SPOJ - DISCRT1
//listings:discrete_root
// Discrete root solver - Given a prime n and integers a, k, we want
// to find all x satisfying x^k = a (mod n). Complexity: O(sqrt(n) log(n))
vector<ll> discrete_root(ll n, ll k, ll a) {
	ll g = primitive_root(n);		// n must be prime
	ll sq = (ll)sqrt(n) + 1;
	vector<pair<ll,ll>> dec(sq);
	for (ll i = 1; i <= sq; ++i)
		dec[i-1] = {expmod(g, i * sq * k % (n - 1), n), i};
	sort(dec.begin(), dec.end());
	ll ans = -1;
	for (ll i=0; i<sq; ++i) {
		ll my = expmod(g, i * k % (n - 1), n) * a % n;
		auto it = lower_bound(dec.begin(), dec.end(), make_pair(my, 0LL));
		if (it != dec.end() && it->first == my) {
			ans = it->second * sq - i;  break;
		}
	}
	// Optional: if you only need one solution, return ans
	vector<ll> res;  if (ans == -1) return res;
	ll delta = (n-1) / __gcd (k, n-1);
	for (ll cur = ans % delta; cur < n - 1; cur += delta)
		res.push_back(expmod(g, cur, n));
	return res;
}
//listings:/discrete_root

ll mult(ll a,ll b,ll mod) {
	return big(a)*b%mod;
}

// Written by Peter
//listings:discrete_log
// Discrete Logarithm. Complexity: O(sqrt(M)log(M))
// Solves a^x == b (mod mod) for integer x. The smallest non-negative x is chosen.
// Returns -1 if there is no such x. x is assumed to be strictly less than M;  
// To optimise set M to
//	 phi(mod/gcd(mod,lcm(a^lg(mod),b^lg(mod)))) + lg(mod)
//	 lg(mod) is the maximum multiplicity of a prime factor. This is length of
//	 path before entering cycle.
// Can also derive extra conditions to determine when solution exists before
// running algorithm.
ll discrete_log(ll a,ll b,ll mod) {
  static pair<ll,ll> seen[5000000]; // Must be at least ceil(sqrt(M))
  ll M=mod, s=0, as=1, bas; // step size, a^s, ba^s
  for (;s*s<M;s++) as=mult(as,a,mod), bas=mult(b,as,mod), seen[s]={bas,s+1};
  sort(seen,seen+s);
	for (ll i=1,ap=1,ct=0,p;i<=s && ct<=s;i++) {
		ap=mult(ap,as,mod); //(ll)ap*as%mod;
		int j=lower_bound(seen,seen+s,pair<ll,ll>{ap+1,0})-seen;
		for (;--j>=0 && seen[j].X==ap && ct<=s;ct++)
			if (expmod(a,p=(ll)i*s-seen[j].Y,mod)==b) return p;
	}
	return -1;
}
//listings:/discrete_log

// ------------------------------------------------------------------
//						TEST
// ------------------------------------------------------------------

// Verdict: AC
void solve_spoj_etf() {
	int T; cin >> T;
	while (T--) {
		int n; cin >> n;
		cout << phi(n) << endl;
	}
}

void test_discrete_log() {
		int maxmod=250;
		for (int mod=1;mod<=maxmod;mod++) for (int a=0;a<mod;a++) for (int b=0;b<mod;b++) {
			bool f=0;
			int p=1%mod;
			for (int ans=0;ans<mod;ans++,p=p*a%mod) if (p==b) {
				f=1;
				if (discrete_log(a,b,mod)!=ans) {
					cerr << a << ' ' << b << ' ' << mod << ' ' << discrete_log(a,b,mod) << ' ' << ans << endl;
					return;
				}
				break;
			}
			if (!f && discrete_log(a,b,mod)!=-1) {
					cerr << a << ' ' << b << ' ' << mod << ' ' << discrete_log(a,b,mod) << endl;
					//return;
			}
		}
		cerr << "success\n";
	}

// Verdict: AC
void solve_UVA10179() {
  ll n;
  while (cin >> n, n > 0) {
    cout << phi(n) << endl;
  }
}

// Verdict: AC
void solve_UVA1230() {
  int c; cin >> c;
  while (c--) {
    ll x, y, n; cin >> x >> y >> n;
    cout << expmod(x, y, n) << endl;
  }
}

// Verdict: AC
void solve_UVA374() {
  int B, P, M;
  while (cin >> B) {
    cin >> P >> M;
    cout << expmod(B, P, M) << endl;
  }
}

// Verdict: AC
void solve_UVA10104() {
  int A, B; ll X, Y, D;
  while (cin >> A) {
    cin >> B;
    D = gcd(A, B, X, Y);
    cout << X << ' ' << Y << ' ' << D << endl;
  }
}

// Verdict: AC
void solve_UVA10090() {
  ll n, c1, c2, n1, n2, x0, y0;
  while (cin >> n, n > 0) {
    cin >> c1 >> n1 >> c2 >> n2;
    ll G = gcd(n1, n2, x0, y0);
    if (n % G != 0) { cout << "failed" << endl; continue; }
    else {
      ll x = x0 * n / G, y = y0 * n / G;
      ll d1 = n2 / G, d2 = n1 / G;
      // Make x non-negative
      if (x < 0) {
        ll q = -(x-d1+1)/d1;
        x += q * d1, y -= q * d2;
        if (y < 0) { cout << "failed" << endl; continue; }
      }
      // Make y non-negative
      if (y < 0) {
        ll q = -(y-d2+1)/d2;
        y += q * d2, x -= q * d1;
        if (x < 0) { cout << "failed" << endl; continue; }
      }
      // Minimise c1 x + c2 y
      if (c1 * d1 > c2 * d2) {
        ll q = x / d1;
        x -= q * d1, y += q * d2;
      } else if (c1 * d1 < c2 * d2) {
        ll q = y / d2;
        y -= q * d2, x += q * d1;
      }
      // Output
      cout << x << ' ' << y << endl;
    }
  }
}

// Brute-force test the modular inverse function
uniform_int_distribution<uint32_t> uint_dist;
mt19937 rng;
void test_inv() {
  int M = 1e9 + 7;
  int num_tests = 10e6;
  rng.seed(time(NULL));
  for (int i = 1; i <= num_tests; i++) {
    int a = (uint_dist(rng) + M) % M;
    int z = inv(a, M);
    int res = (ll(a) * z) % M;
    if (i % 1000 == 0) cout << "Test " << i << " / " << num_tests << '\r' << flush;
    assert(res == 1);
  }
  cout << "All tests passed!                " << endl;
}

// Verdict: AC
const int MOD = 1e9 + 7;
const int MAXN = 1000100;
ll fact[MAXN];
void solve_UVA11904() {
  fact[0] = fact[1] = 1;
  for (int i=2; i<MAXN; i++) fact[i] = (fact[i-1] * i) % MOD;
  int T, c = 1; cin >> T;
  while (T--) {
    int n; cin >> n;  vi k(n);
    for (int i = 0; i < n; i++) cin >> k[i];
    ll sum = 0, ans = 1;
    for (auto x : k) {
      ans = (ans * fact[sum + x - 1]) % MOD;
      ans = (ans * inv((fact[x - 1] * fact[sum]) % MOD, MOD)) % MOD;
      sum += x;
    }
    cout << "Case " << c++ << ": " << ans << endl;
  }
}

// Verdict: AC
void solve_UVA756() {
  const char* msg = "Case %d: the next triple peak occurs in %d days.\n";
  const int g = __gcd(23, __gcd(28, 33)), lcm = 23 * 28 * 33 / g;
  int p, e, i, d, caseno = 1;
  while (cin >> p >> e >> i >> d, p != -1) {
    vi m = {23, 28, 33}, a = {p, e, i};
    int x = cra(a, m);
    while (x <= d) x += lcm;
    printf(msg, caseno++, x - d);
  }
}

// Brute-force primes
void test_primes() {
  const int maxn = 1e7;
  sieve(maxn); cout << "First sieve complete..." << endl;
  fast_sieve(maxn); cout << "Second sieve complete..." << endl;
  for (int p : pr) {
    assert(fac[p] == p);
    assert(isprime[p]);
    assert(is_prime(p));
    assert(trial_division(p));
  }
  for (int x = 1; x <= 1e5; x++) {
    if (isprime[x]) {
      assert(fac[x] == x);
      assert(is_prime(x));
      assert(trial_division(x));
      assert(binary_search(pr.begin(), pr.end(), x));
    } else {
      assert(fac[x] != x);
      assert(!is_prime(x));
      assert(!trial_division(x));
      assert(!binary_search(pr.begin(), pr.end(), x));
    }
  }
  cout << "All tests passed." << endl;
}

// Verdict: AC
void solve_UVA10394() {
  // Pre-compute
  const int maxn = 20000000;
  sieve(maxn);  fast_sieve(maxn);
  vector<pii> twins;  vector<bool> seen(maxn);
  for (int p : pr) seen[p] = true;
  for (int x = 1; x <= maxn; x++)
    if (isprime[x] && isprime[x+2]) {
      assert(is_prime(x) && is_prime(x+2));
      assert(seen[x] && seen[x+2]);
      assert(trial_division(x) && trial_division(x+2));
      twins.emplace_back(x, x+2);
    }
  // Answer input
  int S;
  while (cin >> S) printf("(%d, %d)\n", twins[S-1].X, twins[S-1].Y);
}

// Verdict: AC
void solve_UVA543() {
  // Pre-compute
  const int maxn = 20000000;
  sieve(maxn);  fast_sieve(maxn);
  // Answer queries
  int n;
  while (cin >> n, n != 0) {
    for (int p : pr) {
      if (isprime[n - p]) {
        assert(is_prime(n - p));
        assert(binary_search(pr.begin(), pr.end(), n - p));
        assert(trial_division(n - p));
        printf("%d = %d + %d\n", n, p, n - p);
        break;
      }
    } 
  }
}

// Brute-force test for factorisation
void test_factors() {
  // Pre-compute
  const int maxn = 20000000;
  fast_sieve(maxn);
  for (int x = 1; x <= 1e5; x++) {
    auto F1 = slow_factors(x);  auto F2 = fast_factors(x);
    for (int i = 0; i < (int)F1.size(); i++) assert(F1[i] == F2[i]);
    for (int f : F2) assert(is_prime(f) && x % f == 0);
  }
  cout << "Tests complete..." << endl;
}

// Verdict: AC
void solve_UVA583() {
  ll n;
  while (cin >> n, n != 0) {
    bool neg = false;
    if (n < 0) { n *= -1; neg = true; }
    auto F2 = slow_factors_wdupe(n);
    for (int f : F2) assert(is_prime(f) && n % f == 0);
    sort(F2.begin(), F2.end());
    if (neg) cout << '-';
    cout << n << " =";
    if (neg) cout << " -1 x";
    for (int i=0; i<(int)F2.size(); i++) {
      cout << ' ' << F2[i];
      if (i < (int)F2.size() - 1) cout << " x";
    }
    cout << endl;
  }
}

// Verdict: AC
// NOTE: Only works with __int128. WA if you don't have it. You'll
// also need to choose the biggest test-numbers for is_prime.
void solve_UVA10392() {
  ll n;
  while (cin >> n, n != -1) {
    auto res = factor_huge_wdupe(n);
    sort(res.begin(), res.end());
    for (auto x : res) {
      cout << "    " << x << '\n';
    }
    cout << '\n';
  }
}

// Test whether r is a primitive root modulo p
bool is_proot(ll n, ll r, ll tot, vector<ll>& fact) {
  bool ok = __gcd(n, r) == 1;
  for (int i = 0; i < (int)fact.size() && ok; ++i)
  	ok &= expmod(r, tot / fact[i], n) != 1;
	return ok;
}

// Brute-force test primitive roots
void test_proot() {
  const int maxn = 10010;
  fast_sieve(maxn);
  for (ll n = 2; n < maxn; n++) {
    if (n % 100 == 0) cout << "Testing " << n << " / 10,000" << '\r' << flush;
    auto roots = all_primitive_root(n);
    auto tot = phi(n);  auto fact = fast_factors(tot);
    vector<ll> fact2(fact.begin(), fact.end());
    for (auto root : roots) assert(is_proot(n, root, tot, fact2));
    assert(roots.empty() || (int)roots.size() == phi(phi(n)));
  }
  cout << "All tests passed          " << endl;
}

// Verdict: AC
void solve_SPOJ_PROOT() {
  ll p, n;
  while (cin >> p >> n, p > 0) {
    ll tot = p - 1;
    auto fact = slow_factors(tot);
    for (int i = 0; i < n; i++) {
      ll r; cin >> r;
      if (is_proot(p, r, tot, fact)) cout << "YES" << '\n';
      else cout << "NO" << '\n';
    }
  }
}

// Verdict: AC
void solve_CF284A() {
  ll p; cin >> p;
  cout << phi(phi(p)) << endl;
}

// Verdict: AC
void solve_spoj_discrt1() {
  ios_base::sync_with_stdio(0); cin.tie(0);
  ll n, k, a; cin >> n >> k >> a;
  auto res = discrete_root(n, k, a);
  sort(res.begin(), res.end());
  cout << res.size() << '\n';
  for (auto x : res) cout << x << '\n';
}

// Verdict: AC.
void solve_spoj_mod() {
  ios_base::sync_with_stdio(0); cin.tie(0);
  ll x, z, k;
  while (cin >> x >> z >> k, x > 0) {
    auto ans = discrete_log(x, k, z);
    if (ans == -1) cout << "No Solution" << '\n';
    else cout << ans << '\n';
  }
}


int main() {
  //solve_UVA10179();
  //solve_UVA1230();
  //solve_UVA374();
  //solve_UVA10104();
  //solve_UVA10090();
  //test_inv();
  //solve_UVA11904();
  //solve_UVA756();
  test_primes();
  //solve_UVA10394();
  //solve_UVA543();
  //test_factors();
  //solve_UVA583();
  //solve_UVA10392();
  //solve_SPOJ_PROOT();
  //solve_CF248A();
  //test_proot();
  //solve_spoj_discrt1();
  //test_discrete_log();
  //solve_spoj_mod();
}
