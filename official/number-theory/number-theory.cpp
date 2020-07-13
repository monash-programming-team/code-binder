// All of the number theory algorithms
//
// Author: Daniel Anderson, Peter Whalan, Xin Wei Chow, e-maxx.ru, Darcy's binder
// Date: 10-12-2016
// Reliability: See individual functions
// Tested on: See individual functions
//

typedef long long ll;
typedef vector<int> vi;
typedef pair<int,int> pii;

//listings:int128
typedef __int128 big;  // Use this if necessary. Mainly needed for huge prime testing.
//listings:/int128

// Tested on: UVA1230, UVA374
//listings:expmod
// Binary exponentiation - compute a^b mod m. Complexity O(log(n))
ll expmod(big a, big b, big m) {
  big res=1%m;
  a %= m;
  for(; b; b /= 2) {
    if (b&1)
      res=res*a%m;
    a=a*a%m;
  }
  return res;
}
//listings:/expmod

// Tested on: UVA10104, UVA10090
//listings:euclidean
// Extended Euclidean Algorithm. Finds x,y such that
// ax + by = gcd(a,b). Returns gcd(a,b). Compexity: O(log(min(a,b)))
ll gcd(ll a, ll b, ll& x, ll& y) {
  if (b == 0) {
    y = 0;
    x = (a < 0) ? -1 : 1;
    return (a < 0) ? -a : a;
  } else {
    ll g = gcd(b, a%b, y, x);
    y -= a/b*x;
    return g;
  }
}
//listings:/euclidean

// Tested on: Brute force, UVA11904
//listings:inverse
// Multiplicative inverse of a mod m, for a,m coprime. Complexity: O(log(a))
ll inv(ll a, ll m) {
  ll x, y;
  gcd(m,a,x,y);
  return ((y % m) + m) % m;
}
//listings:/inverse

// Tested on: ECNA07G, UVA756
//listings:cra
// Chinese Remainder Algorithm. Solves x = a[i] mod m[i] for x mod lcm(m)
// for m[i] pairwise coprime. In general x = x0 + t*lcm(m) for all t.
ll cra(vi& a, vi& m) {
  int n = (int)a.size();
  big u = a[0], v = m[0];
  ll p, q, r, t;
  for (int i = 1; i < n; ++i) {
    r = gcd(v, m[i], p, q);
    t = v;
    if ((a[i] - u) % r != 0) {
      return -1;  // no solution!
    }
    v = v/r * m[i];
    u = ((a[i] - u)/r * p * t + u) % v;
  }
  if (u < 0)
    u += v;
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
  for (ll i = 2; i*i <= n; ++i)
    if (n % i == 0) {
      while (n % i == 0)
        n /= i;
      res -= res / i;
    }
  if (n > 1)
    res -= res / n;
  return res;
}
//listings:/phi

// Tested on: Brute force, UVA10394, UVA543
//listings:sieve
// Sieve for primality testing up to 10^8. Complexity: O(n log(log(n)))
vector<bool> isprime;
void sieve(int n) {
  isprime.assign(n + 1, true);
  isprime[0] = isprime[1] = false;
  for (ll i = 2; i * i <= n; ++i)
    if (isprime[i])
      for (ll j = i*i; j <= n; j += i)
        isprime[j] = false;
}

// Sieve for factoring up to 10^7. Complexity: O(n)
// fac contains a prime factor, pr is a list of primes.
vi fac, pr;
void fast_sieve(int n) {
  fac.assign(n + 1, 0);
  for (ll i = 2; i <= n; ++i) {
    if (fac[i] == 0)
      fac[i] = i, pr.push_back(i);
    for (int p : pr)
      if (p > fac[i] || i * p > n)
        break;
      else
        fac[i * p] = p;
  }
}
//listings:/sieve


//listings:mr_seeds
// Deterministic Miller-Rabin primality test. Complexity: O(log(n))
//vi val = {2, 7, 61};                                  // n <= 2^32
//vi val = {2, 13, 23, 1662803};                        // n <= 10^12
vi val = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};  // n <= 2^64 (Needs __int128)
//listings:/mr_seeds


// Tested on: Brute force, UVA10394, UVA543
//listings:prime_test
bool is_prime(ll n) {
  if (n < 2)
    return false;
  ll s = __builtin_ctzll(n-1), d = (n-1) >> s;
  for (int v : val) {
    if (v >= n)
      break;
    ll x = expmod(v, d, n);
    if (x == 1 || x == n - 1)
      continue;
    for (ll r=1; r<s; r++)
      if ((x = ((big(x)*x) % n)) == n - 1)
        goto nextPr;
    return false;
nextPr:
    ;
  }
  return true;
}
//listings:/prime_test

// Tested on: Brute-force, UVA583, UVA10392
//listings:factorise
// Factorises n in O(log(n)) using precomputed fast_sieve(N >= n).
vi fast_factors(int n) {
  vi res;
  while (n > 1) {
    res.push_back(fac[n]);
    n /= fac[n];
  }
  return res;
}

// Factorises n in O(sqrt(n)) with no precomputation.
vector<ll> slow_factors(ll n) {
  vector<ll> res;
  for (ll i = 2; i*i <= n; ++i)
    for( ; n % i == 0 ; n /= i)
      res.push_back(i);
  if (n > 1)
    res.push_back(n);
  return res;
}

// Finds one (not necessarily prime) factor of n.
// Works best on semi-primes (n = pq for p, q distinct primes)
// Does not work well on perfect powers -- check those separately.
// Expected complexity: O(n^(1/4)) (only a heuristic)
ll F(ll x,ll n,ll c) {
  x=big(x)*x%n-c;
  return (x < 0 ? x + n : x);
}
ll pollardRho(ll n) {
  ll i,c,b,x,y,z,g;
  for(g=0, c=3; g%n == 0; c++)
    for (g=b=x=y=z=1; g == 1; b *= 2, g = __gcd(z,n), z = 1, y = x)
      for (i=0; i<b; i++) {
        x = F(x,n,c);
        z = (big)z * abs(x-y) % n;
      }
  return g;
}

// Factorise a huge number (n <= 10^18). Expected Complexity: O(n^(1/3))
vector<ll> factor_huge(ll n) {
  vector<ll> res;
  for (ll i = 2; i*i*i <= n; ++i)
    for( ; n % i == 0 ; n /= i)
      res.push_back(i);
  if (n == 1)
    return res;
  ll sqrt_n = sqrt(n)+0.5;
  if (sqrt_n*sqrt_n == n) {
    res.push_back(sqrt_n);
    res.push_back(sqrt_n);
    return res;
  }
  if (is_prime(n))
    return res.push_back(n), res;
  ll q = pollardRho(n);
  res.push_back(q);
  res.push_back(n/q);
  return res;
}
//listings:/factorise

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
  ll tot = phi(n);			        // if n is prime, can use tot = n - 1
  auto fact = slow_factors(tot);		// use fast_factors if you need
  for (ll res=2; res<n; ++res) {
    bool ok = __gcd(res, n) == 1;
    for (int i = 0; i < (int)fact.size() && ok; ++i)
      ok &= expmod(res, tot / fact[i], n) != 1;
    if (ok)
      return res;
  }
  return -1;
}

// Get a list of all primitive roots
vector<ll> all_primitive_root(ll n) {
  ll tot = phi(n);	          // if n is prime, can use tot = n - 1
  auto fact = slow_factors(tot);  // use fast_factors if you need
  vector<ll> ans;
  for (ll res=2; res<n; ++res) {
    bool ok = __gcd(res,n) == 1;
    for (int i = 0; i < (int)fact.size() && ok; ++i)
      ok &= expmod(res, tot / fact[i], n) != 1;
    if (ok)
      ans.push_back(res);
  }
  return ans;
}

// Test whether r is a primitive root modulo p
bool is_proot(ll n, ll r, ll tot, vector<ll>& fact) {
  bool ok = __gcd(n, r) == 1;
  for (int i = 0; i < (int)fact.size() && ok; ++i)
    ok &= expmod(r, tot / fact[i], n) != 1;
  return ok;
}
//listings:/primitive_root

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
      ans = it->second * sq - i;
      break;
    }
  }
  // Optional: if you only need one solution, return ans
  vector<ll> res;
  if (ans == -1)
    return res;
  ll delta = (n-1) / __gcd (k, n-1);
  for (ll cur = ans % delta; cur < n - 1; cur += delta)
    res.push_back(expmod(g, cur, n));
  return res;
}
//listings:/discrete_root

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
ll mult(ll a,ll b,ll mod) {
  return big(a)*b%mod;
}
ll discrete_log(ll a,ll b,ll mod) {
  static pair<ll,ll> seen[5000000]; // Must be at least ceil(sqrt(M))
  ll M=mod, s=0, as=1, bas; // step size, a^s, ba^s
  for (; s*s<M; s++)
    as=mult(as,a,mod), bas=mult(b,as,mod), seen[s]= {bas,s+1};
  sort(seen,seen+s);
  for (ll i=1,ap=1,ct=0,p; i<=s && ct<=s; i++) {
    ap=mult(ap,as,mod); //(ll)ap*as%mod;
    int j=lower_bound(seen,seen+s,pair<ll,ll> {ap+1,0})-seen;
    for (; --j>=0 && seen[j].first==ap && ct<=s; ct++)
      if (expmod(a,p=(ll)i*s-seen[j].second,mod)==b)
        return p;
  }
  return -1;
}
//listings:/discrete_log
