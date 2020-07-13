#include "test-includes.h"
#include "../official/number-theory/number-theory.cpp"

// Primality check via trial division in O(sqrt(n))
bool trial_division(ll n) {
  if (n <= 1) return false;  if (n <= 3) return true;
  if (n % 2 == 0 || n % 3 == 0) return false;
  for (ll i = 5; i * i <= n; i += 6)
    if (n % i == 0 || n % (i + 2) == 0) return false;
  return true;
}

TEST_CASE( "Test Discrete Logs" , "[discrete_log]" ) {
  int maxmod=250;
  for (int mod=1;mod<=maxmod;mod++)
    for (int a=0;a<mod;a++)
      for (int b=0;b<mod;b++) {
	bool f = false;
	int p=1%mod;
	for (int ans=0;ans<mod;ans++,p=p*a%mod)
	  if (p==b) {
	    f = true;
	    REQUIRE( discrete_log(a,b,mod) == ans );
	    break;
	  }
	if(!f)
	  REQUIRE( discrete_log(a,b,mod) == -1 );
      }
}

TEST_CASE( "Brute-force test the modular inverse function", "[modular_inverse]"){
  uniform_int_distribution<uint32_t> uint_dist;
  mt19937 rng;
  
  int M = 1e9 + 7;
  int num_tests = 10e6;
  rng.seed(time(NULL));
  for (int i = 1; i <= num_tests; i++) {
    int a = (uint_dist(rng) + M) % M;
    if(a == 0) continue;
    int z = inv(a, M);
    int res = (ll(a) * z) % M;
    REQUIRE( res == 1 );
  }
}

TEST_CASE( "Brute-force prime sieve", "[prime_sieves]" ){
  const int maxn = 1e7;
  sieve(maxn);
  fast_sieve(maxn);
  for (int p : pr) {
    REQUIRE( fac[p] == p );
    REQUIRE( isprime[p] );
    REQUIRE( is_prime(p) );
    REQUIRE( trial_division(p) );
  }
  for (int x = 1; x <= 1e5; x++) {
    if (isprime[x]) {
      REQUIRE( fac[x] == x);
      REQUIRE( is_prime(x));
      REQUIRE( trial_division(x));
      REQUIRE( binary_search(pr.begin(), pr.end(), x));
    } else {
      REQUIRE( fac[x] != x );
      REQUIRE( !is_prime(x) );
      REQUIRE( !trial_division(x) );
      REQUIRE( !binary_search(pr.begin(), pr.end(), x) );
    }
  }
}

TEST_CASE( "Brute-force test for factorisation" , "[factorisation]" ){
  // Pre-compute
  const int maxn = 20000000;
  fast_sieve(maxn);
  for (int x = 1; x <= 1e5; x++) {
    auto F1 = slow_factors(x);  auto F2 = fast_factors(x);
    for (int i = 0; i < (int)F1.size(); i++)
      REQUIRE( F1[i] == F2[i] );
    for (int f : F2){
      REQUIRE( is_prime(f) );
      REQUIRE( x % f == 0 );
    }
  }
}

TEST_CASE( "Brute-force test primitive roots", "[primitive-roots]" ){
  const int maxn = 10010;
  fast_sieve(maxn);
  for (ll n = 2; n < maxn; n++) {
    auto roots = all_primitive_root(n);
    auto tot = phi(n);  auto fact = fast_factors(tot);
    vector<ll> fact2(fact.begin(), fact.end());
    for (auto root : roots)
      REQUIRE( is_proot(n, root, tot, fact2) );
    if( !roots.empty() )
      REQUIRE( (int)roots.size() == phi(phi(n)) );
  }
}
