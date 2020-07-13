#include <bits/stdc++.h>
using namespace std;

typedef long long int ll;

const int MOD = 1000000007;

const int N = 2000010;
int fact[N], invfact[N], der[N];

inline int add(int a, int b) { return ((ll)a + b) % MOD; }
inline int mul(int a, int b) { return ((ll)a * b) % MOD; }
inline int power(int a, int b) { int res = 1; while (b > 0) { if (b & 1) { res = mul(res, a); } b >>= 1; a = mul(a, a); } return res; }
inline int inv(int a) { return power(a, MOD - 2); }

inline int binomial(int n, int k) { return mul(fact[n], mul(invfact[k], invfact[n - k])); }
inline int arrangements(int n, int k) { return mul(fact[n], invfact[n - k]); }
inline int catalan(int n) { return mul(fact[2*n], mul(invfact[n+1], invfact[n])); }

void precompute() {
  fact[0] = invfact[0] = der[0] = 1;
  for (int i = 1; i < N; i++) {
    fact[i] = mul(fact[i - 1], i);
    invfact[i] = mul(invfact[i - 1], inv(i));
	der[i] = mul(der[i - 1], i);
    if (i & 1) der[i] = add(der[i], MOD - 1);
    else der[i] = add(der[i], 1);
  }
}

int main() {
  precompute();
  

  return 0;
}
