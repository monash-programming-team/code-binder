#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
typedef vector<int> vi;

ll mult_mod(ll x, ll y, ll m){ return (__int128)x * y % m; }
__int128 fast_exp_mod(__int128 b, __int128 n, __int128 m){
    if (n == 0) return 1 % m;
    if (n % 2 == 0) return fast_exp_mod((b*b)%m, n/2, m);
    return (fast_exp_mod(b, n-1, m) * b) % m;
}

/* use this if __int128 not available
ll q_mod(ll x, ll m) { return (x >= m) ? x-m : x; }
ll mult_mod(ll x, ll y, ll m){
    ll r = 0;
    while (y){
	if(y % 2) r = q_mod(r+x,m);
	y >>= 1; x = q_mod(x << 1, m);
    } return r;
}
ll fast_exp_mod(ll b, ll n, ll m){
    if (n == 0) return 1 % m;
    if (n % 2 == 0) return fast_exp_mod(mult_mod(b, b, m), n/2, m);
    return mult_mod(fast_exp_mod(b, n-1, m), b, m);
}
*/

vector<bool> isprime;
vi primes;

void sieve(int N){
    isprime.assign(N + 1, 1);
    isprime[0] = isprime[1] = 0;
    for (int i = 2; i <= N; i++){
	if (!isprime[i]) continue;
	primes.push_back(i);
	for (ll j = (ll)i * i; j <= N; j += i){
	    isprime[j] = 0;
	}
    }
}

bool trial_division(ll n){
    if (n <= 1) return false;
    if (n <= 3) return true;
    if (n % 2 == 0 || n % 3 == 0) return false;
    for (ll i = 5; i * i <= n; i += 6){
	if (n % i == 0) return false;
	if (n % (i + 2) == 0) return false;
    }
    return true;
}

bool is_prime(ll n){
    if (n <= 1) return false;
    vi val = {2, 7, 61}; // n <= 2^32
    // vi val = {2, 13, 23, 1662803}; // n <= 10^12
    // vi val = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}; // n <= 2^64

    // perform sieve up to cube root of n
    for (int x : primes) { 
	if ((ll)x * x * x > n) break;
	if (n % x == 0) return n == x; 
    }
    ll s = __builtin_ctzll(n-1), d = (n-1) >> s;
    for (int v : val){
	if (v >= n) break;
	ll x = fast_exp_mod(v, d, n);
	if (x == 1 || x == n-1) continue;
	for (ll r=1; r<s; r++) if ((x = mult_mod(x, x, n)) == n-1) goto nextPr;
	return false;
    nextPr:;
    }
    return true;
}

ll F(ll x, ll n, ll c) { x=mult_mod(x,x,n)-c; return (x < 0 ? x + n : x); }

int main(){
    sieve(1000*1000);
    cout << is_prime(30000000001LL) << endl;
    cout << trial_division(30000000001LL) << endl;
    for (ll i = 1; i <= 1000000000; i++){
    	if (is_prime(i * i)) {
    	    cout << i << endl;
    	    assert(false);
    	}
    	if (i % 1000000 == 0) cout << "=====" << i << "=====" << endl;
    }
}
