#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
typedef pair<int, int> pii;
const int MAXN = 100000;
const int NMOD = 5;
const int MOD[] = {1000000103, 1000000321, 1000000447, 1000000637, 1000000891};
const int BASE = 256;

struct Hash {
    int h[NMOD];
    Hash () {}
    Hash(int x) {
	for (int i = 0; i < NMOD; i++) h[i] = x;
    }
    bool operator<(const Hash &other) const {
	for (int i = 0; i < NMOD; i++)
	    if (h[i] != other.h[i]) return h[i] < other.h[i];
	return false;
    }
    Hash operator+(const Hash &other){
	Hash res;
	for (int i = 0; i < NMOD; i++)
	    res.h[i] = (h[i] + other.h[i]) % MOD[i];
	return res;
    }
    Hash operator-(const Hash &other){
	Hash res;
	for (int i = 0; i < NMOD; i++)
	    res.h[i] = (h[i] - other.h[i] + MOD[i]) % MOD[i];
	return res;
    }
    Hash operator*(const Hash &other){
	Hash res;
	for (int i = 0; i < NMOD; i++)
	    res.h[i] = (ll) h[i] * other.h[i] % MOD[i];
	return res;
    }
};

Hash pw[MAXN], pre[MAXN];
string str;
int n;

void precalc(){
    pw[0] = Hash(1);
    for (int i = 1; i < MAXN; i++)
        pw[i] = pw[i-1] * Hash(BASE);
}

void gen_prefix(const string &str){
    int n = str.size();
    for (int i = 0; i < n; i++){
        pre[i] = pw[i] * Hash(str[i]);
        if (i) pre[i] = pre[i] + pre[i-1];
    }
}

Hash get_hash(int L, int R){
    Hash res = pre[R];
    if (L) res = (pre[R] - pre[L - 1]);
    res = res * pw[n - L - 1];
    return res;
}

int main(){
    precalc();
    string str = "Hello World";
    gen_prefix(str);

    return 0;
}
