#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
typedef vector<int> vi;
typedef pair<int, int> pii;
typedef ll big;

char letter(int i) {
	if (i == 27) return ' ';
	else return ('A' - 1 + i);
}


ll gcd(ll a, ll b, ll& x, ll& y) {
	if (b == 0) { y = 0; x = (a < 0) ? -1 : 1; return (a < 0) ? -a : a; }
	else { ll g = gcd(b, a%b, y, x); y -= a/b*x; return g; }
}

big cra(vi& a, vi& m) {
	int n = (int)a.size();
	big u = a[0], v = m[0]; ll p, q, r, t;
	for (int i = 1; i < n; ++i) {
		r = gcd(v, m[i], p, q); t = v;
		if ((a[i] - u) % r != 0) { return -1; }  // no solution!
		v = v/r * m[i];  u = ((a[i] - u)/r * p * t + u) % v;
	}
	if (u < 0) u += v;
	return u;
}

void solve(vi& keys, vi& encoding) {
	string result = "";
	for (int e : encoding) {
		vi r(4);
		for (int i=3; i>=0; i--) {
			r[i] = e % 100;
			e /= 100;
		}
		int decode = cra(r, keys);
		vi l(3);
		for (int i=2; i>=0; i--) {
			l[i] = decode % 100;
			decode /= 100;
		}
		for (int i=0; i<3; i++) result += letter(l[i]);
	}

	string trimmed;
	trimmed.reserve(result.size());
	bool nonspace = false;
	for (int i=result.size()-1; i >= 0; i--) {
		if (result[i] != ' ') nonspace = true;
		if (nonspace) trimmed = result[i] + trimmed;	
	}
	cout << trimmed << endl;
}

int main(){
	int t; cin >> t;
	while (t--) {
		int n; cin >> n;
		vi encoding(n), keys(4);
		for (int i=0; i<4; i++) cin >> keys[i];
		for (int i=0; i<n; i++) cin >> encoding[i];
		solve(keys, encoding);
	}
	return 0;
}