#include <bits/stdc++.h>
using namespace std;

namespace {
#pragma GCC diagnostic ignored "-Wunused-function"

#define all(x) begin(x), end(x)
#define rall(x) rbegin(x), rend(x)
#define len(x) (int)((x).size())
#define X first
#define Y second

#define FOR(i, begin, end) for (__typeof(end) i = (begin) - ((begin) > (end)); i != (end) - ((begin) > (end)); i += 1 - 2 * ((begin) > (end)))

// Types
using ll = long long int;
using ld = long double;
using vi = vector<int>;
using vvi = vector<vi>;
using pii = pair<int,int>;
using pll = pair<ll,ll>;
using vll = vector<ll>;
template<typename T> using minpq = priority_queue<T, vector<T>, greater<T>>;

#ifndef ONLINE_JUDGE
#define DEBUG(args...) cerr << "[Line " << __LINE__ << "]: ";  { vector<string> _v = __split(#args, ','); __ferr(_v.begin(), args); }
#define DEBUG_2D(A) cerr << "[Line " << __LINE__ << "]: " << #A << " = \n";  for (const auto& R : (A)) { cerr << '\t' << R << '\n'; }
vector<string> __split(const string& s, char c) { vector<string> v; stringstream ss(s); string x; while (getline(ss, x, c)) v.emplace_back(x);  return v; }
void __err(vector<string>::iterator it) { cerr << endl; }
template<typename T, typename... Args> void __err(vector<string>::iterator it, T a, Args... args) { cerr << ", " << it -> substr((*it)[0] == ' ', it -> length()) << " = " << a; __err(++it, args...); }
template<typename T, typename... Args> void __ferr(vector<string>::iterator it, T a, Args... args) { cerr << it -> substr((*it)[0] == ' ', it -> length()) << " = " << a; __err(++it, args...); }
#else
#define DEBUG(...) 
#define DEBUG_2D(...)
#endif  // ONLINE_JUDGE

// Printing containers
template<typename U, typename V> ostream& operator<<(ostream &s, const pair<U, V> &x) { s << "(" << x.first << ", " << x.second << ")"; return s; }
template<typename U> ostream& operator<<(ostream &s, const vector<U> &x) { s << "["; bool was = false; for (auto it : x) { if (was) s << ", "; was = true; s << it; } s << "]";  return s; }
template<typename U> ostream& operator<<(ostream &s, const deque<U> &x) { s << "["; bool was = false; for (auto it : x) { if (was) s << ", "; was = true; s << it; } s << "]"; return s; }
template<typename U, typename V> ostream& operator<<(ostream &s, const map<U, V> &x) { s << "{"; bool was = false; for (auto it : x) { if (was) s << ", "; was = true; s << it; }  s << "}"; return s; }
template<typename U> ostream& operator<<(ostream &s, const set<U> &x) { s << "{"; bool was = false; for (auto it : x) { if (was) s << ", "; was = true; s << it; } s << "}"; return s; }
template<typename U> ostream& operator<<(ostream &s, const multiset<U> &x) { s << "{"; bool was = false; for (auto it : x) { if (was) s << ", "; was = true; s << it; } s << "}"; return s; }

// Useful functions
ll expmod(ll a,ll b,ll mod) {ll res=1;a%=mod;for(;b;b>>=1){if(b&1)res=res*a%mod;a=a*a%mod;}return res;}
template<typename T> T sqr(const T& x) { return x*x; }
ll flog(const ll x) { return 63 - __builtin_clzll(x); }
template<typename T> void sort(T& t) { sort(all(t)); }
template<typename T> void undupe(vector<T>& v) { sort(v); v.erase(unique(all(v)), v.end()); }
template<typename T> string binary(T x, int w=8*sizeof(T)) { string r; for (int i=0;i<w;i++) r=(((T(1)<<i)&x)?'1':'0')+r; return r; }
}

const ll MOD = 1000000009LL;
const ll INF = (ll)1e15;
const double EPS = 1e-8;

int main(int argc, char* argv[]) {
    if (argc >= 1) freopen(argv[1], "r", stdin);
    ios::sync_with_stdio(0); cin.tie(0);
    
    
}
