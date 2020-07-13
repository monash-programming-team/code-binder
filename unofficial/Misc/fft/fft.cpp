// Fast convolution using Fast Fourier Transform
//
// Author: Daniel (based on e-maxx.ru)
// Date: 19-01-2016
// Reliability: 5
// Tested on: SPOJ-MUL, SPOJ-VFMUL, SPOJ-POLYMUL, CF286E, UVA12879
//
// Complexity: O(N logN)
#include<bits/stdc++.h>
using namespace std;

typedef long long int ll;
typedef vector<int> vi;

//listings:fft
// Fast convolution using Fast Fourier Transform. Complexity: O(n log(n))
typedef complex<double> comp;
const double PI=acos(-1.0);

void fft(vector<comp> &a, int invert=0) {    // Compute the FFT of the polynomial
  int n=a.size(), i, j, len; comp w, u, v;   // whose coefficients are given by
  for(i=1, j=0;i<n;i++) {                    // the elements of a.
    int bit = n/2; for(; j >= bit; bit /= 2) j-=bit;
    j += bit; if(i < j) swap(a[i], a[j]);
  }
  for(len=2;len<=n;len<<=1) {
    double ang=2*PI/len*(invert?-1:1); comp wlen = polar(1.0, ang);
    for(i=0; i<n; i+=len) for(j=0, w=1; j < len/2; j++)
      u=a[i+j], v=a[i+j+len/2]*w, a[i+j]=u+v, a[i+j+len/2]=u-v, w*=wlen;
  }
  if(invert) for(i=0;i<n;i++) a[i]/=n;
}

// Compute the convolution a * b
template<typename T> vector<T> multiply(const vector<T>& a, const vector<T>& b) {
  int i, n;  vector<comp> fa(a.begin(), a.end()), fb(b.begin(), b.end());
  for(n=1;n<2*(int)max(a.size(), b.size());n*=2);
  fa.resize(n), fb.resize(n), fft(fa), fft(fb);
  for(i=0;i<n;i++) fa[i]*=fb[i];
  fft(fa, 1); vector<T> res(n);   // Remove rounding below if T is non-integral
  for(i=0;i<n;i++) res[i]=(T)(fa[i].real()+0.5);
  return res;
}
//listings:/fft

// Verdict: AC
void solve_SPOJ_MUL() {
  int n; cin >> n;
  while (n--) {
    string a, b; cin >> a >> b;
    vi l1, l2;
    for (char c : a) l1.push_back(c-'0');
    for (char c : b) l2.push_back(c-'0');
    reverse(l1.begin(), l1.end());
    reverse(l2.begin(), l2.end());
    vi res = multiply(l1, l2);
    while (res.back()==0) res.pop_back();
    if (res.empty()) res.push_back(0);
    string answer;
    for (int i=0; i<(int)res.size(); i++) {
      answer.push_back('0'+(res[i]%10));
      if (res[i] >= 10) {
        if (i == (int)res.size() - 1) res.push_back(0);
        res[i+1] += res[i] / 10;
      }
    }
    reverse(answer.begin(), answer.end());
    cout << answer << '\n';
  }
}

// Verdict: AC
void solve_SPOJ_VFMUL() {
  int n; cin >> n;
  while (n--) {
    string a, b; cin >> a >> b;
    vector<ll> l1, l2;
    for (char c : a) l1.push_back(c-'0');
    for (char c : b) l2.push_back(c-'0');
    reverse(l1.begin(), l1.end());
    reverse(l2.begin(), l2.end());
    vector<ll> res = multiply(l1, l2);
    while (res.back()==0) res.pop_back();
    if (res.empty()) res.push_back(0);
    string answer;
    for (int i=0; i<(int)res.size(); i++) {
      answer.push_back('0'+(res[i]%10));
      if (res[i] >= 10) {
        if (i == (int)res.size() - 1) res.push_back(0);
        res[i+1] += res[i] / 10;
      }
    }
    reverse(answer.begin(), answer.end());
    cout << answer << '\n';
  }
}

// Verdict: 
void solve_SPOJ_POLYMUL() {
  int T; cin >> T;
  while (T--) {
    int n; cin >> n;
    vector<ll> a(n+1), b(n+1);
    for (auto& x : a) cin >> x;
    for (auto& x : b) cin >> x;
    while (a.back() == 0) a.pop_back();
    if (a.empty()) a.push_back(0);
    while (b.back() == 0) b.pop_back();
    if (b.empty()) b.push_back(0);
    reverse(a.begin(), a.end());
    reverse(b.begin(), b.end());
    auto res = multiply(a, b);
    while (res.back() == 0) res.pop_back();
    reverse(res.begin(), res.end());
    res.resize(2*n+1);
    for (int i=0; i<2*n+1; i++) cout << res[i] << " \n"[i==2*n];
  }
}

// Verdict: AC
void solve_CF286E() {
  int n, m;
  cin >> n >> m;
  vi w(m);
  for (int i=0; i<n; i++) {
    int a; cin >> a;
    w[a] = 1;
  }
  vi w2 = multiply(w,w);
  vi ans;
  for (int i=1; i<=m; i++) {
    if (w[i] == 1 && w2[i] == 0) ans.push_back(i);
    else if (w2[i] > 0 && w[i] == 0) { cout << "NO\n"; return; }
  }
  cout << "YES\n";
  cout << ans.size() << '\n';
  for (auto& x : ans) cout << x << ' '; cout << '\n';
}

// Verdict: AC
void solve_UVA12879() {
  int N;
  while (cin >> N) {
    vi dist;
    for (int i=0; i<N; i++) {
      int k; cin >> k;
      if ((int)dist.size() <= k) dist.resize(k+1);
      dist[k] = 1;
    }    
    vi dist2 = multiply(dist, dist);
    int M, ans = 0; cin >> M;
    for (int i=0; i<M; i++) {
      int d; cin >> d;
      if ((d < (int)dist.size() && dist[d])
          || (d < (int)dist2.size() && dist2[d])) ans++;
    }
    cout << ans << '\n';
  }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  //solve_SPOJ_MUL();
  solve_SPOJ_VFMUL();
  //solve_SPOJ_POLYMUL();
  //solve_CF286E();
  //solve_UVA12879();
}
