// Inversion index of a sequence
//
// Computes the number of inversions between two sequences a and b.
// a and b must be permutations of one another. The number of inversions
// is the number of adjacent transpositions required to sort b into a
// (the number of swaps that insertion sort would do to sort b into a).
//
// Author: Daniel Anderson
// Date: 18-01-2017
// Reliability: 5
// Tested on: SPOJ-INVCNT, SPOJ-YODANESS, UVA299, UVA11858, UVA10810
// 
// Complexity: O(n log(n))
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;

// Fenwick Tree with point query and ranged updates
struct FenwickTree {
	int N;  vi A;
	FenwickTree(int n): N(n+1), A(N) {}                    // Create tree with n elements
	void adjust(int b,int v) { for (;b;b-=b&-b) A[b]+=v; }             // Add v to A[0,b)
	void adjust(int a,int b,int v) { adjust(b,v), adjust(a,-v); }      // Add v to A[a,b)
	int pq(int i) { int r=0; for (i++;i<N;i+=i&-i) r+=A[i]; return r;	}       // Get A[i]
};

//listings:inversion_index
// Computes the inversion index between sequence a and b. The inversion index
// is the minimum number of required adjacent transpositions required to sort
// b into a. a and b must be a permutation of one and other. Complexity: O(n log(n))
template<typename T> ll inversion_index(const vector<T>& a, const vector<T>& b) {
  int n = (int)a.size(); ll ans = 0;  map<T,int> index;
  FenwickTree ft(n);    // Point query, ranged update Fenwick Tree
  for (int i=0; i<n; i++) ft.adjust(i,i+1,i), index[b[i]] = i;
  for (auto& c : a) ans += ft.pq(index[c]), ft.adjust(index[c],n,-1);
  return ans;
}
//listings:/inversion_index

// Verdict: AC
void solve_SPOJ_INVCNT() {
  int t; cin >> t;
  while (t--) {
    int n; cin >> n;
    vi a(n);  for (auto& x : a) cin >> x;
    vi i = a; sort(i.begin(), i.end());
    cout << inversion_index(i, a) << '\n';
  }
}

// Verdict: AC
void solve_SPOJ_YODANESS() {
  int t; cin >> t;
  while (t--) {
    int n; cin >> n;
    vector<string> a(n), b(n);
    for (auto& x : a) cin >> x;
    for (auto& x : b) cin >> x;
    cout << inversion_index(a, b) << '\n';
  }
}

// Verdict: AC
void solve_UVA299() {
  int N; cin >> N;
  while (N--) {
    int L; cin >> L;
    vi A(L); for (auto& x : A) cin >> x;
    vi I = A; sort(I.begin(), I.end());
    cout << "Optimal train swapping takes " << inversion_index(I, A) << " swaps.\n";
  }
}

// Verdict: AC
void solve_UVA11858() {
  int n;
  while (cin >> n) {
    vector<ll> a(n); for (auto& x : a) cin >> x;
    vector<ll> i = a; sort(i.begin(), i.end());
    cout << inversion_index(i, a) << '\n';
  }
}

// Verdict: AC
void solve_UVA10810() {
  int n;
  while (cin >> n && n) {
    vi a(n); for (auto& x : a) cin >> x;
    vi i = a; sort(i.begin(), i.end());
    cout << inversion_index(i, a) << '\n';
  }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  //solve_SPOJ_INVCNT();
  //solve_SPOJ_YODANESS();
  //solve_UVA299();
  //solve_UVA11858();
  solve_UVA10810();
}