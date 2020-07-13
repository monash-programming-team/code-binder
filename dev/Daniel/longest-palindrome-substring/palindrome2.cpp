// Finds the longest palindromic substrings of a string.
//
// Author      : Daniel Anderson (Based on Ryan's LPS)
// Date        : 07-12-2016
// Reliability : 2
// Tested On   : SPOJ - LPS, Random
//
// Find the location and length of all longest palindromic substrings.

#include<bits/stdc++.h>
using namespace std;

#include "../code-template/debug.h"

#define X first
#define Y second

typedef pair<int,int> pii;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:lps
// Finds the longest palindromic substrings of a. Returns the longest length and a
// vector of positions at which longest palindromes occur. Complexity: O(n)
pair<int, vi> longest_palindrome(const vi& a) {
  int n = 2 * a.size() + 1, b = 0, m = 0, res = 0; vi pos;
  vector<pii> R(n);  vi p(n, -1);  // -1 should be something not in the input
  for (int i = 1; i < n; i += 2) p[i] = a[i/2];
  for (int i = 1; i < n; i++) {
    int w = i < b ? min(R[2 * m - i].Y, b-i) : 0;
    for (int l = i-w-1, u = i+w+1; l >= 0 && u < n && p[l--] == p[u++]; w++);
    R[i] = {(i - w)/2, w};
    if (i + w > b) b = i + w, m = i;
    if (w > res) res = R[i].Y;
  }
  for (auto& x : R) if (x.Y == res) pos.push_back(x.X);
  return {res, pos};
}
//listings:/lps

// Longest palindrome substring at each point in the string. Complexity: O(n log(n))
vi longest_palindromes(const vi& a) {
  int n = 2*a.size()+1, b=0, m=0, h=0; vi r(n), p(n, -1);  // -1 should be something
  for (int i = 1; i < n; i += 2) p[i] = a[i/2];            // not be in the input
  for (int i = 1; i < n; i++) {
    int w = i < b ? min(r[2*m-i], b-i) : 0;
    for (int l = i-w-1, u = i+w+1; l >= 0 && u < n && p[l--] == p[u++]; w++);
    r[i] = w;  if (i+w > b) b = i+w, m = i;  if (w > h) h = r[i];
  }  // If you only need the longest length, just return h here. Complexity: O(n)
  vi res(n/2), st(n/2,-1), cl(n/2,-1); set<int> c;
  for (int i=1;i<n-1;i++) st[(i-r[i])/2] = i, cl[(i+r[i]-1)/2] = i;
  for (int i=0;i<n/2;i++) c.insert(st[i]), res[i]=*c.rbegin() - 2*i, c.erase(cl[i]);
  return res;
}

//brute force palindrome checking for randomized testing
bool is_palindrome(string& in, int a, int b) {
  string x = in.substr(a, b-a+1);
  string y = x; reverse(y.begin(), y.end());
  return (x == y);
}

//brute force lps checking for randomized testing
int slow_lps(string& in) {
  int best = 0, n = in.length();
  for (int i=0; i < n; i++) {
    for (int j=i; j<n; j++) {
      if (is_palindrome(in, i, j)) best = max(best, j - i + 1);
    }
  }
  return best;
}

//Runs randomized tests
void test() {
  srand(time(NULL));
  int num_tests = 20000, len = 100, alphabet = 3;
  for (int t=1; t<=num_tests; t++) {
    // Random string test ---------------------------------------------------
    string ran;
    for (int i=0; i < len; i++) ran.push_back(((char)(rand()%alphabet))+'a');
    vi a(ran.begin(), ran.end());
    cout << "Test " << t << "/" << num_tests << '\r' << flush;
    // Test Manacher's algorithm --------------------------------------------
    auto res = longest_palindrome(a);
    int ans1 = res.X, ans2 = slow_lps(ran);
    if (ans1 != ans2) {
      cout << ran << ": ";
      cout << "FAST: " << ans1 << "\tSLOW: " << ans2 << "\t\t" << (ans1 == ans2 ? "GOOD" : "BAD") << endl;
      return;
    }
    for (auto& x : res.Y) {
      assert(is_palindrome(ran, x, x + ans1 - 1));
    }
    // Test Longest Palindromes --------------------------------------------
    auto improved_res = longest_palindromes(a);
    assert(*max_element(improved_res.begin(),improved_res.end())==ans1);
    for (int i=0;i<len;i++) {
      assert(is_palindrome(ran, i, i + improved_res[i] - 1));
      if (i + improved_res[i] + 1 < len) assert(!is_palindrome(ran,i,i+improved_res[i]));
    }
  }
  cout << "All tests passed.            " << endl;
}

// Verdict: AC
void solve_spoj_lps() {
  int N;  string S;  cin >> N >> S;
  vi a(S.begin(), S.end());
  auto ans = longest_palindrome(a);
  cout << ans.X << endl;
}

int main() {
  test();
  //solve_spoj_lps();
  //string s = "abcba";
  //auto res = improved_lps(vi{s.begin(),s.end()});
  //for (auto x : res) cout << x << ' '; cout << endl;
}
