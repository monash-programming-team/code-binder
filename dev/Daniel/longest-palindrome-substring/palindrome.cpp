#include<bits/stdc++.h>
using namespace std;

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

//Finds the length of the longest palindromic substring of a given string
//If you want to reconstruct the original string, record the index i of best
//and take the string within radius DP[i] around i (and remove padding!)
int lps(string& in) {
  int n = in.length(); string s; s.resize(2*n+4);
  s[0] = '#'; s[2*n+3] = '$'; //Make these two different chars not in the input string
  for (int i=1; i<2*n+3; i++) s[i] = (i % 2 == 1 ? '_' : in[(i-1)/2]); // Make _ something not in the input
  int DP[2*n+3]; fill(DP, DP + 2*n+3, 0);
  int center = 0, end = 0, best = 0;
  for (int i=1; i<2*n+3; i++) {
    if (i > end) center = i;
    else DP[i] = min(DP[2*center - i], end - i);
    int r = max(i+1, (i <= end ? min(end + 1, i + DP[2*center - i] + 1) : 0));
    int l = 2*i - r;
    while (s[r++] == s[l--]) DP[i]++;
    best = max(best, DP[i]);
    if (DP[i] + i > end) { end = DP[i] + i; center = i; }
  }
  return best;
}

//Runs randomized tests
int main() {
  srand(time(NULL));
  for (int t=0; t<100; t++) {
    string ran;
    for (int i=0; i<50; i++) ran.push_back(((char)(rand()%3))+'a');
    cout << ran << ": ";
    cout << "FAST: " << lps(ran) << "\tSLOW: " << slow_lps(ran) << "\t\t" << (lps(ran) == slow_lps(ran) ? "GOOD" : "BAD") << endl;
  }
  return 0;
}
