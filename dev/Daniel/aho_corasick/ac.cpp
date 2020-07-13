// Aho-Corasick Algorithm
//
// Author: Daniel Anderson
// Date: 26-01-2017
// Reliability: 1
// Tested on: Brute-force tests
//
// Reference: http://codeforces.com/blog/entry/14854
//
#include<bits/stdc++.h>
using namespace std;

#include "../code-template/debug.h"

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> pii;

//listings:ac
// Aho-Corasick dictionary matching automaton. Add dictionary words with add_key(key)
// then build_links(). To report every match in the text, report the contents of
// node[u].output for all suffix link ancestors u of v, for each node v on the path
// taken through the automaton by the text.
// Complexity: to build - O(M), to count matches - O(N), to report all matches - O(kN)
// k = no. of keys, M = total key length, N = text length.
struct AhoCorasick {
  struct State { map<char,int> edge; int link, cnt, tot; vi output; };
  int n, k; vector<State> node; vi len;
  int make_node() { node.emplace_back(); return n++; }
  void add_key(const string& y) {  // Add key y to the dictionary
    int v = 0;
    for (char c : y) {
      if(!node[v].edge[c]) node[v].edge[c] = make_node();
      v = node[v].edge[c];
    }
     node[v].cnt++, node[v].output.push_back(k++), len.push_back((int)y.size());
  }
  void build_links() {  // Call this once all keys have been inserted
    node[0].link = -1, node[0].tot = 0;  queue<int> q; q.push(0);
    while (!q.empty()) {
      int v = q.front(); q.pop(); node[v].tot = node[v].cnt;
      if (node[v].link != -1) node[v].tot += node[node[v].link].tot;
      for (auto it: node[v].edge) {
        int c = it.first, u = it.second, j = node[v].link;
        while (j != -1 && !node[j].edge[c]) j = node[j].link;
        if (j != -1) node[u].link=node[j].edge[c];
        q.push(u);
      }
    }
  }  // Create an empty Aho-Corasick automaton
  AhoCorasick() : n(1), k(0), node(1) { }  
  ll count_matches(const string& x) {  // Count the number of substrings of the given
    ll ans = 0;  int v = 0;                // text that match a key: Complexity: O(N)
    for (int i=0; i<(int)x.size(); i++) {
      while (v && !node[v].edge[x[i]]) v = node[v].link;
      v = node[v].edge[x[i]]; ans += node[v].tot;
    }
    return ans;
  }
};
//listings:/ac

/*
vector<pii> find_matches(const string& x) {  // Find all matches in the text x
  vector<pii> ans;  int v = 0;               // in the form {position, key_id}
  for (int i=0; i<(int)x.size(); i++) {      //             Complexity: O(N^2)
    while (v && !node[v].edge[x[i]]) v = node[v].link;
    v = node[v].edge[x[i]];
    for (int u = v; u != -1; u = node[u].link)
      for (auto& p : node[u].output) ans.emplace_back(i-len[p]+1,p);
  }
  return ans;
}
*/

// Brute-force counting the number of matches in a set of random keys
namespace brute_force_test {
  const int T = 1000; // number of tests
  const int k = 1000; // number of patterns
  const int N = 1000; // text length
  const int M1 = 1, M2 = 20; // pattern length [M1,M2]
  const int alpha = 5; // alphabet size
  void test_count() {
    srand(time(NULL));
    cout << "Testing count_matches()" << endl;
    for (int t=1; t<=T; t++) {
      cout << "Running test " << t << "/" << T << "            \r" << flush;
      AhoCorasick ac; vector<string> keys;
      // Build random string ----------------------------------------------------------
      string txt; for (int i=0; i < N; i++) txt.push_back(((char)(rand()%alpha))+'a');
      // Build random patterns --------------------------------------------------------
      ll num_matches = 0;
      for (int i=0; i<k; i++) {
        int len = M1 + (rand() % (M2-M1+1));
        string key; for (int i=0; i < len; i++) key.push_back(((char)(rand()%alpha))+'a');
        ac.add_key(key); keys.push_back(key);
        for (int i=0; i<N-len+1; i++) if (txt.substr(i,len) == key) num_matches++;
      }
      // Execute machine ---------------------------------------------------------------
      ac.build_links();
      ll ans = ac.count_matches(txt);
      assert(ans == num_matches);
    }
  }
}

int main() {
  brute_force_test::test_count();
}
