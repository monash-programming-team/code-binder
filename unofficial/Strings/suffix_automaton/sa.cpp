// Linear time online Suffix Automaton
//
// Author: Daniel (Based on code from National Taiwan University's code binder)
// Date: 23-01-2017
// Reliability: 5
// Tested on: SPOJ-DISUBSTR, SPOJ-SUBST1, SPOJ-SUB_PROB, SPOJ-NHAY, CC-SUBQUERY
//            SPOJ-LCS
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:sa
// Linear time online suffix automaton. node stores the states of the automaton.
// node[1] is the root. tail is the value of run(S).  Complexity: O(n) to build.
// State: par -- parent suffix link (edges of the suffix tree of the reverse of S)
//        pos -- length of prefix of S such that run(S[0..pos)) = node
//        edge[x] -- index of node following edge with character x (0 if none)
// Terminal states are all suffix link ancestors of tail (including tail).
// Useful facts: Each node corresponds to an equivalence class of strings w such
// that run(w) = node. Every string in this equivalence class is a suffix of the
// longest string W in the equivalence class. The suffix link leads to the state 
// corresponding to the equivalence class of the longest suffix of W that is not
// in the same equivalence class.
struct SuffixAutomaton{
	struct State{
    int par, pos; map<char,int> edge;
		State (int v) : par(0), pos(v) { }
	};
	vector<State> node;	int root, tail;
  SuffixAutomaton(const string& S) : root(1), tail(1) { // Create an automaton from S
    node.assign(2, State(0));  for (char c : S) extend(c);
  }
	void extend(char w) {  // Add a character to the string and extend the automaton
		int p = tail, np = node.size();	node.emplace_back(node[p].pos+1);
		for (; p && node[p].edge[w]==0; p=node[p].par)	node[p].edge[w] = np;
		if (p == 0)	node[np].par = root;
		else {
			if (node[node[p].edge[w]].pos == node[p].pos+1) node[np].par = node[p].edge[w];
			else {
				int q = node[p].edge[w], r = node.size();	node.push_back(node[q]);
				node[r].pos = node[p].pos+1, node[q].par = node[np].par = r;
				for (; p && node[p].edge[w] == q; p=node[p].par) node[p].edge[w] = r;
			}
		}
		tail = np;
	}                                  
  int run(const string& pat) {  // Return the node reached by running the machine
    int n = root;               // on the input pat, or 0 if pat is not a substring
    for (char c : pat) if ((n = node[n].edge[c]) == 0) return 0;
    return n;
  }
};
//listings:/sa

// ----------------------------------------------------------------------------
//                                   TEST PROBLEMS
// ----------------------------------------------------------------------------
namespace problems {
  // Verdict: AC
  // Note: Also solves SPOJ - SUBST1. Same input/output format, just larger bounds.
  namespace SPOJ_DISUBSTR {
    vector<ll> DP;
    ll f(int state, SuffixAutomaton& sa) {
      if (DP[state] != -1) return DP[state];
      ll& ans = DP[state] = 1;
      for (auto& e : sa.node[state].edge) if (e.Y) ans += f(e.Y, sa);
      return ans;
    }
    void solve() {
      int T; cin >> T;
      while (T--) {
        string s; cin >> s;
        SuffixAutomaton sa(s);
        DP.assign(sa.node.size(), -1);
        cout << f(1, sa) - 1 << '\n';
      }
    }
  }
  
  // Verdict: AC
  namespace CC_SUBQUERY {
    vector<ll> NP, SP, LP; vi terminal; vvi radj;
    ll np(int state, SuffixAutomaton& sa) {  // Number of paths that begin in state and go to a terminal
      if (NP[state] != -1) return NP[state];
      ll& ans = NP[state] = terminal[state];
      for (const auto& e : sa.node[state].edge) if (e.Y) ans += np(e.Y, sa);
      return ans;
    }
    ll sp(int state) {  // Shortest path to state
      if (state == 1) return SP[state] = 0;
      if (SP[state] != -1) return SP[state];
      ll& ans = SP[state] = LLONG_MAX;
      for (int v : radj[state]) ans = min(ans, sp(v) + 1);
      return ans;
    }
    ll lp(int state) {  // Longest path to state
      if (state == 1) return LP[state] = 0;
      if (LP[state] != -1) return LP[state];
      ll& ans = LP[state] = LLONG_MIN;
      for (int v : radj[state]) ans = max(ans, lp(v) + 1);
      return ans;
    }
    void solve() {
      // Build suffix automaton
      string S; cin >> S;
      SuffixAutomaton sa(S);
      int n = sa.node.size();
      // Determine terminal states
      terminal.resize(n);
      for (int p = sa.tail; p > 0; p = sa.node[p].par) terminal[p] = 1;
      // Build the reverse DAG of the suffix automaton
      radj.resize(n);
      for (int u=1; u<n; u++) {
        for (const auto& e : sa.node[u].edge) if (e.Y) {
          radj[e.Y].push_back(u);
        }
      }
      // Determine the intervals of occurrences for each node
      NP.assign(n,-1), SP.assign(n,-1), LP.assign(n,-1);
      vector<vector<ll>> start(S.size()+1), close(S.size()+1);
      for (int u=2; u<n; u++) {
        start[np(u,sa)].push_back(sp(u));
        close[np(u,sa)].push_back(lp(u));
      }
      for (int i=0; i<=(int)S.size(); i++) {
        sort(start[i].begin(), start[i].end());
        sort(close[i].begin(), close[i].end());
      }
      // Answer queries
      int N; cin >> N;
      while (N--) {
        int L, P; cin >> L >> P;
        if (P > (int)S.size()) cout << 0 << '\n';
        else {
          int ans = distance(lower_bound(close[P].begin(), close[P].end(), L), close[P].end())
                - distance(upper_bound(start[P].begin(), start[P].end(), L), start[P].end());
          cout << ans << '\n';
        }
      }
    }
  }
  
  // Verdict: AC
  namespace SPOJ_NHAY {
    void dfs(int u, vi& vis, const vvi& st) {
      vis[u] = 1; for (auto& v : st[u]) dfs(v,vis,st);
    }
    void solve() {
      int n; bool first = true;
      while (cin >> n) {
        if (!first) cout << '\n'; first = false;
        string pat, str;
        cin >> pat >> str;
        str = "$" + str;
        SuffixAutomaton sa(str);
        // Build the reverse suffix tree from the suffix automaton
        vvi suffix_tree(sa.node.size());
        for (int u=1; u<(int)sa.node.size(); u++) 
          if (sa.node[u].par) suffix_tree[sa.node[u].par].push_back(u);
        // Find the state corresponding to running pat
        int node = sa.run(pat);
        if (node == 0) continue;  // pat not found
        // DFS the suffix tree to find locations of pat
        vi vis(sa.node.size()), loc;
        dfs(node, vis, suffix_tree);
        for (int u=0; u<(int)sa.node.size(); u++)
          if (suffix_tree[u].empty() && vis[u]) 
            loc.push_back(sa.node[u].pos - pat.size() - 1);
        sort(loc.begin(), loc.end());
        for (int x : loc) cout << x << '\n';
      }
    }
  }
  
  // Verdict: AC
  namespace SPOJ_SUB_PROB {
    void solve() {
      string M; cin >> M;
      SuffixAutomaton sa(M);
      int N; cin >> N;
      while (N--) {
        string S; cin >> S;
        cout << (sa.run(S) ? 'Y' : 'N') << '\n';
      }
    }
  }
  
  // Verdict: AC
  namespace SPOJ_LCS {
    void solve() {
      string S, T;
      cin >> S >> T;
      SuffixAutomaton sa(S);
      int best = 0, len = 0, node = 1;
      for (char c : T) {
        auto it = sa.node[node].edge.find(c);
        if (it == sa.node[node].edge.end()) {
          while (node > 1 && it == sa.node[node].edge.end())
            node = sa.node[node].par, it = sa.node[node].edge.find(c);
          len = sa.node[node].pos;
        }
        if (it != sa.node[node].edge.end()) len++, node = it->Y;
        best = max(best, len);
      }
      cout << best << '\n';
    }
  }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  //problems::SPOJ_DISUBSTR::solve();
  //problems::SPOJ_NHAY::solve();
  //problems::SPOJ_SUB_PROB::solve();
  //problems::CC_SUBQUERY::solve();
  problems::SPOJ_LCS::solve();
}

