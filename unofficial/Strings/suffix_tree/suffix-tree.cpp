// A linear time online suffix tree using Ukkonen's algorithm
//
// Author: Daniel (modified from http://codeforces.com/blog/entry/16780)
// Date: 11/10/2016
// Reliability: 5
// Tested on: divisions16-I, SPOJ-NHAY, SPOJ-DISUBSTR, SPOJ-SUBST1, SPOJ-LCS
//
// Complexity: O(n) to build, O(m) to search.
#include <bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef pair<int,int> pii;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:suffix_tree
// Linear time online suffix tree. Complexity: O(n) to build.
// Each node in the tree is indexed by an integer, starting from 0 as the root.
// to[u][c] is the node pointed to by node u along an edge beginning with char c.
// len[u] is the length of the parent edge of u (NOTE: len[u] may be greater than the
// length of the string if u is a leaf, ie. true length is min(len[u],n-fpos[u]))
// fpos[u] is an index of s containing the substring on the parent edge of u.
struct SuffixTree {
	const int INF = 1e9;	int node = 0, pos = 0, cap, sz = 1, n = 0;
	string s;  vi len, fpos, link;	vector<map<int, int>> to;
	int make_node(int _pos, int _len) {	fpos[sz] = _pos, len [sz] = _len; return sz++; }
	void add_letter(int c) {
		int last = 0;  s += c; n++; pos++;
		while(pos > 0) {
			while(pos > len[to[node][s[n - pos]]]) node=to[node][s[n-pos]], pos-=len[node];
			int edge = s[n - pos], &v = to[node][edge], t = s[fpos[v] + pos - 1];
			if (v == 0) v = make_node(n - pos, INF), link[last] = node, last = 0;
			else if (t == c) { link[last] = node;	return;	}
			else {
				int u = make_node(fpos[v], pos - 1);
				to[u][c] = make_node(n - 1, INF),	to[u][t] = v;
				fpos[v] += pos - 1,	len [v] -= pos - 1;
				v = u, link[last] = u, last = u;
			}
			if(node == 0) pos--;
			else node = link[node];
		}
	}
	SuffixTree(const string& S) : cap(2*S.size()), len(cap), fpos(cap), link(cap),
    to(cap) { len[0] = INF; s.reserve(S.size()); for (char c : S) add_letter(c); }
	SuffixTree(int N) : cap(2*N), len(cap), fpos(cap), // Create an empty suffix tree with
    link(cap), to(cap) { len[0] = INF; s.reserve(N); }      // capacity for N characters
	// Find the longest substring of the given pattern beginning at idx that matches a
	// substring in the tree. Returns {position, length} of the match. Complexity: O(m)
	pii longest_match(const string& pat, int idx) {
		int node = 0, jump = 0, ans = 0, m = (int)pat.size();
		if (to[node][pat[idx]] == 0) return {-1, 0};
		while (to[node][pat[idx]] > 0) {
			jump = 0;  node = to[node][pat[idx]];
			for (int i = fpos[node]; i < n && idx + jump < m
				&& jump < len[node] && pat[idx + jump] == s[i]; i++, jump++, ans++);
			if (jump < len[node]) break;
			idx += jump;
		}
		return {fpos[node] + jump - ans, ans};
  }
};
//listings:/suffix_tree

namespace problems {
  // Verdict: AC
  namespace DIVIS16_I {
    void solve() {
      string s, t;
      cin >> s >> t;
      SuffixTree tree(s);
      int n = (int)t.size();
      int pos = 0, ans = 0;
      while (pos < n) {
        auto match = tree.longest_match(t, pos);
        pos += max(match.Y, 1);
        ans++;
      }
      cout << ans << endl;
    }
  }
  
  // Verdict: AC
  namespace SPOJ_NHAY {
    string pat, str;
    void dfs(int u, int h, vi& depth, vi& isleaf, SuffixTree& st) {
      depth[u] = h;
      bool leaf = true;
      for (const auto& e : st.to[u]) if (e.Y) {
        int len = min(st.len[e.Y], (int)str.size() - st.fpos[e.Y]);
        dfs(e.Y, h + len, depth, isleaf, st);
        leaf = false;
      }
      if (leaf) isleaf[u] = 1;
    }
    void dfs2(int u, vi& leaves, const vi& isleaf, SuffixTree& st) {
      if (isleaf[u]) leaves.push_back(u);
      else for (const auto& e : st.to[u]) if (e.Y)
        dfs2(e.Y, leaves, isleaf, st);
    }
    void solve() {
      int n; bool first = true;
      while (cin >> n) {
        if (!first) cout << '\n'; first = false;
        cin >> pat >> str;
        // Construct the suffix tree and pre-process depths and leaves
        str += "$";
        SuffixTree st(str);
        vi depth(2*str.size()), isleaf(2*str.size()), leaves, loc;
        dfs(0, 0, depth, isleaf, st);
        // Traverse the suffix tree with the string pat
        if (st.to[0][pat[0]] == 0) continue;
        int node = st.to[0][pat[0]], pos = 1; bool good = true;
        for (int i=1; i<(int)pat.size(); i++) {
          if (pos == st.len[node]) {
            if (st.to[node][pat[i]] == 0) { good = false; break; }
            else node = st.to[node][pat[i]], pos = 1;
          } else {
            if (str[st.fpos[node]+pos] != pat[i]) { good = false; break; }
            else pos++;
          }
        }
        if (!good) continue;  // pat is not a substring of str
        // Find the leaves reachable from node
        dfs2(node, leaves, isleaf, st);
        // Compute positions relative to leaves
        for (int leaf : leaves) loc.push_back((int)str.size() - depth[leaf]);
        sort(loc.begin(), loc.end());
        for (int x : loc) cout << x << '\n';
      }
    }
  }
  
  // Verdict: AC
  // Note: Also solves SPOJ - SUBST1. Same input/output format, just larger bounds.
  namespace SPOJ_DISUBSTR {
    void solve() {
      int T; cin >> T;
      while (T--) {
        string S; cin >> S;
        int n = (int)S.size();
        SuffixTree st(S);
        int ans = 0;
        for (int u=1; u<2*n; u++) if (st.len[u] > 0)
          ans += min(st.len[u], n - st.fpos[u]);
        cout << ans << '\n';
      }
    }
  }
  
  // Verdict: TLE because judge time-limit is impossible. AC on random cases.
  namespace SPOJ_LCS {
    SuffixTree st(1000009); int ans = 0; int split;
    int depth[1000009]{}, left[1000009]{}, right[1000009]{}; 
    void dfs(int u, int h) {
      depth[u] = h;
      if (u && st.fpos[u] <= split && st.fpos[u] + st.len[u] > split) left[u] = 1;
      else if (u && st.fpos[u] > split && st.len[u] == st.INF) right[u] = 1;
      else for (const auto& e : st.to[u]) if (e.Y) {
        dfs(e.Y, h + st.len[e.Y]);
        if (left[e.Y]) left[u] = 1;
        if (right[e.Y]) right[u] = 1;
        if (u && left[u] && right[u]) ans = max(ans, depth[u]);
      }
    }
    void solve() {
      // Read strings and build composite suffix tree
      string S, T; cin >> S >> T; split = (int)S.size();
      for (char c :S) st.add_letter(c);  st.add_letter('$');
      for (char c : T) st.add_letter(c);  st.add_letter('#');
      // Compute depths and leaf reachability
      dfs(0, 0);
      // Find the deepest internal node that can reach '$' and '#'
      cout << ans << '\n';
    }
  }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  //problems::DIVIS16_I::solve();
  //problems::SPOJ_NHAY::solve();
  //problems::SPOJ_DISUBSTR::solve();
  problems::SPOJ_LCS::solve();
}
