// 2SAT
//
// Author      : Xin Wei Chow
// Date        : 8th Sept 2016
// Reliability :
// Tested on   :
#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef pair<int, int> pii;
typedef vector<vi> vvi;
#define X first
#define Y second

// SCC -----------------------------------------------------
class SCC {
private:
    int n, comp;
    vvi g, gt;
    vi seq, vis;
    void dfs(int u, const vvi &adj){
	for (int v : adj[u]){
            if (vis[v] == -1){
                vis[v] = comp;
                dfs(v, adj);
            }
        }
        seq.push_back(u);
    }
public:
    SCC() {}
    SCC(int _n){
        n = _n;
        g.assign(n, vi()); gt.assign(n, vi());
    }
    void add_edge(int u, int v){
        g[u].push_back(v); gt[v].push_back(u);
    }
    pair<int, vi> find_SCC(){
        vis.assign(n, -1); comp = 0;
        for (int i = 0; i < n; i++){
            if (vis[i] == -1){
                vis[i] = comp;
                dfs(i, g);
            }
        }
        vis.assign(n, -1); comp = 0;
        for (int i = n-1; i >= 0; i--){
            int u = seq[i];
            if (vis[u] == -1){
                vis[u] = comp;
                dfs(u, gt);
                comp++;
            }
        }
        return {comp, vis};
    }
};

// ---------------------------------------------------------
int VAR(int x) { return 2*x; }
int NOT(int x) { return x^1; }

class twoSAT {
private:
    int n;
    SCC scc;

public:
  twoSAT(int _n) {
        n = _n;
        scc = SCC(2 * n);
    }

    void add_clause(int u, int v){
        if (u == NOT(v)) return;
        scc.add_edge(NOT(u), v);
        scc.add_edge(NOT(v), u);
    }
    void add_true(int u){ add_clause(u, u); }
    void add_false(int u) { add_clause(NOT(u), NOT(u)); }
    void add_xor(int u, int v){
	add_clause(u, v);
	add_clause(NOT(u), NOT(v));
    }

    pair<bool, vector<bool>> solve(){
        vi comp = scc.find_SCC().Y;
        vector<bool> val(n);
        for (int i = 0; i < 2 * n; i += 2){
            if (comp[i] == comp[i + 1]) return {false, val};
            val[i/2] = (comp[i] > comp[i + 1]);
        }
        return {true, val};
    }
};

void add(int a, int b, int c, int d, twoSAT &sat){
    sat.add_clause(a, c);
    sat.add_clause(b, c);
    sat.add_clause(b, d);
    sat.add_clause(a, d);
}

int main(){
    int tc; cin >> tc;
    while (tc--){
        int S, A, m; cin >> S >> A >> m;
        twoSAT sat(S + A);
        for (int i = 0; i < m; i++){
            int s1, a1, s2, a2;
            cin >> s1 >> a1 >> s2 >> a2;
            s1--; a1--; s2--; a2--;
            int a = VAR(s1), d = VAR(s2), c = VAR(a1 + S), b = VAR(a2 + S);
            if (s1 < s2){
                if (a1 < a2){
                    add(a, b, c, d, sat);
                } else if (a1 > a2) {
                    add(NOT(a), b, c, NOT(d), sat);
                } else {
		    sat.add_true(b);
                    //sat.add_clause(b, b);
                }
            } else if (s1 > s2){
                if (a1 < a2){
                    add(a, NOT(b), NOT(c), d, sat);
                } else if (a1 > a2) {
                    add(NOT(a), NOT(b), NOT(c), NOT(d), sat);
                } else {
		    sat.add_false(b);
                    //sat.add_clause(NOT(b), NOT(b));
                }
            } else {
                if (a1 < a2){
		    sat.add_true(a);
                    // sat.add_clause(a, a);
                } else if (a1 > a2) {
		    sat.add_false(a);
                    // sat.add_clause(NOT(a), NOT(a));
                } else {

                }
            }
        }
        auto p = sat.solve().X;
        if (p) cout << "Yes" << endl;
        else cout << "No" << endl;
    }
}
