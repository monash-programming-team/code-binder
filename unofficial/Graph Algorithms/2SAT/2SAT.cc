// 2SAT
//
// Author      : Xin Wei Chow
// Date        : 8th Sept 2016
// Reliability : 5
// Tested on   : UVA 10319 (Manhattan)
//             : CF 468B (Two Sets)
//             : SPOJ TORNJEVI
//             : NWERC2011 D (Piece it Together)
//             : SPOJ SUPSUP
//
// Determines if 2SAT problem is satisfiable
// If satisfiable: returns a possible variable assignment
//
// Complexity: O(V + E)
//
// NOTE:
// use VAR(x) to refer to a variable
// use NOT(VAR(x)) to refer to negation of a variable

#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef pair<int, int> pii;
typedef vector<vi> vvi;
#define X first
#define Y second

// ====================== SCC ===========================
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
// ====================== SCC ===========================

//listings:twosat
// 2-SAT solver. Include SCC code from graph algorithms. VAR(x) is variable x,
// NOT(VAR(x)) is the negation of variable x. Complexity: O(n + m)
int VAR(int x) { return 2*x; }
int NOT(int x) { return x^1; }

struct TwoSAT {
	int n;  SCC scc;
	// Create a 2-SAT equation with n variables
	TwoSAT(int n) : n(n), scc(2 * n) { }
	void add_or(int u, int v) {
		if (u == NOT(v)) return;
		scc.add_edge(NOT(u), v);  scc.add_edge(NOT(v), u);
	}
	void add_true(int u){ add_or(u, u); }
	void add_false(int u) { add_or(NOT(u), NOT(u)); }
	void add_xor(int u, int v) { add_or(u, v); add_or(NOT(u), NOT(v)); }
	pair<bool, vector<bool>> solve() {
		vi comp = scc.find_SCC().Y;  vector<bool> val(n);
		for (int i = 0; i < 2 * n; i += 2){
			if (comp[i] == comp[i + 1]) return {false, val};
			val[i/2] = (comp[i] > comp[i + 1]);
		}
		return {true, val};
	}
};
//listings:/twosat

int main(){
    // satisfiable
    TwoSAT cc(2);
    cc.add_or(VAR(0), VAR(1));
    cc.add_or(NOT(VAR(0)), NOT(VAR(1)));
    cout << cc.solve().X << endl;

    // not satisfiable
    cc = TwoSAT(3);
    cc.add_or(VAR(0), VAR(1));
    cc.add_or(NOT(VAR(0)), VAR(1));
    cc.add_or(NOT(VAR(1)), VAR(2));
    cc.add_or(NOT(VAR(1)), NOT(VAR(2)));
    cout << cc.solve().X << endl;
}
