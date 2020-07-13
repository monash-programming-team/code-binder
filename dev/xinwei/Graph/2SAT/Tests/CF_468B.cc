// 2SAT
//
// Author      : Xin Wei Chow
// Date        : 8th Sept 2016
// Reliability : 2
// Tested on   : UVA 10319 (Manhattan)
//             : CF 468B (Two Sets)
//
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

int main(){
    ios::sync_with_stdio(false); cin.tie(0);
    int n, a, b; cin >> n >> a >> b;
    vi arr(n);
    for (int i = 0; i < n; i++) {
	cin >> arr[i];
    }
    map<int, int> mmap;
    for (int i = 0; i < n; i++){
	mmap[arr[i]] = i;
    }
    twoSAT sat(n);
    bool poss = true;
    for (int i = 0; i < n; i++){
	int x = -1, y = -1;
	if (mmap.count(a - arr[i])) x = mmap[a - arr[i]];
	if (mmap.count(b - arr[i])) y = mmap[b - arr[i]];
	if (x == -1 && y == -1){
	    poss = false;
	} else if (x == -1){
	    sat.add_true(VAR(i));
	    sat.add_true(VAR(y));
	} else if (y == -1){
	    sat.add_false(VAR(i));
	    sat.add_false(VAR(x));
	} else {
	    sat.add_clause(VAR(i), NOT(VAR(x)));
	    sat.add_clause(NOT(VAR(i)), VAR(y));
	}
    }
    auto p = sat.solve();
    if (p.X && poss) {
	cout << "YES" << endl;
	for (int x : p.Y)
	    cout << x << ' ';
	cout << endl;
    } else {
	cout << "NO" << endl;
    }

}
