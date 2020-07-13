// 2SAT
//
// Author      : Xin Wei Chow
// Date        : 8th Sept 2016
// Reliability : 3
// Tested on   : UVA 10319 (Manhattan)
//             : CF 468B (Two Sets)
//             : SPOJ TORNJEVI
//             : NWERC2011 D (Piece it Together)
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

    void add_or(int u, int v){
        if (u == NOT(v)) return;
        scc.add_edge(NOT(u), v);
        scc.add_edge(NOT(v), u);
    }
    void add_true(int u){ add_or(u, u); }
    void add_false(int u) { add_or(NOT(u), NOT(u)); }
    void add_xor(int u, int v){
	add_or(u, v);
	add_or(NOT(u), NOT(v));
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

string yes = "March onward";
string no = "Coordination issue";
void add(vi &v, twoSAT &sat){
    int n = (int)v.size();
    for (int i = 0; i < n; i++){
        for (int j = i + 1; j < n; j++){
            sat.add_or(v[i], v[j]);
        }
    }
}
void solve(int n){
    twoSAT sat(n * 100);
    int nodes = 0;
    map<int, int> mmap;
    for (int i = 0; i < n; i++){
        int m; cin >> m;
        vi v;
        for (int j = 0; j < m; j++){
            string str; cin >> str;
            char c = str.back();
            str.pop_back();
            int num = stoi(str);
            if (!mmap.count(num)) mmap[num] = nodes++;
            int x = VAR(mmap[num]);
            if (c == 'S') x = NOT(x);
            v.push_back(x);
        }
        add(v, sat);
    }
    if (sat.solve().X) cout << yes << endl;
    else cout << no << endl;
}

int main(){
    ios::sync_with_stdio(false); cin.tie(0);

    int n;
    while (cin >> n && n){
        solve(n);
    }
}
