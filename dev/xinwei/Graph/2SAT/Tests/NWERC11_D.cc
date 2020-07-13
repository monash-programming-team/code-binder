#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef pair<int, int> pii;
typedef vector<vi> vvi;
char grid[512][512];
int node[512][512];
int dx[] = {-1, 0, +1, 0};
int dy[] = {0, -1, 0, +1};
int rev[] = {2, 3, 0, 1};

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

void add_nodes(const vi &nd, twoSAT &sat){
    int sz = (int)nd.size();
    if (sz == 1) {
        int A = nd[0];
        sat.add_true(VAR(A));
    } else if (sz == 2){
        int A = nd[0], B = nd[1];
        sat.add_clause(VAR(A), VAR(B));
    }
    for (int i = 0; i < sz; i++){
        for (int j = i + 1; j < sz; j++){
            sat.add_clause(NOT(VAR(nd[i])), NOT(VAR(nd[j])));
        }
    }
}

void solve(){
    int n, m; cin >> n >> m;
    int NNODES = n * m * 4;
    twoSAT sat(NNODES);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            cin >> grid[i][j];
    int b = 0, w = 0;
    int nodes = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++){
            b += (grid[i][j] == 'B');
            w += (grid[i][j] == 'W');
            if (grid[i][j] == 'W'){
                node[i][j] = nodes;
                nodes += 4;
            }
        }
    if (w != 2*b) {
        cout << "NO" << endl;
        return;
    }
    // process black squares
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            if (grid[i][j] == 'B'){
                // vert
                vector<int> nd;
                for (int k = 0; k < 4; k += 2){
                    int i2 = i + dx[k], j2 = j + dy[k];
                    if (i2 < 0 || i2 >= n || j2 < 0 || j2 >= m) continue;
                    if (grid[i2][j2] != 'W') continue;
                    nd.push_back(node[i2][j2] + rev[k]);
                }
                if ((int)nd.size() == 0){
                    cout << "NO" << endl;
                    return;
                }
                add_nodes(nd, sat);
                // hort
                nd.clear();
                for (int k = 1; k < 4; k += 2){
                    int i2 = i + dx[k], j2 = j + dy[k];
                    if (i2 < 0 || i2 >= n || j2 < 0 || j2 >= m) continue;
                    if (grid[i2][j2] != 'W') continue;
                    nd.push_back(node[i2][j2] + rev[k]);
                }
                if ((int)nd.size() == 0){
                    cout << "NO" << endl;
                    return;
                }
                add_nodes(nd, sat);
            }

    // process white squares
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            if (grid[i][j] == 'W'){
                vi nd;
                for (int k = 0; k < 4; k++){
                    int i2 = i + dx[k], j2 = j + dy[k];
                    if (i2 < 0 || i2 >= n || j2 < 0 || j2 >= m) continue;
                    if (grid[i2][j2] != 'B') continue;
                    nd.push_back(node[i][j] + k);
                }
                if ((int)nd.size() == 0){
                    cout << "NO" << endl;
                    return;
                }
                add_nodes(nd, sat);
            }
    vector<bool> val;
    if (sat.solve().X) cout << "YES" << endl;
    else cout << "NO" << endl;
}

int main(){
    ios::sync_with_stdio(false);
    cin.tie(0);
    int tc; cin >> tc;
    while (tc--)
        solve();
}
