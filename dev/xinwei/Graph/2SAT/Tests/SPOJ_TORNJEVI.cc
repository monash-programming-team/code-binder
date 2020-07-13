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

char grid[128][128];
int tower[128][128];
int dx[] = {-1, 0, 0, 1};
int dy[] = {0, -1, 1, 0};
int dk[] = {3, 2, 1, 0};
int chosen[128][128][4];
int main(){
    ios::sync_with_stdio(false); cin.tie(0);
    int r, s; cin >> r >> s;
    for (int i = 0; i < r; i++) 
	for (int j = 0; j < s; j++)
	    cin >> grid[i][j];
    int n = 0;
    for (int i = 0; i < r; i++){
	for (int j = 0; j < s; j++){
	    if (grid[i][j] == 'T'){
		tower[i][j] = 4 * n;
		n += 1;
	    }	
	}
    }
    twoSAT sat(4 * n);
    for (int i = 0; i < r; i++){
	for (int j = 0; j < s; j++){
	    if (grid[i][j] == 'T'){
		sat.add_xor(VAR(tower[i][j]), VAR(tower[i][j] + 3));
		sat.add_xor(VAR(tower[i][j] + 1), VAR(tower[i][j] + 2));
	    }	
	}
    }
    memset(chosen, -1, sizeof chosen);
    for (int i = 0; i < r; i++)
	for (int j = 0; j < s; j++){
	    if (grid[i][j] == 'T') {
	    // check if collide with other tower
		for (int k = 0; k < 4; k++){
		    int i2 = i, j2 = j;
		    while (1){
			i2 += dx[k]; j2 += dy[k];
			if (i2 < 0 || i2 >= r || j2 < 0 || j2 >= s) break;
			if (grid[i2][j2] == '#') break;
			if (grid[i2][j2] == 'n'){
			    chosen[i2][j2][dk[k]] = tower[i][j] + k;
			}
			if (grid[i2][j2] == 'T'){
			    sat.add_false(VAR(tower[i][j] + k));
			    break;
			}
		    }
		}
	    } 
	}
    for (int i = 0; i < r; i++){
	for (int j = 0; j < s; j++){
	    if (grid[i][j] != 'n') continue;
	    vector<bool> b(4);
	    for (int k = 0; k < 4; k++) b[k] = (chosen[i][j][k] != -1);
	    if (b[0] ^ b[3]){
		if (b[1] ^ b[2]){
		    int x, y;
		    if (b[0]) x = chosen[i][j][0];
		    else x = chosen[i][j][3];
		    if (b[1]) y = chosen[i][j][1];
		    else y = chosen[i][j][2];
		    sat.add_clause(VAR(x), VAR(y));
		} else {
		    int x;
		    if (b[0]) x = chosen[i][j][0];
		    else x = chosen[i][j][3];
		    sat.add_true(VAR(x));
		}
	    } else {
		if (b[1] ^ b[2]){
		    int y;
		    if (b[1]) y = chosen[i][j][1];
		    else y = chosen[i][j][2];
		    sat.add_true(VAR(y));		
		} else {
		    assert(false);
		}	    
	    }
	}
    }
    auto ans = sat.solve();
    assert(ans.X);
    for (int i = 0; i < r; i++){
	for (int j = 0; j < s; j++){
	    if (grid[i][j] != 'T') continue;
	    vi val(4, 0);
	    for (int k = 0; k < 4; k++){
		val[k] = ans.Y[tower[i][j] + k];
	    }
	    if (val[0]){
		if (val[1]) grid[i][j] = '4';
		else grid[i][j] = '3';
	    } else {
		if (val[1]) grid[i][j] = '1';
		else grid[i][j] = '2';
	    }
	}
    }
    for (int i = 0; i < r; i++){
	for (int j = 0; j < s; j++)
	    cout << grid[i][j];
	cout << endl;
    }

    return 0;
}
