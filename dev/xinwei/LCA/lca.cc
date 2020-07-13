#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef vector<vi> vvi;
const int MAXN = 100010;
const int LOGN = 20;
vi adj[MAXN];

class LCA {
private:
    int n;
    vvi anc;
    vi par, lvl;
    void dfs(int u, int p, int depth){
	par[u] = p;
	lvl[u] = depth;
	for (int v : adj[u]){
	    if (v != p)
		dfs(v, u, depth + 1);
	}
    }

public:
    LCA(int n, int rt = 0) : n(n), par(n), lvl(n) {
	dfs(rt, -1, 0);
	anc.assign(n, vi(LOGN, -1));
	for (int i = 0; i < n; i++) anc[i][0] = par[i];
	for (int k = 1; k < LOGN; k++){
	    for (int i = 0; i < n; i++){
		if (anc[i][k-1] != -1) anc[i][k] = anc[anc[i][k-1]][k-1];
	    }
	}
    }

    int query(int u, int v){
	if (lvl[u] > lvl[v]) swap(u, v);
	for (int k = LOGN - 1; k >= 0; k--){
	    if (lvl[v] - (1 << k) >= lvl[u])
		v = anc[v][k];
	}
	if (u == v) return u;
	for (int k = LOGN - 1; k >= 0; k--){
	    if (anc[u][k] == anc[v][k]) continue;
	    u = anc[u][k]; v = anc[v][k];
	}
	return par[u];
    }

    int dist(int u, int v){
	int p = query(u, v);
	return lvl[u] + lvl[v] - 2 * lvl[p];
    }
};

int main(){


    return 0;
}
