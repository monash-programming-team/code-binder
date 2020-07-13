#include <bits/stdc++.h>
using namespace std;

#define MAXN 100010
typedef vector<int> vi;
typedef vector<vi> vvi;

// ==================== BEGIN Centroid Decomposition =====================
vi adj[MAXN], cadj[MAXN];
int marked[MAXN], par[MAXN], vis[MAXN], sz[MAXN];
int cnt = 0;

void dfs(int u){
    sz[u] = 1;
    for (int v : adj[u]){
	if (vis[v] == cnt || marked[v]) continue;
	vis[v] = cnt;
	dfs(v);
	sz[u] += sz[v];
    }
}

int getCentroid(int src, int n){
    for (int v : adj[src]){
	if (vis[v] == cnt || marked[v]) continue;
	if (sz[v] * 2 > n) {
	    vis[v] = cnt;
	    return getCentroid(v, n);
	}
    }
    return src;
}

int getCentroid(int src){
    vis[src] = ++cnt;
    dfs(src);
    vis[src] = ++cnt;
    int centroid = getCentroid(src, sz[src]);
    marked[centroid] = 1;
    return centroid;
}

int decomposeTree(int root){
    int cend_tree = getCentroid(root);
    for (int v : adj[cend_tree]){
	if (marked[v]) continue;
	int cend_subtree = decomposeTree(v);
	cadj[cend_tree].push_back(cend_subtree);
	cadj[cend_subtree].push_back(cend_tree);
	par[cend_subtree] = cend_tree;
    }
    return cend_tree;
}
// ==================== END Centroid Decomposition =====================

int main(){
    memset(par, -1, sizeof par);
    decomposeTree(0);


}
