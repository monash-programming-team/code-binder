// Strongly Connected Components (Kosaraju's)
//
// Author      : Xin Wei Chow
// Date        : Jun 15, 2016
// Reliability : 5
// Tested on   : UVA 247
//               SPOJ CAPCITY
//               CF 427C
//               LightOJ 1034
//               UVA 10731
//
// Finds the number of SCCs
// OPTIONAL: returns the DAG of SCCs
//
// Complexity: O(V + E) - find_SCC()
//             O(E log E) - get_dag()
//
// Return:
//    pair<int, vi> find_SCC():
//        Returns the number of SCCs
//        and vector<int> containing component number of each node
//    vvi get_dag():
//        Returns the DAG of contracted SCCs

#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef pair<int, int> pii;
typedef vector<vi> vvi;

class SCC {
private:
    int n, comp_sz;
    vvi g, gt;
    vi seq, vis;
    void dfs(int u, const vvi &adj){
	for (int v : adj[u]){
            if (vis[v] == -1){
                vis[v] = comp_sz;
                dfs(v, adj);
            }
        }
        seq.push_back(u);
    }
public:
    SCC(int _n){
        n = _n;
        g.assign(n, vi()); gt.assign(n, vi());
    }
    void add_edge(int u, int v){
        g[u].push_back(v); gt[v].push_back(u);
    }
    pair<int, vi> find_SCC(){
        vis.assign(n, -1); comp_sz = 0;
        for (int i = 0; i < n; i++){
            if (vis[i] == -1){
                vis[i] = comp_sz;
                dfs(i, g);
            }
        }
        vis.assign(n, -1); comp_sz = 0;
        for (int i = n-1; i >= 0; i--){
            int u = seq[i];
            if (vis[u] == -1){
                vis[u] = comp_sz;
                dfs(u, gt);
                comp_sz++;
            }
        }
        return {comp_sz, vis};
    }
    // OPTIONAL: find_SCC() must be called first
    vvi get_dag(){
        map<pii, int> mmap;
        vvi dag(comp_sz, vi());
        for (int u = 0; u < n; u++){
            for (int v : g[u]){
                if (vis[u] == vis[v]) continue;
                if (!mmap.count(pii(vis[u], vis[v]))){
                    dag[vis[u]].push_back(vis[v]);
                    dag[vis[v]].push_back(vis[u]); // modification to base code
                    mmap[pii(vis[u], vis[v])] = 1;
                }
            }
        }
        return dag;
    }
};
vi vis, comp_sz;
bool large = false;
int dfs(const vvi &dag, int u){
    int sz = 0;
    if (comp_sz[u] == 1) sz = 1;
    else if (large == false){
        large = true;
        sz = 1;
    }
    for (int v : dag[u]){
        if (!vis[v]){
            vis[v] = 1;
            sz += dfs(dag, v);
        }
    }
    return sz;
}

int main(){
    int n, m; cin >> n >> m;
    SCC scc(n);
    for (int i = 0; i < m; i++){
	int u, v; cin >> u >> v;
	u--; v--;
	scc.add_edge(u, v);
    }
    int comps;
    vi cmp;
    tie(comps, cmp) = scc.find_SCC();
    comp_sz.assign(comps, 0);
    for (int i = 0; i < n; i++){
	comp_sz[cmp[i]]++;
    }
    int res = 0;
    for (int i = 0; i < comps; i++){
	res += ((comp_sz[i] > 1) ? comp_sz[i] : 0);
    }
    vvi dag = scc.get_dag();
    vis.assign(comps, 0);
    for (int i = 0; i < comps; i++){
    	if (!vis[i]){
            large = false;
    	    vis[i] = 1;
    	    int sz = dfs(dag, i);
    	    res += sz - 1;
    	}
    }
    cout << res << endl;


    return 0;
}
