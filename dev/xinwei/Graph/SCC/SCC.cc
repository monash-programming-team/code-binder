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
//               CF 505D
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
    // OPTIONAL: find_SCC() must be called first
    vvi get_dag(){
        map<pii, int> mmap;
        vvi dag(comp, vi());
        for (int u = 0; u < n; u++){
            for (int v : g[u]){
                if (vis[u] == vis[v]) continue;
                if (!mmap.count(pii(vis[u], vis[v]))){
                    dag[vis[u]].push_back(vis[v]);
                    mmap[pii(vis[u], vis[v])] = 1;
                }
            }
        }
        return dag;
    }
};

int main(){
    int n = 8;
    SCC scc(n);
    vector<pii> edges = {{0, 1}, {1, 3}, {2, 1}, {3, 2}, {3, 4}, {4, 5},
                         {5, 7}, {6, 4}, {7, 6}};
    for (pii p : edges)
        scc.add_edge(p.first, p.second);

    int comps;
    tie(comps, ignore) = scc.find_SCC();
    cout << comps << endl;
    vvi g = scc.get_dag();
    for (int u = 0; u < comps; u++){
        for (int v : g[u]){
            cout << u << "-->" << v << endl;
        }
    }

    return 0;
}
