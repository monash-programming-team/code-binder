// Tested on: UVA 247
//            SPOJ CAPCITY
//            CF 427C
//            LightOJ 1034
//            
#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef pair<int, int> pii;
typedef vector<vi> vvi;

class SCC {
private:
    int n, comp;
    vvi g, gt;
    vi seq;
    void dfs(int u, const vvi &adj){
	for (int i = 0; i < (int) adj[u].size(); i++){
	    int v = adj[u][i];
            if (vis[v] == -1){
                vis[v] = comp;
                dfs(v, adj);
            }
        }
        seq.push_back(u);
    }
public:
    vi vis; // stores component for each node
    SCC(int _n){
        n = _n;
        g.assign(n, vi()); gt.assign(n, vi());
    }
    void add_edge(int u, int v){
        g[u].push_back(v); gt[v].push_back(u);
    }
    int find_SCC(){
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
        return comp;
    }
    // OPTIONAL: find_SCC() must be called first
    vvi get_dag(){
        map<pii, int> mmap;
        vvi dag(comp, vi());
        for (int u = 0; u < n; u++){

            for (int i = 0; i < (int) g[u].size(); i++){
		int v = g[u][i];
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
vi vis;
void dfs(const vvi &dag, int u){
    for (int i = 0; i < (int) dag[u].size(); i++){
	int v = dag[u][i];
	if (!vis[v]){
	    vis[v] = 1;
	    dfs(dag, v);
	}
    }
}
int main(){
    int tc; cin >> tc;
    for (int cs = 1; cs <= tc; cs++){
	int n, m; cin >> n >> m;
	SCC scc(n);
	for (int i = 0; i < m; i++){
	    int u, v; cin >> u >> v;
	    u--; v--;
	    scc.add_edge(u, v);
	}
	int comps = scc.find_SCC();
	int res = 0;
	vis.assign(comps, 0);
	vvi dag = scc.get_dag();
	for (int i = 0; i < comps; i++){
	    if (!vis[i]){
		vis[i] = 1;
		dfs(dag, i);
		res++;
	    }
	}
	cout << "Case " << cs << ": " << res << endl;
    }
    return 0;
}
