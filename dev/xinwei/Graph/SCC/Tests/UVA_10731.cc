// Strongly Connected Components (Kosaraju's)
//
// Author      : Xin Wei Chow
// Date        : Jun 15, 2016
// Reliability : 5
// Tested on   : UVA 247
//               SPOJ CAPCITY
//               CF 427C
//               LightOJ 1034
//
// Finds the number of SCCs
// OPTIONAL: returns the DAG of SCCs
//
// Complexity: O(V + E) - find_SCC()
//             O(E log E) - get_dag()
//
// Return: 
//    int find_SCC():
//        Returns the number of SCCs
//        Component of each node is stored in 'vis'
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
    int n; 
    bool first = true;
    while (cin >> n && n){
	vi used(26, 0);
	if (first) first = false;
	else cout << endl;
	SCC scc(26);
	for (int i = 0; i < n; i++){
	    vi v(5);
	    for (int j = 0; j < 5; j++){
		char c; cin >> c;
		v[j] = (int) c - 'A';
		used[v[j]] = 1;
	    }
	    char c; cin >> c;
	    int ans = (int) c - 'A';
	    for (int j = 0; j < 5; j++){
		if (v[j] != ans){
		    scc.add_edge(ans, v[j]);
		}
	    }
	}

	int comps;
	vi cmp;
	tie(comps, cmp) = scc.find_SCC();
	vector<vi> res(comps, vi());
	for (int i = 0; i < 26; i++){
	    res[cmp[i]].push_back(i);
	}
	for (int i = 0; i < comps; i++)
	    sort(res[i].begin(), res[i].end());
	sort(res.begin(), res.end());
	for (int i = 0; i < comps; i++){
	    //if (res[i].size() == 0) continue;
	    if (used[res[i][0]] == 0) continue;
	    for (int j = 0; j < (int) res[i].size(); j++){
		char c = (char) res[i][j] + 'A';
		cout << c;
		if (j != (int) res[i].size() - 1) cout << ' ';
	    }
	    cout << endl;
	}
    }

    return 0;
}
