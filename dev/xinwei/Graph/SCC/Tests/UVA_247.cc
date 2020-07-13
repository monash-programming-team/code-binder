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
        for (int v : adj[u]){
            if (vis[v] == -1){
                vis[v] = comp;
                dfs(v, adj);
            }
        }
        seq.push_back(u);
    }
public:
    vi vis;
    SCC(int _n){
        n = _n;
        g.assign(n, vi()); gt.assign(n, vi());
    }
    void add_edge(int u, int v){
        g[u].push_back(v);
        gt[v].push_back(u);
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
    // find_SCC() must be called first
    vvi get_dag(){
        map<pii, int> mmap;
        vvi dag(comp, vi());
        for (int u = 0; u < n; u++){
            for (int v : g[u]){
                if (vis[u] == vis[v]) continue;
                if (!mmap.count({vis[u], vis[v]})){
                    dag[vis[u]].push_back(vis[v]);
                    mmap[{vis[u], vis[v]}] = 1;
                }
            }
        }
        return dag;
    }
};

int main(){
    int n, m, cs = 1;
    while (cin >> n >> m && (n || m)){
        if (cs != 1) cout << endl;
        int id = 0;
        SCC scc(n);
        map<string, int> mmap;
        vector<string> names;
        for (int i = 0; i < m; i++){
            string s1, s2; cin >> s1 >> s2;
            if (!mmap.count(s1)) mmap[s1] = id++, names.push_back(s1);
            if (!mmap.count(s2)) mmap[s2] = id++, names.push_back(s2);
            scc.add_edge(mmap[s1], mmap[s2]);
        }
        cout << "Calling circles for data set " << cs++ << ":" << endl;
        int comps = scc.find_SCC();
        for (int i = 0; i < comps; i++){
            bool first = true;
            for (int j = 0; j < n; j++){
                if (scc.vis[j] == i){
                    if (first) first = false;
                    else cout << ", ";
                    cout << names[j];
                }
            }
            cout << endl;
        }
    }

    return 0;
}
