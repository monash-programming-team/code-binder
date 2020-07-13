#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef vector<vi> vvi;

pair<bool, vi> topological_sort(const vvi &adj) {
  int n = (int) adj.size();
  vi indeg(n, 0), order;
  for (int i = 0; i < n; i++){
    for (int v : adj[i]) indeg[v]++;
  }
  queue<int> q;
  for (int i = 0; i < n; i++) 
    if (indeg[i] == 0) q.push(i);
  while (!q.empty()){
    int u = q.front(); q.pop();
    order.push_back(u);
    for (int v : adj[u]){
      if (--indeg[v] == 0) q.push(v);
    }
  }
  return { (int) order.size() == n, order };
}

int main(){
  int n, m;
  while (cin >> n >> m && n){
    vector<vi> adj(n, vi());
    for (int i = 0; i < m; i++){
      int u, v; cin >> u >> v;
      u--; v--;
      adj[u].push_back(v);
    }
    auto ans = topological_sort(adj);
    for (int i = 0; i < n; i++){
      cout << ans.second[i] + 1;
      if (i != n-1) cout << ' ';
    }
    cout << endl;
  }


  return 0;
}
