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
  priority_queue<int, vector<int>, greater<int>> q;
  for (int i = 0; i < n; i++) 
    if (indeg[i] == 0) q.push(i);
  while (!q.empty()){
    int u = q.top(); q.pop();
    order.push_back(u);
    for (int v : adj[u]){
      if (--indeg[v] == 0) q.push(v);
    }
  }
  return { (int) order.size() == n, order };
}

int main(){
  int n, m;
  cin >> n >> m;
  vvi adj(n);
  for (int i = 0; i < m; i++){
    int target; cin >> target; target--;
    int k; cin >> k;
    for (int j = 0; j < k; j++){
      int v; cin >> v; v--;
      adj[v].push_back(target);
    }
  }
  auto p = topological_sort(adj);
  for (int x : p.second){
    cout << x + 1 << ' ';
  }
  cout << endl;



  return 0;
}
