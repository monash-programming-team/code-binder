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
  int cs = 1;
  while (cin >> n){
    vector<string> names(n);
    map<string, int> idx;
    for (int i = 0; i < n; i++) {
      cin >> names[i];
      idx[names[i]] = i;
    }
    cin >> m;
    vvi adj(n);
    for (int i = 0; i < m; i++){
      string s1, s2; cin >> s1 >> s2;
      adj[idx[s1]].push_back(idx[s2]);
    }
    auto p = topological_sort(adj);
    cout << "Case #" << cs++ << ": Dilbert should drink beverages in this order:";
    for (int x : p.second){
      cout << ' ' << names[x];
    }
    cout << "." << endl;
    cout << endl;
  }



  return 0;
}
