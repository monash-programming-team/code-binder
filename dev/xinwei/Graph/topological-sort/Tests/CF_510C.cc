// UVA - 10305
// UVA - 11060
// SPOJ - PFDEP

#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int, int> pii;

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

pii compare(string &A, string &B){
  int len = min((int) A.size(), (int) B.size());
  for(int i = 0; i < len; i++){
    if (A[i] == B[i]) continue;
    return { (int) A[i] - 'a', (int) B[i] - 'a' };
  }
  if ((int) A.size() <= (int) B.size()) return {-1, 1};
  return {-1, -1}; // invalid
}

int main(){
  int n; cin >> n;
  vvi adj(26);
  vector<string> arr(n);
  for(int i = 0; i < n; i++) cin >> arr[i];
  for (int i = 0; i < n; i++){
    for (int j = i + 1; j < n; j++){
      if (arr[i] == arr[j]) continue;
      pii p = compare(arr[i], arr[j]);
      if (p.first != -1){
	adj[p.first].push_back(p.second);
      } else {
	if (p.second == -1) {
	  cout << "Impossible" << endl;
	  return 0;
	}
      }
    }
  }
  auto ans = topological_sort(adj);
  if (ans.first == false) {
    cout << "Impossible" << endl;
    return 0;
  }
  for (int i = 0; i < 26; i++){
    int c = ans.second[i];
    char cc = (char) c + 'a';
    cout << cc;
  }
  cout << endl;
  return 0;
}
