#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef vector<vi> vvi;

// Topological Sort
//
// Author      : Xin Wei Chow
// Date        : 29 Jun 2016
// Reliability : 5
// Tested On   : UVA - 10305
//               UVA - 11060
//               SPOJ - PFDEP
//               UVA - 200
//               CF - 510C
// 
// Finds and returns a valid topological ordering.
//
// Complexity: O( V + E )
// 
// Return:
//  pair of bool and vector<int>. The first is 'bool'. The second is order.
//  'bool' :  returns true if a valid ordering exists
//  order  :  order contains the index of vertices in topological order
// 
// Parameters:
//  adj    : The directed graph
//
// Note: if a specific ordering is needed (e.g. lexicographically least),
//       queue will have to be changed to priority_queue
//       Complexity: O( V log V + E )

//listings:topsort
// Returns whether a topological ordering exists and a valid topological
// ordering if one does. Complexity: O(V + E)
pair<bool, vi> topological_sort(const vvi &adj) {
  int n = (int) adj.size();  vi indeg(n, 0), order;
  for (int i = 0; i < n; i++) for (int v : adj[i]) indeg[v]++;
  queue<int> q; // priority_queue can be used here if you want tie-breaking
  for (int i = 0; i < n; i++) if (indeg[i] == 0) q.push(i);
  while (!q.empty()){
    int u = q.front(); q.pop(); // q.front()->q.top() if priority_queue
    order.push_back(u);
    for (int v : adj[u]) if (--indeg[v] == 0) q.push(v);
  }
  return { (int)order.size() == n, order };
}
//listings:/topsort

int main(){
  int n = 4;
  vvi adj(n);
  adj[0].push_back(1);
  adj[0].push_back(2);
  adj[1].push_back(3);
  adj[3].push_back(2);
  auto p = topological_sort(adj);
  if (p.first){
    cout << "VALID TOPOLOGICAL ORDERING: " << endl;
    // 0 1 3 2
    for (int x : p.second) {
      cout << x << ' '; 
    }
    cout << endl;
  } else 
    cout << "NO VALID TOPOLOGICAL ORDERING" << endl;

  return 0;
}
