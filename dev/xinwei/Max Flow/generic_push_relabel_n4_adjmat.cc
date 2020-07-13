#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const ll INF = LLONG_MAX;
vector<vector<ll> > c, f;
vector<ll> e;
vector<int> h;
int n, s, t;
void push(int u, int v){
  ll r = min(e[u], c[u][v] - f[u][v]);
//  cout << u << ' ' << v << ' ' << r << endl;
  e[u] -= r;
  e[v] += r;
  f[u][v] += r;
  f[v][u] -= r;
}

ll push_relabel(int _n, int _s, int _t){
  n = _n; s = _s; t = _t;
  // initialize
  e.assign(n, 0);
  h.assign(n, 0); h[s] = n;

  for (int i = 0; i < n; i++){
    if (i == s) continue;
    if (c[s][i] == 0) continue;
    e[i] += c[s][i];
    f[s][i] += c[s][i];
    f[i][s] -= c[s][i];
  }

  while (1){
    bool change = false;
    for (int i = 0; i < n; i++){
      if (i == s || i == t) continue;
      if (e[i] == 0) continue;
      // overflowing
      change = true;
      int minH = 2 * n;
      bool pushed = false;
      for (int j = 0; j < n; j++){
	if (c[i][j] == f[i][j]) continue;
	if (h[i] == h[j] + 1){
	  push(i, j);
	  pushed = true;
	  break;
	}
	minH = min(minH, h[j]);
      }
      if (pushed == false){
	assert(minH != 2 * n);
	h[i] = 1 + minH;
      }
    }
    if (!change) break;
  }
  return e[t];
}

int main(){
  int n,m,s,t;
  scanf("%d%d%d%d",&n,&m,&s,&t);
  c.assign(n, vector<ll> (n, 0));
  f.assign(n, vector<ll> (n, 0));

  for(int i=0;i<m;i++){
    int u,v,cap;
    scanf("%d%d%d",&u,&v,&cap); u--; v--;
    if (u == v) continue;
    c[u][v] += cap;
  }


  ll mf = push_relabel(n, s-1, t-1);
  cout << mf << endl;
  return 0;
}
