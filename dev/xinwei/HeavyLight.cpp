// Heavy Light Decomposition
// Solves QTREE - SPOJ 
// http://www.spoj.com/problems/QTREE/
// Complexity: O(n + q log log n)
#include <bits/stdc++.h>

using namespace std;
const int INF = int(1e9);
const int MAXN = 100010;
typedef vector<int> vi;
typedef pair<int, int> pii;

class SegmentTree {
private:
    int n;
    vi arr, st, lazy;
    int merge(int p1, int p2){
	return max(p1, p2);
    }
    void build(int p, int L, int R){
	if (L == R) st[p] = arr[L];
	else {
	    int mid = (L + R) / 2;
	    build(p * 2, L, mid);
	    build(p * 2 + 1, mid + 1, R);
	    st[p] = merge(st[p * 2], st[p * 2 + 1]);
	}
    }
    int query(int p, int L, int R, int i, int j){
	if (i > R || j < L) return -INF;
	if (i <= L && j >= R) return st[p];
	int mid = (L + R) / 2;
	int p1 = query(p * 2, L, mid, i, j);
	int p2 = query(p * 2 + 1, mid + 1, R, i, j);
	return merge(p1, p2);
    }
    void update(int p, int L, int R, int pos, int val){
	if (pos > R || pos < L) return;
	if (pos == L && pos == R) st[p] = val;
	else {
	    int mid = (L + R) / 2;
	    update(p * 2, L, mid, pos, val);
	    update(p * 2 + 1, mid + 1, R, pos, val);
	    st[p] = merge(st[p * 2], st[p * 2 + 1]);
	}
    }

public:
    SegmentTree(int n, const vi &arr) : n(n), arr(arr) {
	st.assign(4 * n, 0);
	lazy.assign(4 * n, 0);
	build(1, 0, n-1);
    }
    int query(int i, int j){
	return query(1, 0, n-1, i, j);
    }
    void update(int pos, int val){
	update(1, 0, n-1, pos, val);
    }
};

int n;
int heavy[MAXN], head[MAXN], par[MAXN], pos[MAXN], level[MAXN], cost[MAXN];
vector<pii> adj[MAXN];

int dfs(int u, int p = -1, int d = 0){
  int size = 1, maxSubTree = 0, maxSubTree_id;
  par[u] = p;
  level[u] = d;
  for (int id = 0; id < (int) adj[u].size(); id++){
    pii v = adj[u][id];
    if (v.first != p){
      cost[v.first] = v.second;
      int subTree = dfs(v.first, u, d + 1);
      if (subTree > maxSubTree){
	maxSubTree = subTree;
	maxSubTree_id = v.first;
      }
      size += subTree;
    }
  }
  if (maxSubTree * 2 >= size){
    heavy[u] = maxSubTree_id;
  }
  return size;
}

void HeavyLight(vi &arr, int root = 0){
  memset(heavy, -1, sizeof heavy);
  dfs(root);
  arr.resize(n);
  // re-number vertices
  int curPos = 0;
  for (int i = 0; i < n; i++){
    if (par[i] == -1 || heavy[par[i]] != i){ // find head of heavy path
      int cur = i;
      while (cur != -1) {
	arr[curPos] = cost[cur];
	pos[cur] = curPos++;
	head[cur] = i;
	cur = heavy[cur];
      }
    }
  }
}

int solve(int u, int v, SegmentTree &cc){
  int res = -INF;
  while (head[u] != head[v]){
    if (level[head[u]] > level[head[v]]) swap(u, v);
    res = max(res, cc.query(pos[head[v]], pos[v]));
    v = par[head[v]];
  }
  if (u == v) return res;
  res = max(res, cc.query(min(pos[u], pos[v]) + 1, max(pos[u], pos[v])));
  return res;
}

int main(){
  int t; scanf("%d", &t);
  while (t--){
    scanf("%d", &n);
    for (int i = 0; i < MAXN; i++) adj[i].clear();
    vector<pii> edges(n-1);
    for (int i = 0; i < n-1; i++){
      int u, v, c; 	scanf("%d%d%d", &u, &v, &c);
      u--; v--;
      adj[u].push_back(pii(v, c));
      adj[v].push_back(pii(u, c));
      edges[i] = pii(u, v);
    }

    vi arr;
    HeavyLight(arr);
    SegmentTree cc(n, arr);

    char type[10];
    while (scanf("%s", type) && type[0] != 'D'){
      if (type[0] == 'C'){
	int id, c; scanf("%d%d", &id, &c);
	id--;
	int u = edges[id].first, v = edges[id].second;
	if (level[u] > level[v]) swap(u, v);
	assert(par[v] == u);
	cost[v] = c;
	cc.update(pos[v], c);
	continue;
      } 
      int u, v; scanf("%d%d", &u, &v);
      u--; v--;
      int res = solve(u, v, cc);
      cout << res << "\n";
    }
  }

  return 0;
}