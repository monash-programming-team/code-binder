#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
int n, m;
vi arr;

typedef vector<vi> vvi;
template<typename T>
class RMQ {
  typedef pair<T, int> pti;
private:
  int n;
  vector<T> v;
  vector<pti> st;

  pti merge(pti p1, pti p2) { 
    if (p1.second == -1) return p2;
    if (p2.second == -1) return p1;
    return max(p1, p2); // change this to max for max query
  }
  void build(int p, int L, int R){
    if (L == R) st[p] = {v[L], L};
    else {
      int mid = (L + R) / 2;
      build(p * 2, L, mid);
      build(p * 2 + 1, mid + 1, R);
      st[p] = merge(st[p * 2], st[p * 2 + 1]);
    }
  }
  pti query(int p, int L, int R, int i, int j){
    if (i > R || j < L) return {-1, -1};
    if (i <= L && j >= R) return st[p];
    int mid = (L + R) / 2;
    pti p1 = query(p * 2, L, mid, i, j);
    pti p2 = query(p * 2 + 1, mid + 1, R, i, j);
    return merge(p1, p2);
  }
  void update(int p, int L, int R, int pos, T val){
    if (pos > R || pos < L) return;
    if (pos == L && pos == R) st[p].first = val;
    else {
      int mid = (L + R) / 2;
      update(p * 2, L, mid, pos, val);
      update(p * 2 + 1, mid + 1, R, pos, val);
      st[p] = merge(st[p * 2], st[p * 2 + 1]);
    }
  }

public:
  RMQ(vector<T> &v) : v(v) {
    n = v.size();
    st.resize(4 * n);
    build(1, 0, n-1);
  }
  pti query(int i, int j){ return query(1, 0, n-1, i, j); }
  void update(int pos, T val){ update(1, 0, n-1, pos, val); }
};

int main(){
  std::ios_base::sync_with_stdio(false);
  scanf("%d%d", &n, &m);  
	
  arr.resize(n);
  for (int i = 0; i < n; i++){
    scanf("%d", &arr[i]);
  }
  RMQ<int> cc(arr);
	
  while (m--){
    int type; scanf("%d", &type);
    if (type == 1){
      int id; scanf("%d", &id);
      int ans1 = -1, ans2 = n;
      int L = 0, R = id - 1;
      if (cc.query(L, R).first >= arr[id]){
	while (L < R){
	  int mid = (L + R) / 2;
	  int v = cc.query(mid, id - 1).first;
	  if (v >= arr[id]) L = mid + 1;
	  else R = mid;
	}
	if (arr[L] >= arr[id]) ans1 = L;
	else ans1 = L - 1;
      }
      L = id + 1, R = n - 1;
      if (cc.query(L, R).first >= arr[id]){
	while (L < R){
	  int mid = (L + R) / 2;
	  int v = cc.query(id + 1, mid).first;
	  if (v < arr[id]) L = mid + 1;
	  else R = mid;
	}
	if (arr[L] >= arr[id]) ans2 = L;
	else ans2 = L + 1;
      }
      cout << ans1 << ' ' << ans2 << "\n";
    } else {
      int id, val; scanf("%d%d", &id, &val);
      arr[id] = val;
      cc.update(id, val);
    } 
  }

  return 0;
}

