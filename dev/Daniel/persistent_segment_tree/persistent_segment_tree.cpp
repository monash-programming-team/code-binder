// Persistent segment tree
//
// Author: Daniel Anderson
// Reliability: 3

#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> pii;

// Memory-safe Persistent Segment Tree
// T must define:
//   typedef typename type: the value type
//   type merge(type a, type b): merge two child values
//   type I(): an identity value for the merge operation
template<typename T> struct PersistentSegmentTree {
  typedef typename T::type value_type;
  struct Node {
    value_type val;  int lo, hi;
    shared_ptr<Node> left, right;
    Node(value_type val, int lo, int hi, shared_ptr<Node> left, shared_ptr<Node> right)
      : val(val), lo(lo), hi(hi), left(move(left)), right(move(right)) { }
  };
  int n;  shared_ptr<Node> root;
  PersistentSegmentTree(const vector<value_type>& A) : n(A.size()) {
    root = build(A, 0, n-1);
  }
  shared_ptr<Node> build(const vector<value_type>& A, int lo, int hi) {
    if (lo > hi)
      return nullptr;
    if (lo == hi)
      return make_shared<Node>(A[lo], lo, hi,nullptr, nullptr);
    else {
        int mid = (lo + hi)/2;
        auto left = build(A, lo, mid);
        auto right = build(A, mid+1, hi);
        value_type val = T::merge(left ? left->val : T::I(), right ? right->val : T::I());
        return make_shared<Node>(val, lo, hi, move(left), move(right));
      }
  }
  void update(int i, value_type val) {
    root = update(i, val, root);
  }
  value_type query(int L, int R) {
    return query(L, R, *root);
  }
  value_type query(int L, int R, Node& node) {
    if (node.lo >= L && node.hi <= R) return node.val;
    else if (node.lo > R || node.hi < L) return T::I();
    else return T::merge(node.left ? query(L,R,*node.left) : T::I(), node.right ? query(L,R,*node.right) : T::I());
  }
  shared_ptr<Node> update(int i, value_type val, const shared_ptr<Node>& node) {
    if (node->lo > i || node->hi < i) return node;
    if (node->lo == i && node->hi == i) {
      return make_shared<Node>(val, i, i, nullptr, nullptr);
    } else {
      shared_ptr<Node> left = node->left ? update(i, val, node->left) : nullptr;
      shared_ptr<Node> right = node->right ? update(i, val, node->right) : nullptr;
      value_type merged = T::merge(left ? left->val : T::I(), right ? right->val : T::I());
      return make_shared<Node>(merged, node->lo, node->hi, move(left), move(right));
    }
  }
};

// Persistent segment tree traits for ranged sum with point-value assignment.
template<typename T> struct SumAssign {
    typedef T type;
    static T merge(T a, T b) { return a + b; }
    static T I() { return 0; };
};

namespace SPOJ_MKTHNUM {
  void solve() {
    int n, m;
    cin >> n >> m;
    vi a(n);
    for (auto& x : a) cin >> x;
    
    vi keys = a;
    sort(keys.begin(), keys.end());
    keys.erase(unique(keys.begin(), keys.end()), keys.end());
    
    int maxa = (int)keys.size();
    
    PersistentSegmentTree<SumAssign<int>> tree(vi(maxa, 0));
    auto dummy = tree;
    
    vector<PersistentSegmentTree<SumAssign<int>>> history;
    
    auto key = [&](int x) { 
      return lower_bound(keys.begin(), keys.end(), x) - keys.begin();
    };
    
    for (int i=0; i<n; i++) {
      int c = key(a[i]);
      int cur = tree.query(c,c);
      tree.update(c, cur+1);
      history.push_back(tree);
    }
    
    auto query = [&](int i, int j, int k) {
      auto node1 = history[j].root;
      auto node2 = i ? history[i-1].root : dummy.root;
      
      while (node1->left) {
        if (node1->left->val - node2->left->val < k) {
          k -= (node1->left->val - node2->left->val);
          node1 = node1->right;
          node2 = node2->right;
        } else {
          node1 = node1->left;
          node2 = node2->left;
        }
      }
      
      return keys[node1->lo];
    };
    
    while (m--) {
      int i, j, k;
      cin >> i >> j >> k;
      i--, j--;
      cout << query(i,j,k) << '\n';
    }
  }
}

namespace SPOJ_DQUERY {
  const int MAXA = 1000000;
  
  void solve() {
    // Read input
    int n;
    scanf("%d", &n);
    vi a(n);
    for (auto& x : a) scanf("%d", &x);
    
    // Preprocess
    PersistentSegmentTree<SumAssign<int>> tree(vi(MAXA + 1, 0));
    vi last_pos(MAXA + 1, -1);
    
    vector<PersistentSegmentTree<SumAssign<int>>> history;
    
    for (int i=0; i<n; i++) {
      if (last_pos[a[i]] != -1) tree.update(last_pos[a[i]], 0);
      last_pos[a[i]] = i;
      tree.update(i, 1);
      history.push_back(tree);
    }
    
    // Solve one query
    auto query = [&](int i, int j) {
      return history[j].query(i, j);
    };
    
    // Process queries
    int q;
    scanf("%d", &q);
    while (q--) {
      int i, j;
      scanf("%d %d", &i, &j);
      i--, j--;
      printf("%d\n", query(i,j));
    }
  }
}

namespace SPOJ_COT {

  template<typename T = int> struct LCA {
	  const int LOGN = 20;	// works for n <= 10^6. Change appropriately.
	  int n;  vi par, lvl;  vvi anc;  vector<T> len;  vector<vector<pair<int,T>>> adj;
	  void dfs(int u, int p, int l, T d) {
		  par[u] = p, lvl[u] = l, len[u] = d;
		  for (auto v : adj[u]) if (v.X != p) dfs(v.X, u, l+1, d+v.Y);
	  }
	  // Create a tree with n nodes. Add edges then call build(root).
	  LCA(int n) : n(n), par(n), lvl(n), len(n), adj(n) { }
    void add_edge(int u, int v, T w = 1) {
      adj[u].emplace_back(v, w), adj[v].emplace_back(u, w);
    }
    void build(int root = 0) {  // Call this before making queries
      dfs(root,-1,0,0), anc.assign(n, vi(LOGN, -1));
		  for (int i = 0; i < n; i++) anc[i][0] = par[i];
		  for (int k = 1; k < LOGN; k++) for (int i = 0; i < n; i++)
			  if (anc[i][k-1] != -1) anc[i][k] = anc[anc[i][k-1]][k-1];
    }
	  int query(int u, int v) {  // LCA with respect to original root
		  if (lvl[u] > lvl[v]) swap(u, v);
		  for (int k = LOGN - 1; k >= 0; k--) 
			  if (lvl[v] - (1 << k) >= lvl[u]) v = anc[v][k];
		  if (u == v) return u;
		  for (int k = LOGN - 1; k >= 0; k--) {
			  if (anc[u][k] == anc[v][k]) continue;
			  u = anc[u][k]; v = anc[v][k];
		  }
		  return par[u];
	  }
    int query(int u, int v, int root) {  // OPTIONAL: LCA with respect to any root
      int a = query(u, v), b = query(u, root), c = query(v, root);
      if (a == c && c != b) return b;
      else if (a == b && c != b) return c;
      else return a;
    }
	  T dist(int u, int v) { return len[u] + len[v] - 2 * len[query(u,v)]; }
  };

  typedef PersistentSegmentTree<SumAssign<int>> pst;

  void solve() {
    int N, M;
    cin >> N >> M;
    vector<ll> w(N);
    for (auto& x : w) cin >> x;
    
    auto keys = w;
    sort(keys.begin(), keys.end());
    keys.erase(unique(keys.begin(), keys.end()), keys.end());
    int num_keys = (int)keys.size();
    
    auto key = [&](ll x) {
      return lower_bound(keys.begin(), keys.end(), x) - keys.begin();
    };
    
    LCA<int> lca(N);
    for (int i=0; i<N-1; i++) {
      int u, v;
      cin >> u >> v;
      u--, v--;
      lca.add_edge(u, v);
    }
    lca.build();
    
    pst tree(vi(num_keys, 0));
    pst dummy = tree;
    vector<pst> history(N, pst(vi()));
   
    function<void(int,int,pst)> dfs = [&](const int u, const int p, pst cur_tree) {
      int c = key(w[u]);
      int cur = cur_tree.query(c,c);
      cur_tree.update(c, cur+1);
      history[u] = cur_tree;
      for (const auto& v : lca.adj[u]) if (v.X != p) {
        dfs(v.X, u, cur_tree);
      }
    };
    dfs(0, -1, tree);
   
    auto query = [&](int u, int v, int k) {
      int anc = lca.query(u, v);
      
      auto an = history[u].root, bn = history[v].root;
      auto ln = history[anc].root;
      auto pln = (lca.par[anc] == -1 ? dummy.root : history[lca.par[anc]].root);
    
      while (an->left) {
        int total = an->left->val + bn->left->val - ln->left->val - pln->left->val;
        if (total >= k) {
          an = an->left;
          bn = bn->left;
          ln = ln->left;
          pln = pln->left;
        } else {
          an = an->right;
          bn = bn->right;
          ln = ln->right;
          pln = pln->right;
          k -= total;
        }
      }
    
      return keys[an->lo];
    };
    
    while (M--) {
      int u, v, k;
      cin >> u >> v >> k;
      u--, v--;
      cout << query(u, v, k) << '\n';
    }
  }
}

int main() {
  ios::sync_with_stdio(0); cin.tie(0);
  SPOJ_COT::solve();
}


