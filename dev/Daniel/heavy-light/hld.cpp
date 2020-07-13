#include<bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;

//Segment tree implementing range-update & range-max-query by default
struct segment_tree {
  int n;
  vi t, lazy;

  //build tree full of initial values
  void build(int root, int left, int right, vi& a) {
    if (left == right) t[root] = a[left];
    else {
      int mid = (left + right) / 2;
      build(2*root+1, left, mid, a);
      build(2*root+2, mid+1, right, a);
      t[root] = max(t[2*root+1], t[2*root+2]);
    }
  }

  segment_tree() { }

  segment_tree(vi& a) : n(a.size()), t(4*n), lazy(4*n) {
    build(0, 0, n - 1, a);
  }

  int query(int from, int to) {
    if (from > to) swap(from, to);
    return query(from, to, 0, 0, n - 1);
  }

  void update(int from, int to, int delta) {
    update(from, to, delta, 0, 0, n - 1);
  }

  void push(int v) {
    t[2*v+1] += lazy[v];
    t[2*v+2] += lazy[v];
    lazy[2*v+1] += lazy[v];
    lazy[2*v+2] += lazy[v];
    lazy[v] = 0;
  }

  int query(int from, int to, int root, int left, int right) {
    if (from == left && to == right) return t[root];
    push(root);
    int mid = (left + right) / 2;
    int res = INT_MIN;
    if (from <= mid) res = max(res, query(from, min(to, mid), 2*root+1, left, mid));
    else if (to > mid) res = max(res, query(max(from, mid+1), to, 2*root+2, mid+1, right));
    return res;
  }

  void update(int from, int to, int delta, int root, int left, int right) {
    if (from == left && to == right) {
      t[root] += delta;
      lazy[root] += delta;
      return;
    }
    push(root);
    int mid = (left + right) / 2;
    if (from <= mid) update(from, min(to, mid), delta, 2*root+1, left, mid);
    if (to > mid) update(max(from, mid+1), to, delta, 2*root+2, mid+1, right);
    t[root] = max(t[2*root+1], t[2*root+2]);
  }
};

//HLD supporting range-minimum query on the tree
//Small test cases are successful so far
//TODO: Add ranged updates
struct hld {
  int n, root;
  vi size, parent, level, head, len, chainHeads, val;
  map<int, segment_tree> trees;

  //Finds the heaviest child of a node, -1 if it has no children
  int heavy_child(vvi& tree, int cur) {
    int sc = -1; int sc_size = -1;
    for (int i=0; i<tree[cur].size(); i++)
      if (size[tree[cur][i]] > sc_size) {
	sc_size = size[tree[cur][i]];
	sc = i;
      }
    return sc;
  }

  //Build/perform the HLD
  void dfs(vvi& tree, int cur) {
    if (head[cur] == -1) {
      head[cur] = cur;
      if (!tree[cur].empty()) chainHeads.push_back(cur);
    }
    len[head[cur]]++;
    //Find the heaviest child
    int sc = heavy_child(tree, cur);
    //Traverse the heavy child
    if (sc > -1) {
      head[tree[cur][sc]] = head[cur];
      dfs(tree, tree[cur][sc]);
    }
    //Traverse every other child
    for (int i=0; i<tree[cur].size(); i++) {
      if (i != sc) {
	dfs(tree, tree[cur][i]);
      }
    }
  }

  //Precomputes subtree size, parents and levels
  int pre(vvi& tree, int v, int l) {
    level[v] = l;
    if (size[v] > 0) return size[v];
    if (tree[v].empty()) return (size[v] = 1);
    int s = 1;
    for (int c : tree[v]) {
      parent[c] = v;
      s += pre(tree, c, l+1);
    }
    return (size[v] = s);
  }

  //Builds the segment trees
  void make_st(vvi& tree) {
    for (int h : chainHeads) {
      vi chain = {val[h]}; chain.reserve(len[h]);
      for (int i = 0, v = h; i < len[h] - 1; i++)
	chain.push_back(val[(v = tree[v][heavy_child(tree, v)])]);
      trees.insert({h, segment_tree(chain)});
    }
  }

  //Find the minimum weight node on the path from u to v inclusive
  int query(int u, int v) {
    if (level[head[u]] < level[head[v]]) swap(u, v);
    int upos = level[u] - level[head[u]];
    int vpos = level[v] - level[head[v]];
    if (u == v) return max(val[u], val[v]);
    else if (head[u] == head[v]) return trees[head[u]].query(upos, vpos);
    else {
      int local;
      if (trees.count(head[u])) local = trees[head[u]].query(0, upos);
      else local = val[u];
      return max(local, query(parent[head[u]], v));
    }
  }

  hld(vvi& tree, vi& val, int root) : n(tree.size()), root(root), size(n),
				      parent(n, -1), level(n), head(n, -1), len(n), val(val) {
    pre(tree, root, 0);
    dfs(tree, root);
    make_st(tree);
  }
};

int main() {
  vvi tree(10, vi());
  tree[0] = {1, 6};
  tree[1] = {2, 5};
  tree[2] = {3, 4};
  tree[6] = {7};
  tree[7] = {8, 9};

  vi values = {1, 2, 3, 4, 5, 6, 15, 16, 9, 10};

  hld h(tree, values, 0);

  cout << "QUERY : " << h.query(3, 9) << endl;

  return 0;
}
