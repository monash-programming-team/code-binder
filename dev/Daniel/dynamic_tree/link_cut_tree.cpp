// Link-Cut Tree
// Based on code from https://github.com/saadtaame/link-cut-tree
// 
// Author: Daniel, based on code from https://github.com/saadtaame
// Reliability: 2
// Tested on: SPOJ-DYNALCA, SPOJ-DYNACON1
//
#include<bits/stdc++.h>

using namespace std;

//listings:lct
// Link-Cut Tree for dynamic connectivity on a forest of trees, dynamic lowest common
// ancestor queries and dynamic aggregate statistics for root-to-node paths. Default
// aggregate is node depths, can be customised. Complexity: O(log(N)) amortized queries
struct LinkCutTree {
  struct Node { int sz,i,f; Node *p,*pp,*l,*r; Node() : f(0),p(0),pp(0),l(0),r(0) {} };
  // Initialise: Create an initially disconnected forest of n isolated vertices --------
  LinkCutTree(int n) : V(n) { for(int i=0; i<n; i++) V[i].i = i, update(&V[i]); }
  // Update operations -----------------------------------------------------------------
  void link(int u, int v) { _link(&V[u], &V[v]); }              // Make u a subtree of v
  void cut(int u) { _cut(&V[u]); }                       // Disconnect u from its parent
  void make_root(int u) {                  // Make u the root of its connected component
    Node* x = &V[u]; access(x);
    if (x->l) x->l->p = 0, x->l->f ^= 1, x->l->pp = x, x->l = 0, update(x);
  }
  // Query operations ------------------------------------------------------------------
  int parent(int u) { access(&V[u]); return V[u].l ? V[u].l->i : -1; }    // Parent of u
  int root(int u) { return _root(&V[u])->i; }       // The root of the tree containing u
  bool connected(int u, int v) { return root(u) == root(v); }  // Are u and v connected?
  int lca(int u, int v) { return _lca(&V[u], &V[v])->i; }     // Find the LCA of u and v
  int query(int u) { return _query(&V[u]); }   // Aggregate path statistic query (depth)
   // OPTIONAL: Customise the aggregate path query below (default is node depth) ------
  void update(Node* x) { x->sz = 1 + (x->l ? x->l->sz : 0) + (x->r ? x->r->sz : 0); }
  int _query(Node* x) { access(x); return x->sz-1; }
  // Internal node operations (probably don't modify below here) -----------------------
  vector<Node> V;  
  Node* _root(Node* x) { access(x); while(x->l) { x=x->l;push(x); } splay(x); return x;}
  void _cut(Node* x) { access(x); x->l->p = 0; x->l = 0; update(x); }
  void _link(Node* x, Node* y) { access(x); access(y); x->l = y; y->p = x; update(x); }
  Node* _lca(Node* x, Node* y) { access(x); return access(y); }
  void push(Node* x) {  // Push lazy subtree flipping down the auxillary tree
    if (x->f == 0) return;  x->f = 0;  swap(x->l, x->r);
    if (x->l) x->l->f ^= 1; if (x->r) x->r->f ^= 1; update(x);
  }  // Splay tree right rotation for the auxillary trees
  void rotr(Node* x) {
    Node* y = x->p; Node* z = y->p;
    if((y->l = x->r)) y->l->p = y;
    x->r = y, y->p = x;
    if((x->p = z)) { if(y == z->l) z->l = x; else z->r = x; }
    x->pp = y->pp, y->pp = 0, update(y);
  }  // Splay tree left rotation for the auxillary trees
  void rotl(Node* x) {
    Node* y = x->p; Node* z = y->p;
    if((y->r = x->l)) y->r->p = y;
    x->l = y, y->p = x;
    if((x->p = z)) { if(y == z->l) z->l = x; else z->r = x; }
    x->pp = y->pp, y->pp = 0, update(y);
  }
  void splay(Node* x) {  // Rotates x to become the root of its auxillary tree
    for (Node* y = x->p; y; y = x->p) {
      if (x->p->p) push(x->p->p);  push(x->p); push(x);  // Push flips down the tree
      if(y->p == 0) { if (x == y->l) rotr(x); else rotl(x); }
      else { 
        if(y == y->p->l) { if(x == y->l) rotr(y), rotr(x); else rotl(x), rotr(x); }
        else { if(x == y->r) rotl(y), rotl(x); else rotr(x), rotl(x); }
      }
    }
    push(x), update(x);
  }    // Makes the root-to-v path preferred and makes v the root of its auxillary tree.
  Node* access(Node* x) {      // Returns the lowest ancestor of x in the root auxillary
    Node* last = x;  splay(x);        // tree (LCA with the most recently accessed node)
    if(x->r) x->r->pp = x, x->r->p = 0, x->r = 0, update(x);
    while(x->pp) {
      Node* y = x->pp; last = y; splay(y);
      if(y->r) y->r->pp = y, y->r->p = 0;
      y->r = x, x->p = y, x->pp = 0, update(y), splay(x);
    }
    return last;
  }
};
//listings:/lct

namespace problems {
  // Verdict: AC
  void solve_SPOJ_DYNALCA() {
    ios::sync_with_stdio(0); cin.tie(0);
    int N, M; cin >> N >> M;
    LinkCutTree tree(N);
    while (M--) {
      string q; cin >> q;
      if (q == "lca") {
        int A, B; cin >> A >> B; A--, B--;
        cout << tree.lca(A,B)+1 << '\n';
      } else if (q == "link") {
        int A, B; cin >> A >> B; A--, B--;
        tree.link(A,B);
      } else if (q == "cut") {
        int A; cin >> A; A--;
        tree.cut(A);
      }
    }
  }
  // Verdict: AC
  void solve_SPOJ_DYNACON1() {
    ios::sync_with_stdio(0); cin.tie(0);
    int N, M; cin >> N >> M;
    LinkCutTree tree(N);
    set<pair<int,int>> edges;
    int qq = 0;
    while (M--) {
      string q; int A, B;
      cin >> q >> A >> B; A--, B--;
      if (q[0] == 'c') puts(tree.connected(A,B) ? "YES" : "NO");
      else if (q[0] == 'a') {
        if (!tree.connected(A,B)) {
          assert(!edges.count(minmax(A,B)));
          tree.make_root(A), tree.link(A,B), edges.insert(minmax(A,B));
        }
      }
      else {
        if (edges.count(minmax(A,B))) {
          assert(tree.connected(A,B));
          tree.make_root(B), tree.cut(A), edges.erase(minmax(A,B));
        }
      }
    }
  }
}

namespace testing {
  void test() {
    LinkCutTree lct(10);
    for (int i=1; i<10; i++) lct.link(i, i-1);
    auto roots = [&]() {
      cout << "** ROOTS **" << endl;
      for (int u=0; u<10; u++) cout << "u=" << u << ", root=" << lct.root(u) << endl;
    };
    roots();
    lct.make_root(5);
    roots();
    lct.cut(6);
    roots();
  }
}

int main() {
  //problems::solve_SPOJ_DYNACON1();
  problem::solve_SPOJ_DYNALCA();
  //testing::test();
}
