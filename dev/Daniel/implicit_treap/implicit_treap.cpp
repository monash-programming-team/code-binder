// Implicit Treap (implicit Cartesian tree)
//
// Author: Daniel Anderson
// Date: 21-04-2017
// Reliability: 2
// Tested on: Moscow09-C, Brute-Force
//
#include<bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:treap
// Implicit Treap data structure. Supports ranged substring, erase, insert, reverse.
// Complexity: expected O(log(N)) per query, O(N) to build from or convert to vector
template<typename T> struct Node {
  typedef unique_ptr<Node> node_p; T val; node_p l, r; int size; bool rev;
  Node(T val) : val(val), rev(0) { update(); }
  node_p left() { auto t = move(l); update(); return t; }
  node_p right() { auto t = move(r); update(); return t; }
  void left(node_p t) { l = move(t); update(); }
  void right(node_p t) { r = move(t); update(); }
  void update() {
    size = 1 + (l ? l->size : 0) + (r ? r->size : 0);
    if (rev) { rev=0; swap(l, r); if (l) l->rev ^= 1; if (r) r->rev ^= 1; }
  }
};
template<typename T> struct ImplicitTreap {
  typedef Node<T> node; typedef ImplicitTreap<T> treap; 
  typedef unique_ptr<Node<T>> node_p;  node_p root;
  ImplicitTreap(const vector<T>& A) {  // Build an ImplicitTreap containing A
    function<node_p(int,int)> build = [&](int l, int r) {
      node_p v;  int m = (l+r)/2; if (l >= r) return v; v = make_unique<node>(A[m]);
      v->left(build(l, m)), v->right(build(m+1,r)); return v;
    }; root = build(0, A.size());
  }
  int size() { return root ? root->size : 0; }
  T& operator[](int i) { return lookup(root, i); } // Return the element at position i
  T& lookup(node_p& t, int key) {
    t->update(); int cur = (t->l ? t->l->size : 0); if (cur==key) return t->val;
    if (cur > key) return lookup(t->l, key); return lookup(t->r, key-cur-1);
  }
  ImplicitTreap(node_p root) : root(move(root)) { }  // Create a tree rooted at root
  treap cut(int l, int r) {  // Cut out and return the substring [l, r]
    node_p t1,t2,t3; tie(t1,t2)=split(move(root),l), tie(t2,t3)=split(move(t2),r-l+1);
    root = merge(move(t1), move(t3)); return treap(move(t2));
  }
  void insert(int i, treap&& other) {  // Insert the contents of 'other' at position i
    node_p t1, t2; tie(t1, t2) = split(move(root), i);
    root = merge(move(t1), move(other.root)); root = merge(move(root), move(t2));
  }
  void insert(int i, treap& other) { insert(i, move(other)); }
  void reverse(int l, int r) {  // Reverse the contents of [l, r]
    node_p t1,t2,t3; tie(t1,t2)=split(move(root),l); tie(t2,t3)=split(move(t2),r-l+1);
    t2->rev ^= 1; root = merge(move(t1), move(t2)), root = merge(move(root),move(t3));
  }
  pair<node_p, node_p> split(node_p t, int key, int add=0) {
    if (!t) return {nullptr, nullptr};
    t->update(); int cur = add + (t->l ? t->l->size : 0);
    if (key <= cur) {  // Recursively split the left subtree
      node_p left,right,tmp = t->left(); tie(left,right)=split(move(tmp),key,add);
      return {move(left), merge(move(right), move(t))};
    } else {  // Recursively split the right subtree
      node_p left,right,tmp = t->right(); tie(left,right) = split(move(tmp),key,cur+1);
      return {merge(move(t), move(left)), move(right)};
    }
  }
  node_p merge(node_p l, node_p r) {  // Merge the trees rooted at l and r
    if (l) l->update(); if (r) r->update();  if (!l || !r) return l ? move(l) : move(r);
    bool left = (1.0*rand()/RAND_MAX) < (1.0 * l->size) / (l->size + r->size);
    if (left) { l->right(merge(l->right(), move(r))); return l; }   // Merge randomly to
    else { r->left(merge(move(l), r->left())); return r; }  // maintain expected balance
  }
  vector<T> to_vector() {  // Convert the contents of the tree into a vector<T>
    vector<T> res; res.reserve(size()); function<void(node_p&)> go = [&](node_p& v) {
      if (!v) return; v->update(); go(v->l); res.push_back(v->val); go(v->r);
    }; go(root); return res;
  }
};
//listings:/treap

void solve_CODING() {
  #ifdef ONLINE_JUDGE
    freopen("coding.in", "r", stdin);
    freopen("coding.out", "w", stdout);
  #endif
	ios::sync_with_stdio(0); cin.tie(0);
  
  string s;	cin >> s;
	ImplicitTreap<char> txt({s.begin(),s.end()});
	
	int M; cin >> M;
	
	vvi queries;
	for (int m = 0; m < M; m++) {
		int i, j, k;
		cin >> i >> j >> k;
		i--; j--;
		queries.push_back({i,j,k});
	}
	reverse(queries.begin(), queries.end());
	
	for (auto& q : queries) {
		int i = q[0], j = q[1], k = q[2];
    txt.insert(j-k+1, txt.cut(i,i+k-1));
	}
  auto vec = txt.to_vector();
  cout << string(vec.begin(), vec.end()) << endl;
}

typedef ImplicitTreap<int> rope;

void test_reverse() {
  // Randomised testing
  srand(time(0));
  vi A(1000); for (auto& x : A) x = rand();
  rope t(A);
  
  int testno=1;
  while(1) {
    if (testno++ % 1000 == 0) cout << "Test " << testno << "                  \r" << flush;
    int l = rand() % 1000, r = rand() % 1000;
    if (r < l) swap(l, r);
    reverse(A.begin()+l, A.begin()+r+1);
    t.reverse(l, r);
    if (rand() % 10 == 0) assert(t.to_vector() == A);
  }
}

void test_cut() {
  srand(time(0));
  vi A(1000); for (auto& x : A) x = rand();
  rope t(A);
  
  int testno=1;
  while(1) {
    if (testno++ % 1000 == 0) cout << "Test " << testno << "                  \r" << flush;
    int l = rand() % 1000, r = rand() % 1000;
    if (r < l) swap(l, r);
    vi B(A.begin()+l, A.begin()+r+1);
    A.erase(A.begin()+l, A.begin()+r+1);
    int pos = rand() % (A.size()+1);
    A.insert(A.begin()+pos, B.begin(), B.end());
    auto sub = t.cut(l, r);
    t.insert(pos, sub);
    if (rand() % 10 == 0) assert(t.to_vector() == A);
  }
}

void test_all() {
  srand(time(0));
  vi A(1000); for (auto& x : A) x = rand();
  rope t(A);
  
  int testno=1;
  while(1) {
    if (testno++ % 1000 == 0) cout << "Test " << testno << "                  \r" << flush;
    if (rand() % 2 == 0) {  // Test cut + insert
      int l = rand() % 1000, r = rand() % 1000;
      if (r < l) swap(l, r);
      vi B(A.begin()+l, A.begin()+r+1);
      A.erase(A.begin()+l, A.begin()+r+1);
      int pos = rand() % (A.size()+1);
      A.insert(A.begin()+pos, B.begin(), B.end());
      auto sub = t.cut(l, r);
      t.insert(pos, sub);
    } else {  // test reverse
      int l = rand() % 1000, r = rand() % 1000;
      if (r < l) swap(l, r);
      reverse(A.begin()+l, A.begin()+r+1);
      t.reverse(l, r);
    }
    if (rand() % 10 == 0) {
      assert(t.to_vector() == A);
      for (int i=0; i<t.size(); i++) assert(A[i] == t[i]);
    }
  }
}


int main() {
  solve_CODING();
  //test_reverse();
  //test_cut();
  //test_all();
}

