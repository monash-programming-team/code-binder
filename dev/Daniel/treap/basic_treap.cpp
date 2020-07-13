// Customisable Treap
//
// Author: Daniel Anderson
// Date: 21-04-2017
// Reliability: 2
// Tested on: SPOJ-TREAP, SPOJ-ORDERSET
//
#include<bits/stdc++.h>
using namespace std;

//listings:treap
// Customisable Treap data structure. Complexity: expected O(log(N)) per query.
template<typename K, typename V> struct Node {
  typedef unique_ptr<Node<K,V>> node_p; 
  K key; V val;  int p, size=1;  node_p l=0, r=0;
  Node(K key, V val) : key(key), val(val), p(rand()) { update(); }
  node_p left() { auto t = move(l); update(); return t; }
  node_p right() { auto t = move(r); update(); return t; }
  void left(node_p t) { l = move(t); update(); }
  void right(node_p t) { r = move(t); update(); }
  void update() { size = 1 + (l ? l->size : 0) + (r ? r->size : 0); }
};
template<typename K, typename V> struct Treap {
  typedef Node<K,V> node; typedef Treap<K,V> treap; typedef unique_ptr<node> node_p;
  node_p root;  Treap() { }  // Construct an empty Treap
  // Constructs a Treap by merging the Treaps t1 and t2 where t1 < t2
  Treap(treap& t1, treap& t2) : root(merge(move(t1.root), move(t2.root))) { }
  void insert(K key, V val) {  // Insert the (key,value) pair into the Treap
    if (!root) root = make_unique<node>(key, val);
    else root = insert(move(root), key, val);
  } // Remove the item with the given key from the Treap if it exists
  void remove(K key) { if (root) root = remove(move(root), key); }
  // Split the Treap into two Treaps, containing all keys < key and >= key respectively
  pair<treap, treap> split(K key) {       
    node_p left, right; tie(left, right) = split(move(root), key);
    return {treap(move(left)), treap(move(right))};
  }  // Create a Treap owning the given root
  Treap(node_p root) : root(move(root)) { }  
  pair<node_p, node_p> split(node_p t, K key) {  // Split the subtree t at key
    if (!t) return {nullptr, nullptr};
    if (t->key < key) {  // Change < to <= if you want a {<=, >} split
      node_p left,right,tmp = t->right(); tie(left, right) = split(move(tmp), key);
      return {merge(move(t), move(left)), move(right)};
    } else {
      node_p left,right,tmp = t->left(); tie(left, right)=split(move(tmp), key);
      return {move(left), merge(move(right), move(t))};
    }
  }
  node_p merge(node_p a, node_p b) {  // Merge the subtrees a and b where a < b
    if (!a) return b; if (!b) return a;
    if (a->p < b->p) { a->right(merge(a->right(), move(b))); return a; }
    else { b->left(merge(move(a), b->left())); return b; }
  }
  node_p insert(node_p t, K key, V val) {  // Insert(key,val) into the given subtree
    if (!t) return make_unique<node>(key, val);
    if (key < t->key) t->left(insert(t->left(), key, val));
    else if (key > t->key) t->right(insert(t->right(), key, val));
    else t->val = val;
    return normalise(move(t));
  }
  node_p remove(node_p t, K key) {  // Remove key from the given subtree
    if (!t) return t;
    if (key < t->key) { t->left(remove(t->left(), key)); return t; }
    if (key > t->key) { t->right(remove(t->right(), key)); return t; }
    return merge(t->left(), t->right());
  }
  node_p normalise(node_p t) {  // Ensure that the heap-ordering of the p's is correct
    if (t->l && t->l->p < t->p && (!t->r || t->l->p < t->r->p)) {
      auto tmp = t->left(); t->left(tmp->right());
      tmp->right(move(t)); return tmp;
    } else if (t->r && t->r->p < t->p) {
      auto tmp = t->right(); t->right(tmp->left());
      tmp->left(move(t)); return tmp;
    } else return t;
  }
};
//listings:/treap

typedef Treap<int,int> treap;

int main() {
  treap t;
  t.insert(5,5);
  t.insert(3,3);
  t.insert(1,1);
  t.insert(6,1);
  t.insert(8,10);
  t.remove(5);
  t.insert(11,1);
  t.insert(50,2);
  t.insert(-10,1);
  
  treap l,r; tie(l,r) = t.split(6);
  t = treap(l, r);

}