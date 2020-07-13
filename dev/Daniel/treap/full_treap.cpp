#include<bits/stdc++.h>
using namespace std;

typedef long long ll;

const ll INF = LLONG_MAX/2;

std::mt19937 rng;
std::uniform_int_distribution<std::mt19937::result_type> dist(0,INT_MAX); 

// Customisable Treap data structure using randomised priorities. This example
// can compute for any range [i,j], the distance between closest pair of elements.
// Complexity: expected O(log(N)) per query.
template<typename K, typename V> struct Node {
  typedef unique_ptr<Node<K,V>> node_p; K key; V val; int p, size=1;  node_p l=0, r=0;
  Node(K key, V val) : key(key), val(val), p(dist(rng)) { update(); }
  node_p left() { auto t = move(l); update(); return t; }
  node_p right() { auto t = move(r); update(); return t; }
  void set_left(node_p t) { l = move(t); update(); }
  void set_right(node_p t) { r = move(t); update(); }
  K minv = INF, maxv = -INF, mind = INF;
  void update() {
    size = 1 + (l ? l->size : 0) + (r ? r->size : 0);
    minv = key, maxv = key, mind = INF;
    if (l) minv = l->minv, mind = min(min(mind, abs(key - l->maxv)), l->mind);
    if (r) maxv = r->maxv, mind = min(min(mind, abs(r->minv - key)), r->mind);
  }
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
    if (!a) return move(b); if (!b) return move(a);
    if (a->p < b->p) { a->set_right(merge(a->right(), move(b))); return a; }
    else { b->set_left(merge(move(a), b->left())); return b; }
  }
  node_p insert(node_p t, K key, V val) {  // Insert(key,val) into the given subtree
    if (!t) return make_unique<node>(key, val);
    if (key < t->key) t->set_left(insert(t->left(), key, val));
    else if (key > t->key) t->set_right(insert(t->right(), key, val));
    else t->val = val;
    return normalise(move(t));
  }
  node_p remove(node_p t, K key) {  // Remove key from the given subtree
    if (!t) return t;
    if (key < t->key) { t->set_left(remove(t->left(), key)); return t; }
    if (key > t->key) { t->set_right(remove(t->right(), key)); return t; }
    return merge(t->left(), t->right());
  }
  node_p normalise(node_p t) {  // Ensure that the heap-ordering of the p's is correct
    if (t->l && t->l->p < t->p && (!t->r || t->l->p < t->r->p)) {
      auto tmp = t->left(); t->set_left(tmp->right());
      tmp->set_right(move(t)); return move(tmp);
    } else if (t->r && t->r->p < t->p) {
      auto tmp = t->right(); t->set_right(tmp->left());
      tmp->set_left(move(t)); return move(tmp);
    } else return t;
  }
  K find_by_order(node* n, int k) {  // Find the k'th order statistic of the set
    int ord = (n->l ? n->l->size : 0);
    if (ord == k) return n->val;
    else if (ord > k) return find_by_order(n->l.get(), k);
    else return find_by_order(n->r.get(), k-ord-1);
  }
  K min_diff(node* n, int i, int j) {  // Find the minimum difference in elements [i,j]
    K ans = INF, ord = (n->l ? n->l->size : 0);
    if (i <= 0 && j >= n->size-1) return n->mind;
    if (n->l && i < ord) {
      if (j >= ord) ans = min(ans, n->val - n->l->maxv);
      ans = min(ans, min_diff(n->l.get(), i, j));
    }
    if (n->r && j > ord) {
      if (i <= ord) ans = min(ans, n->r->minv - n->val);
      ans = min(ans, min_diff(n->r.get(), i-ord-1, j-ord-1));
    }
    return ans;
  }
  int order_of_key(node* n, int key) {  // Count the elements < the given key
    if (key > n->maxv) return n->size;
    if (key <= n->minv) return 0;
    int ord = (n->l ? n->l->size : 0);
    if (n->key == key) return ord;
    else if (key < n->key) return order_of_key(n->l.get(), key);
    else return ord+1 + order_of_key(n->r.get(), key);
  }
  K find_by_order(int k) { return find_by_order(root.get(), k); }
  K min_diff(int i, int j) { return min_diff(root.get(), i, j); }
  int order_of_key(int key) { return size() ? order_of_key(root.get(), key) : 0; }
  int size() { return root ? root->size : 0; }
};


void solve_SPOJ_TREAP() {
  rng.seed(std::random_device()());
  Treap<ll,ll> t;
  int Q; cin >> Q;
  while (Q--) {
    char type; cin >> type;
    if (type == 'I') {
      int k; cin >> k; t.insert(k, k);
    } else if (type == 'D') {
      int k; cin >> k; t.remove(k);
    } else if (type == 'N') {
      int i, j; cin >> i >> j;
      if (i == j) cout << -1 << '\n';
      else cout << t.min_diff(i, j) << '\n';
    } else {
      int i, j; cin >> i >> j;
      if (i == j) cout << -1 << '\n';
      else cout << t.find_by_order(j) - t.find_by_order(i) << '\n';
    }
  }
}

void solve_SPOJ_ORDERSET() {
  rng.seed(std::random_device()());
  Treap<ll,ll> t;
  int Q; cin >> Q;
  while (Q--) {
    char type; cin >> type;
    if (type == 'I') {
      int k; cin >> k; t.insert(k, k);
    } else if (type == 'D') {
      int k; cin >> k; t.remove(k);
    } else if (type == 'K') {
      int k; cin >> k; k--;
      if (k >= t.size()) cout << "invalid\n";
      else cout << t.find_by_order(k) << '\n';
    } else {
      int x; cin >> x;
      cout << t.order_of_key(x) << '\n';
    }
  }
}

int main() {
  solve_SPOJ_ORDERSET();
}