// GNU Extra Builtin Data Structures
//
// Author: Daniel
// Date: 10-01-2017
// Reliability: 0
// Tested on:
//
#include<bits/stdc++.h>
using namespace std;

//listings:prefix_trie
// GNU Policy-Based Data Structures --------------------------------------------------
// prefix_trie:: A Patricia (compact) trie that implements fast prefix searches.
// Insertion syntax matches std::set::insert, returns pair{iterator, success}
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/trie_policy.hpp>
#include <ext/pb_ds/tag_and_trait.hpp>
using namespace __gnu_pbds;

typedef trie<string, null_type, trie_string_access_traits<>, pat_trie_tag,
	     trie_prefix_search_node_update> prefix_trie;
//listings:/prefix_trie

void test_prefix_trie() {
//listings:prefix_trie_example
// Usage example for Patricia trie
prefix_trie t;                              // Create an empty prefix trie
t.insert("Banana");                         // Insert an element
auto match_range = t.prefix_range("Ban");   // Get all strings matching "Ban*"

for (auto it = match_range.first; it != match_range.second; ++it) cout << *it << ' ';
//listings:/prefix_trie_example
cout << endl;
}

//listings:rope
// GNU Policy-Based Data Structures --------------------------------------------------
// rope:: An Implicit Cartesian Tree; a data structure that allows for 
// fast [O(log(n)] insertion and deletion of arbitrarily long blocks of data.
// Uses most of the same syntax as vector. See examples.
#include <ext/rope>
using namespace __gnu_cxx;
//listings:/rope

void test_rope() {
int n = 10;
int pos = 5;
int length = 5;
//listings:rope_example
// Usage example for rope
rope<int> v;                                // create an empty rope.
for (int i=0; i<n; ++i) v.push_back(i);     // insert into rope
rope<int> cur = v.substr(pos, length);      // get substring from [pos, pos+length)
v.erase(pos, length);                       // erase substring from [pos, pos+length)
v.insert(v.mutable_begin(), cur);           // use mutable_begin for non-const iterator
for (const auto& x : v) cout << x << ' ';   // iterate over the contents of the rope
//listings:/rope_example
cout << endl;
}


/*
//listings:ordered_set
// GNU Policy-Based Data Structures --------------------------------------------------
// ordered_set:: A red-black tree maintaining node order-statistics, allowing for
// fast [O(log(n)] order-statistics queries. Uses the same syntax as std::set.
#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;

typedef tree<int, null_type, less<int>, rb_tree_tag,
		tree_order_statistics_node_update> ordered_set;
//listings:/ordered_set
*/

#include <ext/pb_ds/tree_policy.hpp>

typedef tree<int, null_type, less<int>, rb_tree_tag,
		tree_order_statistics_node_update> ordered_set;

void test_ordered_set() {
int n = 10;
//listings:ordered_set_example
// Usage example for ordered_set
ordered_set s;                          // Create an empty ordered set
for (int i=0; i<n; i++) s.insert(i);    // Insert into the set
cout << *s.find_by_order(3) << endl;    // Find the 3rd element
cout << s.order_of_key(5) << endl;      // Find the order-statistic of 5
//listings:/ordered_set_example
}

int main() {
  test_prefix_trie();
  test_rope();
  test_ordered_set();
}