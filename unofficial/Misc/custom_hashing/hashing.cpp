// Hashing custom types example
//
// Author: Daniel. Based on example in Darcy's code binder
// Date: 19-01-2017
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;

//listings:hashing
// Custom type hashing example. Your type must implement equality. The hash
// function must be consistent with ==, that is (a==b) => (hash(a)==hash(b))
struct MyType {
  int a; string b;
  bool operator==(const MyType& r) const { return a == r.a && b == r.b; }
};
namespace std {
  template<> struct hash<MyType> {
    size_t operator()(const MyType& x) const {
      return hash<int>()(x.a) ^ hash<string>()(x.b);
    }
  };
}

unordered_map<MyType,string> my_map;
//listings:/hashing

int main() {
  my_map[{10,"hello"}] = "my first entry";
  my_map[{20,"bye"}] = "my second entry";
  
  cout << my_map[{10,"hello"}] << endl;
  cout << my_map[{20,"bye"}] << endl;
  cout << my_map[{30,"uh oh"}] << endl;
}