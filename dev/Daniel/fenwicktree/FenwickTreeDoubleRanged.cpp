// Fenwick Tree supporting ranged sum and ranged addition
//
// Uses two Fenwick Trees internally:
//  - For each addition of VAL on [0, b]:
//    t1[b] stores -VAL, t2[b] stores VAL * b
//    
//  - To query [0, b]:
//    Return sum(t1[0,b]) * b + sum(t2[0,b])
//
//
#include<bits/stdc++.h>
using namespace std;

template<typename T> struct FenwickTree {
  int n; vector<T> t1, t2;
  FenwickTree(int n) : n(n), t1(n), t2(n) { }
  FenwickTree(const vector<T>& A) : n(A.size()), t1(n), t2(n) {
    for (int i=0; i<n; i++) { t2[i] += A[i]; if (i|(i+1) < n) t2[i|(i+1)] += t2[i]; }
  }
  void add(int a, int b, T val) { add(t1,a,val), add(t1,b,-val), add(t2,a,-val*(a-1)), add(t2,b,val*b); }
  T query(int a, int b) { return query(b) - (a ? query(a - 1) : 0); }
  T query(int b) { return sum(t1, b) * b + sum(t2, b); }
  void add(vector<T>& t, int i, int val) { for (; i < n; i |= i+1) t[i] += val; }
  T sum(vector<T>& t, int i) { T res = 0; for (; i >= 0; i = (i&(i+1))-1) res += t[i]; return res; }
};

int main() {
  FenwickTree<int> ft(10);
  ft.add(3,8,10);
  cout << ft.query(5,9) << endl;
}