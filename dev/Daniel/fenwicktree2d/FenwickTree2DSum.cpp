// 2D Fenwick Tree
#include<bits/stdc++.h>
using namespace std;

#include "../code-template/debug.h"

namespace PointAdditionRangeSum {

//listings:point_add_range_query
// 2D Fenwick Tree supporting point addition and ranged sum. Complexity: O(log(n)log(m))
template<typename T> struct FenwickTree2D {
  int N, M; vector<vector<T>> A;
  FenwickTree2D(int n, int m) : N(n+1), M(m+1), A(N, vector<T>(M)) { }
  void add(int r, int c, T v) {  // A[r][c] += v
    for (int i=r+1;i<N;i+=i&-i) for (int j=c+1;j<M;j+=j&-j) A[i][j]+=v; }
  T sum(int r, int c) {  // sum(A[0..r)[0..c))
    T s=0; for (int i=r;i;i-=i&-i) for (int j=c;j;j-=j&-j) s+=A[i][j]; return s; }  
  T sum(int r1, int c1, int r2, int c2) {  // sum(A[r1..r2)[c1..c2))
    return sum(r2,c2)-sum(r1,c2+1)-sum(r2+1,c1)+sum(r1,c1); }
};
//listings:/point_add_range_query

}

namespace RangeAdditionPointQuery {

//listings:range_add_point_query
// 2D Fenwick Tree supporting range add and point queries. Complexity: O(log(n)log(m))
template<typename T> struct FenwickTree2D {
  int N, M; vector<vector<T>> A;
  FenwickTree2D(int n, int m) : N(n+1), M(m+1), A(N, vector<T>(M)) { }
  void add(int r, int c, int v) {  // A[0..r)[0..c) += v
    for (int i=r;i;i-=i&-i) for (int j=c;j;j-=j&-j) A[i][j] += v; }
  void add(int r1, int c1, int r2, int c2, T v) {  // A[r1..r2)[c1..c2) += v
    add(r2,c2,v); add(r1,c2+1,-v); add(r2+1,c1,-v); add(r1,c1,v); }
  T query(int r, int c) {  // get A[r][c]
    T s=0; for (int i=r+1;i<N;i+=i&-i) for (int j=c+1;j<M;j+=j&-j) s+=A[i][j]; return s; }
};
//listings:/range_add_point_query

}

int main() {
  PointAdditionRangeSum::FenwickTree2D<int> ft(5, 5);
  ft.add(1, 1, 1);
  ft.add(2,2,2);
  cout << ft.sum(4,4) << endl;

  RangeAdditionPointQuery::FenwickTree2D<int> ft2(5,5);
  ft2.add(0,0,5,5,10);
  cout << ft2.query(3,3) << endl;
  ft2.add(2,2,3,3,5);
  cout << ft2.query(2,2) << endl;

}