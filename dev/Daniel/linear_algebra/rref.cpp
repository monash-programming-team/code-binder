// Reduction to reduced row echelon form
//
// Author      : Daniel Anderson (Based on Stanford's ACM ICPC Code Binder)
// Date        : 06/01/2017
// Reliability : 1
// Tested on   : LightSwitches
//
#include <bits/stdc++.h>
using namespace std;

#include "../code-template/debug.h"

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<double> vd;
typedef vector<vd> vvd;

//listings:rref
// Reduces the given matrix to reduced row-echelon form using Gaussian Elimination.
// Returns the rank of A. T must be a floating-point type. Complexity: O(n^3).
const double EPS = 1e-10;

template<typename T> int rref(vector<vector<T>>& A) {
  int n = (int)A.size(), m = (int)A[0].size(), r = 0;
  for (int c=0; c<m && r<n; c++) {
    int j = r;
    for (int i=r+1; i<n; i++) if (abs(A[i][c]) > abs(A[j][c])) j = i;
    if (abs(A[j][c]) < EPS) continue;
    swap(A[j], A[r]);  T s = 1.0 / A[r][c];
    for (int j=0; j<m; j++) A[r][j] *= s;
    for (int i=0; i<n; i++) if (i != r) {
      T t = A[i][c];
      for (int j=0; j<m; j++) A[i][j] -= t * A[r][j];
    }
    r++;
  }
  return r;
}

//listings:/rref



int main() {
  vvd A = {{2,1,0}, {3,0,1}, {4,10,-0.5}};
  rref(A);
  DEBUG(A);
  
  A = {{1,1,1}, {1,1,1}, {1,1,1}};
  rref(A);
  DEBUG(A);
}
