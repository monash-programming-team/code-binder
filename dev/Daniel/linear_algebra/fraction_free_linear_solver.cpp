// Fraction-free linear system solver
//
// Author: Daniel (Based on Darcy's Code Binder)
// Date: 10-01-2017
// Reliability: 1
// Tested on: Randomly generated matrices
//
//
#include <bits/stdc++.h>
using namespace std;

#include "../code-template/debug.h"

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:fflinsolve
// solves Ax = b exactly. Returns {det, x_star}, solution is x_star[i] / det.
// T must be an integral data type (int, long long, etc.) Complexity: O(n^3)
template<typename T> pair<T,vector<T>> fflinsolve(vector<vector<T>> A, vector<T> b) {
  int k_c, k_r, pivot, sign = 1, n = (int)A.size(); T d = 1;
  for (k_c = k_r = 0; k_c < n; k_c++) {
    for (pivot = k_r; pivot < n && !A[pivot][k_r]; pivot++);
    if (pivot < n) {
      if (pivot != k_r) {
        for (int j = k_c; j < n; j++) swap(A[pivot][j], A[k_r][j]);
        swap(b[pivot], b[k_r]), sign *= -1;
      }
      for (int i = k_r + 1; i < n; i++) {
        for (int j = k_c + 1; j < n; j++)
          A[i][j] = (A[k_r][k_c] * A[i][j] - A[i][k_c] * A[k_r][j]) / d;
        b[i] = (A[k_r][k_c] * b[i] - A[i][k_c] * b[k_r]) / d,  A[i][k_c] = 0;
      }
      if (d) d = A[k_r][k_c];
      k_r++;
    } else d = 0;
  }
  if (!d) {
    for (int k = k_r; k < n; k++) if (b[k]) return {0,{}}; // inconsistent system
    return {0,{}};     // multiple solutions
  }
  vector<T> x_star(n);
  for (int k = n - 1; k >= 0; k--) {
    x_star[k] = sign * d * b[k];
    for (int j = k + 1; j < n; j++) x_star[k] -= A[k][j] * x_star[j];
    x_star[k] /= A[k][k];
  }
  return {sign * d, x_star};
}
//listings:/fflinsolve

int main() {
  vvi A = {{1,0},{0,1}};
  vi b = {1,1};
  int det; vi x_star;
  tie(det, x_star) = fflinsolve(A, b);
  DEBUG(det, x_star);
}
