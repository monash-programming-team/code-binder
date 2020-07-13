// Linear system solver using LU decomposition
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

//listings:lu_decomp
// LU-Decomposition. Can be used to solve Ax = b in floating-point
// Returns {determinant, pivot, LU}. Complexity: O(n^3)
// - Call LU_solve(LU, pivot, b) to solve linear system Ax = b
const double EPS = 1e-9;

template<typename T> tuple<T, vi, vector<vector<T>>> LU_decomp(vector<vector<T>> A) {
  int n = (int)A.size(); vi pivot(n); vector<T> s(n); T c, t, det = 1.0;
  for (int i = 0; i < n; i++) {
    s[i] = 0.0;
    for (int j = 0; j < n; j++) s[i] = max(s[i], fabs(A[i][j]));
    if (s[i] < EPS) return make_tuple(0, pivot, A);  // Singular
  }
  for (int k = 0; k < n; k++) {
    c = fabs(A[k][k] / s[k]), pivot[k] = k;
    for (int i = k + 1; i < n; i++) if ((t = fabs(A[i][k] / s[i])) > c)
      c = t, pivot[k] = i;
    if (c < EPS) return make_tuple(0, pivot, A);  // Singular
    if (k != pivot[k]) {
      det *= -1.0;  swap(s[k], s[pivot[k]]);
      swap_ranges(A[k].begin() + k, A[k].end(), A[pivot[k]].begin() + k);
    }
    for (int i = k + 1; i < n; i++) {
      A[i][k] /= A[k][k];
      for (int j = k + 1; j < n; j++) A[i][j] -= A[i][k] * A[k][j];
    }
    det *= A[k][k];
  }
  return make_tuple(det, pivot, A);
}

// Solve Ax = b in floating-point using the LU-decomposition of A.
// T must be a floating-point type (double, long double). Complexity: O(n^2)
template<typename T> vector<T> LU_solve(vector<vector<T>>& LU, vi& piv, vector<T>& b) {
  int n = (int)LU.size();  vector<T> x = b;
  for (int k = 0; k < n - 1; k++) {
    if (k != piv[k]) swap(x[k], x[piv[k]]);
    for (int i = k + 1; i < n; i++) x[i] -= LU[i][k] * x[k];
  }
  for (int i = n - 1; i >= 0; i--) {
    for (int j = i + 1; j < n; j++) x[i] -= LU[i][j] * x[j];
    x[i] /= LU[i][i];
  }
  return x;
}
//listings:/lu_decomp


int main() {
  vector<vector<double>> A = {{1.0, 0}, {0, 1.0}};
  vector<double> b = {1.0, 1.0};
  
  double det; vi pivot; vector<vector<double>> LU;
  tie(det, pivot, LU) = LU_decomp(A);
  
  vector<double> sol = LU_solve(LU, pivot, b);
  DEBUG(sol);
}

