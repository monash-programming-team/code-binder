// Linear Diophantine Solver
//
// Author: Daniel (Based on Darcy's Code Binder)
// Date: 10-01-2017
// Reliability: 1
// Tested on: ECNA10:H
//
//
#include <bits/stdc++.h>
using namespace std;

#include "../code-template/debug.h"

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:diophantine_linsolve
// Integral matrix triangulation. Used by linear diophantine solver below.
template<typename T> int triangulate(vector<vector<T>>& A, int m, int n, int cols) {
  lldiv_t d;  int ri = 0, ci = 0;
  while (ri < m && ci < cols) {
    int pi = -1;
    for (int i = ri; i < m; i++)
      if (A[i][ci] && (pi == -1 || abs(A[i][ci]) < abs(A[pi][ci]))) pi = i;
    if (pi == -1) ci++;
    else {
      int k = 0;
      for (int i = ri; i < m; i++) if (i != pi) {
        d = lldiv(A[i][ci], A[pi][ci]);
        if (d.quot) { for (int j = ci; j < n; j++) A[i][j] -= d.quot * A[pi][j]; k++; }
      }
      if (!k) { for (int i=ci; i<n && ri!=pi; i++) swap(A[ri][i],A[pi][i]); ri++,ci++; }
    }
  }
  return ri;
}

// System of linear diophantine equations A*x = b. T must be an integral type.
// Returns dim(null space), or -1 if there is no solution, or -2 if inconsistent.
// xp: a particular solution
// basis: an n x n matrix whose first dim columns form a basis of the nullspace.
// All solutions are obtained by adding integer multiples the basis elements to xp.
// Complexity: O(n^3)
template<typename T> tuple<int,vector<T>,vector<vector<T>>> diophantine_linsolve(
    vector<vector<T>>& A, vector<T>& b) {                     
  int m = (int)A.size(), n = (int)A[0].size(), i, j, rank;  T d;
  vector<vector<T>> mat(n + 1, vector<T>(m + n + 1));                       
  for (i = 0; i < m; i++) mat[0][i] = -b[i];
  for (i = 0; i < m; i++) for (j = 0; j < n; j++) mat[j + 1][i] = A[i][j];
  for (i = 0; i < n + 1; i++) for (j = 0; j < n + 1; j++) mat[i][j + m] = (i == j);
  rank = triangulate(mat, n + 1, m + n + 1, m + 1), d = mat[rank - 1][m];
  vector<vector<T>> basis(n, vector<T>(n));  vector<T> xp(n);
  if (d != 1 && d != -1) return make_tuple(-1, xp, basis);
  for (i = 0; i < m; i++) if (mat[rank - 1][i]) return make_tuple(-2, xp, basis);
  for (i = 0; i < n; i++) {
    xp[i] = d * mat[rank - 1][m + 1 + i];
    for (j = 0; j < n + 1 - rank; j++) basis[i][j] = mat[rank + j][m + 1 + i];
  }
  return make_tuple(n + 1 - rank, xp, basis);
}
//listings:/diophantine_linsolve

int main() {
  vvi A = {{1,0},{0,1}};
  vi b = {1,1};
  
  int dim; vi xp; vvi basis;
  tie(dim, xp, basis) = diophantine_linsolve(A, b);
  
  DEBUG(dim, xp, basis);
}
