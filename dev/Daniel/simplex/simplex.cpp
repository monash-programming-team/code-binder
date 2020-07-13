// Two-phase simplex algorithm for solving linear programs.
//
// Author: Stanford ICPC Team
// Date: 10-01-2017
// Reliability: 1
// Tested on: CF605C
//
#include <bits/stdc++.h>
using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:simplex
// Simplex algorithm for solving linear programs of the form
//     maximize     c^T x
//     subject to   Ax <= b
//                  x >= 0
// solve() returns INF if unbounded, NaN if infeasible.
// T must be a floating-point type. Complexity is unbounded in general.
const double EPS = 1e-9;

template<typename T> struct LPSolver {
  const T INF = numeric_limits<T>::infinity(), NaN = numeric_limits<T>::quiet_NaN();
  int m, n;  vi N, B;  vector<vector<T>> D;
  void pivot(int r, int s) {
    T inv = 1.0 / D[r][s];
    for (int i = 0; i < m + 2; i++) if (i != r)
      for (int j = 0; j < n + 2; j++) if (j != s)
        D[i][j] -= D[r][j] * D[i][s] * inv;
    for (int j = 0; j < n + 2; j++) if (j != s) D[r][j] *= inv;
    for (int i = 0; i < m + 2; i++) if (i != r) D[i][s] *= -inv;
    D[r][s] = inv;  swap(B[r], N[s]);
  }
  bool simplex(int phase) {
    int x = phase == 1 ? m + 1 : m, s = -1, r = -1;
    for (; ; s=-1,r=-1) {
      for (int j = 0; j <= n; j++) if (!(phase == 2 && N[j] == -1))
        if (s == -1 || D[x][j] < D[x][s] || (D[x][j] == D[x][s] && N[j] < N[s])) s = j;
      if (D[x][s] > -EPS) return true;
      for (int i = 0; i < m; i++) if (!(D[i][s] < EPS))
        if (r == -1 || D[i][n + 1] / D[i][s] < D[r][n + 1] / D[r][s] ||
          ((D[i][n + 1] / D[i][s]) == (D[r][n + 1] / D[r][s]) && B[i] < B[r])) r = i;
      if (r == -1) return false;
      pivot(r, s);
    }
  }
  // Create a solver for max(c^T x) st. Ax <= b, x >= 0.
  LPSolver(const vector<vector<T>>& A, const vector<T>& b, const vector<T>& c) :
    m(b.size()), n(c.size()), N(n + 1), B(m), D(m + 2, vector<T>(n + 2)) {
    for (int i = 0; i < m; i++) for (int j = 0; j < n; j++) D[i][j] = A[i][j];
    for (int i = 0; i < m; i++) { B[i] = n + i; D[i][n] = -1; D[i][n + 1] = b[i]; }
    for (int j = 0; j < n; j++) { N[j] = j; D[m][j] = -c[j]; }
    N[n] = -1; D[m + 1][n] = 1;
  }
  pair<T, vector<T>> solve() {  // Returns {objective_value, optimial_solution}
    int r = 0; vector<T> x = vector<T>(n);
    for (int i = 1; i < m; i++) if (D[i][n + 1] < D[r][n + 1]) r = i;
    if (D[r][n + 1] < -EPS) {
      pivot(r, n);
      if (!simplex(1) || D[m + 1][n + 1] < -EPS) return {NaN, x};
      for (int i = 0; i < m; i++) if (B[i] == -1) {
        int s = -1;
        for (int j = 0; j <= n; j++) if (s == -1 || D[i][j] < D[i][s]
          || (D[i][j] == D[i][s] && N[j] < N[s])) s = j;
        pivot(i, s);
      }
    }
    if (!simplex(2)) return {INF, x};
    for (int i = 0; i < m; i++) if (B[i] < n) x[B[i]] = D[i][n + 1];
    return {D[m][n + 1], x};
  }
};
//listings:/simplex

namespace problems {
  // Verdict: AC
  namespace CF605C {
    void solve() {
      std::cout.sync_with_stdio(0); cin.tie(0);
      int n, p, q;
      
      // objective function
      cin >> n >> p >> q;
      vector<double> c(n, -1);  
     
      // constraints
      vector<vector<double>> A(2, vector<double>(n));
      vector<double> b = {-1.0*p, -1.0*q};
      for (int i=0; i<n; i++) {
        int ai, bi; cin >> ai >> bi;
        A[0][i] = -ai; A[1][i] = -bi;
      }
      
      // run the simplex
      LPSolver<double> lp(A, b, c);
      double res = -lp.solve().first;
      cout << fixed << setprecision(15) << res << endl;
    }
  }
}

int main() {
  problems::CF605C::solve();
}
