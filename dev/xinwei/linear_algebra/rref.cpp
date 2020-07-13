// Reduction to reduced row echelon form
//
// Author      : Daniel Anderson (Based on Stanford's ACM ICPC Code Binder)
// Date        : 06/01/2017
// Reliability : 0
// Tested on   :
//
#include <bits/stdc++.h>
using namespace std;

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

class LightSwitches {
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
                    for (int j=0; j<m; j++) {
                        A[i][j] -= t * A[r][j];
                        if (abs(A[i][j] + 1) < EPS) A[i][j] = 1.0;
                    }
                }
            r++;
        }
        return r;
    }
public:
    ll countPossibleConfigurations(vector<string> switches){
        int n = (int)switches.size();
        int m = (int)switches[0].size();
        vector<vector<double>> A(n, vector<double>(m, 0.0));
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                A[i][j] = (switches[i][j] == 'Y');
        int rank = rref(A);
        return 1LL << rank;
    }
};

int main() {


}
