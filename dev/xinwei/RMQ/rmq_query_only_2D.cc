// 2D version - Range Minimum/Maximum Query w/o Update
//
// Author      : Xin Wei Chow
// Date        : 06/06/2016
// Reliability : 1
// Tested On   : https://www.codechef.com/JUNE16/problems/CHSQARR
//
// Computes the min/max value in a 2D array
// Default is max, requires 2 changes to convert to max
//
// Complexity: O(nm lg(n)lg(m)) - build sparse table, O(1) - query
//
// Return:
//    int query(int i1, int j1, int i2, int j2)
//         Returns min/max value in rectangle (i1, j1) to (i2, j2)

#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef vector<vi> vvi;

int dp[11][11][1024][1024];

template<typename T>
class RMQ_2D {
private:
    vi lg;

public:
    RMQ_2D(const vector<vector<T>> &grid){
        int n = (int) grid.size(), m = (int) grid[0].size();
        lg.assign(max(n, m) + 1, 0);
        for (int i = 2; i <= max(n, m); i++)
            lg[i] = lg[i / 2] + 1;

        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                dp[0][0][i][j] = grid[i][j];

        for (int p = 0; (1 << p) <= n; p++){
            for (int q = 0; (1 << q) <= m; q++){
                if (p + q == 0) continue;
                for (int i = 0; i + (1 << p) <= n; i++){
                    for (int j = 0; j + (1 << q) <= m; j++){

                        if (p) dp[p][q][i][j] = max(dp[p-1][q][i][j],
                                                    dp[p-1][q][i+(1 << (p-1))][j]);
                        else dp[p][q][i][j] = max(dp[p][q-1][i][j],
                                                  dp[p][q-1][i][j+(1 << (q-1))]);
                    }
                }
            }
        }
    }

    T query(int i1, int j1, int i2, int j2){
        int p = lg[i2 - i1 + 1], q = lg[j2 - j1 + 1];
        i2 -= (1 << p) - 1;
        j2 -= (1 << q) - 1;
        return max({dp[p][q][i1][j1], dp[p][q][i1][j2],
                    dp[p][q][i2][j1], dp[p][q][i2][j2]});
    }
};

int main(){
    std::ios_base::sync_with_stdio(false);
    int n, m; cin >> n >> m;
    vector<vi> grid(n, vi(m));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
            cin >> grid[i][j];

    vector<vi> sum = grid;
    for (int i = 0; i < n; i++){
        for (int j = 0; j < m; j++){
            if (i) sum[i][j] += sum[i - 1][j];
            if (j) sum[i][j] += sum[i][j - 1];
            if (i && j) sum[i][j] -= sum[i - 1][j - 1];
        }
    }

    RMQ_2D<int> cc(grid);
    int q; cin >> q;
    while (q--){
        int h, w; cin >> h >> w;
        int res = INT_MAX;
        for (int i = 0; i + h <= n; i++){
            for (int j = 0; j + w <= m; j++){
                int i2 = i + h - 1, j2 = j + w - 1;
                int Ma = cc.query(i, j, i2, j2);
                int c = sum[i2][j2];
                if (i) c -= sum[i - 1][j2];
                if (j) c -= sum[i2][j - 1];
                if (i && j) c += sum[i - 1][j - 1];
                int cur = Ma * h * w - c;
                res = min(res, cur);
            }
        }
        cout << res << endl;
    }


    return 0;
}
