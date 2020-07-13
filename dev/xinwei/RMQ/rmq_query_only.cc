// Range Minimum/Maximum Query w/o Update
//
// Author      : Xin Wei Chow
// Date        : April 13, 2016
// Reliability : 5
// Tested On   : http://www.spoj.com/problems/RMQSQ/
//               http://www.spoj.com/problems/RPLN/
//               https://www.codechef.com/problems/MSTICK
//               UVA - 11491
//               2016 Anzac Round 3: G (Non-negative Partial Sums)
//               random cases - against Darcy's code-binder code
//
// Computes the min/max value in an array
// Default is min, requires 2 changes to convert to max
//
// Complexity: O(n lg n) - build sparse table, O(1) - query
//
// Return:
//    pair<T, int> query(int L, int R):
//         Returns a pair of (min/max value, index) from [L..R] inclusive

#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;
typedef vector<vi> vvi;

template<typename T>
class RMQ {
    vector<vector<pair<T, int>>> sptable;
    vi lg;

public:
    RMQ(vector<T> &v)  {
	int n = (int) v.size();
	lg.assign(n + 1, 0);
	for (int i = 2; i <= n; i++)
	    lg[i] = lg[i/2] + 1;

	sptable.assign(lg[n] + 1, vector<pair<T, int>>(n));
	for (int i = 0; i < n; i++)
	    sptable[0][i] = {v[i], i};

	for (int i = 1; i <= lg[n]; i++){
	    for (int j = 0; j + (1 << i) - 1 < n; j++){
		// change min to max for maximum query
		sptable[i][j] = min(sptable[i-1][j],
                                    sptable[i-1][j + (1 << (i-1))]);
	    }
	}
    }
    pair<T, int> query(int L, int R){
	int k = lg[R - L + 1];
	// change min to max for maximum query
	return min(sptable[k][L], sptable[k][R - (1 << k) + 1]);
    }
};

int main(){
    vi v = { 1, 3, 2, 7, 10, -5 };
    RMQ<int> cc(v);
    cout << cc.query(0, 5).second << endl; // 5
    cout << cc.query(0, 3).second << endl; // 0
    cout << cc.query(1, 2).second << endl; // 2


    return 0;
}
