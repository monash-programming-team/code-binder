// Coin Change
//
// Author      : Daniel Anderson
// Date        : 3 May 2016
// Reliability : 1
// Tested On   : MCPC2016 F
//
// Finds the smallest way that change can be made for all amounts up and
// including the given maximum price using the given denominations.
// 
// Complexity: O(max(NM, M^2)) -- N = number of coins, M = max price

#include<bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;

const int INF = INT_MAX / 2;

// Computes the minimum number of coins needed to make all values up to max_price
// using the given coins. Returns INF if a given value can not be made.
vi min_coins(const vi& coins, int max_price) {
	int n = coins.size();
	vvi DP(n + 1, vi(max_price + 1, INF));
	for (int i = 0; i <= n; i++) DP[i][0] = 0;
	for (int i = 1; i <= n; i++) {
		for (int c = 1; c <= max_price; c++) {
			if (c < coins[i - 1]) DP[i][c] = DP[i - 1][c];
			else DP[i][c] = min(DP[i - 1][c], DP[i][c - coins[i - 1]] + 1);
		}
	}
	return DP[n];
}

// Computes the minimum number of coins and the configuration to make all values
// up to max_price using the given coins
vector<pair<int, vi>> coin_change(const vi& coins, int max_price) {
	int n = coins.size();
	vvi DP(n + 1, vi(max_price + 1, INF));
	vvi used(n + 1, vi(max_price + 1, 0));
	for (int i = 0; i <= n; i++) DP[i][0] = 0;
	for (int i = 1; i <= n; i++) {
		for (int c = 1; c <= max_price; c++) {
			if (c < coins[i - 1]) DP[i][c] = DP[i - 1][c];
			else if (DP[i][c - coins[i - 1]] < DP[i - 1][c])
				DP[i][c] = DP[i][c - coins[i - 1]] + 1, used[i][c] = 1;
			else DP[i][c] = DP[i - 1][c];
		}
	}
	vector<pair<int, vi>> ans(max_price+1);
	for (int c = 1; c <= max_price; c++) {
		int i = n, j = c;
		ans[c].first = DP[n][c];
		while (i) {
			if (!used[i][j]) i--;
			else ans[c].second.push_back(coins[i - 1]), j -= coins[i - 1];
		}
	}
	return ans;
}

// Usage example
int main() {
	vi denominations = {1, 3, 5, 10};
	auto cnt = min_coins(denominations, 10);
	auto change = coin_change(denominations, 10);
	
	return 0;
}