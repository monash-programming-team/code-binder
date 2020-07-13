// The Hungarian Algorithm
//
// Solves the weighted assignment problem, or in graph terms, the minimum
// weight maximum bipartite matching problem.
//
// Author      : Daniel Anderson
// Date        : 15/11/2016
// Reliability : 5
// Tested on   : SPOJ - Baby, SPOJ - Selfish Cities, CF491C, UVA10746, UVA12637
//
// Usage:
//  pair<T, vi> hungarian(A)
//   Computes the minimum weighted assignment given the weights A[n][m]
//   IMPORTANT: n <= m or you will infinite loop. Returns the minimum weight
//   and a vector matching[n] containing each persons assigned match. If you
//   want a maximum weight matching, just make all of the weights negative.
//
#include <bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef pair<int,int> pii;
typedef vector<int> vi;
typedef vector<vi> vvi;

template<typename T>
pair<T, vi> hungarian(const vector<vector<T>>& A) {
	int n = (int) A.size(), m = (int) A[0].size();
	T inf = numeric_limits<T>::max() / 2;
	vi way(m + 1), p(m + 1), used(m + 1);
	vector<T> u(n + 1), v(m + 1), minv(m + 1);
	for (int i = 1; i <= n; i++) {
		p[0] = i;  int j0 = 0;
		minv.assign(m + 1, inf),  used.assign(m + 1, 0);
		do {
			int i0 = p[j0], j1 = 0; T delta = inf;
			used[j0] = true;
			for (int j = 1; j <= m; j++) 
				if (!used[j]) {
					T cur = A[i0 - 1][j - 1] - u[i0] - v[j];
					if (cur < minv[j]) {
						minv[j] = cur;
						way[j] = j0;
					}
					if (minv[j] < delta) {
						delta = minv[j];
						j1 = j;
					}
				}
			for (int j = 0; j <= m; j++) 
				if (used[j]) {
					u[p[j]] += delta;
					v[j] -= delta;
				} else minv[j] -= delta;
			j0 = j1;
		} while (p[j0] != 0);
		do {
			int j1 = way[j0];
			p[j0] = p[j1];
			j0 = j1;
		} while (j0 != 0);
	}
	vi matching(n);
	for (int i = 1; i <= m; i++) {
		if (p[i] > 0) matching[p[i] - 1] = i - 1;
	}
	return {-v[0], matching};
}

const double EPS = 1e-8;

int dist(const pii& a, const pii& b) {
	return abs(a.X - b.X) + abs(a.Y - b.Y);
}

int main() {
	int R, N;
	while (cin >> R >> N) {
		vvi A(N, vi(R, 0));
		vector<pii> dr(R), dn(N);
		for (int i = 0; i < R; i++) cin >> dr[i].X >> dr[i].Y;
		for (int j = 0; j < N; j++) cin >> dn[j].X >> dn[j].Y;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < R; j++)
				A[i][j] = dist(dr[j], dn[i]);
		auto res = hungarian(A);
		cout << res.first << '\n';
	}
}
