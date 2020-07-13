// Hungarian Algorithm test on SPOJ - Baby
#include <bits/stdc++.h>
using namespace std;

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

int dist(int x1, int y1, int x2, int y2) {
	return abs(x1 - x2) + abs(y1 - y2);
}

int main() {
	int n; 
	while (cin >> n, n > 0) {
		vi a(n), b(n);
		
		for (int i = 0; i < n; i++) cin >> a[i];
		for (int i = 0; i < n; i++) cin >> b[i];
		
		vvi A(n, vi(n));
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				A[i][j] = dist(i, a[i], j, b[j]);
			}
		}
		
		auto res = hungarian(A);
		cout << res.first << endl;
	}
}