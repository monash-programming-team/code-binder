// Hungarian Algorithm test on Codeforces 491C
// This test requires you to print the matching
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

char c(int k) {
	if (k < 26) return 'a' + k;
	else return 'A' + k - 26;
}

int idx(char a) {
	if ('a' <= a && a <= 'z') return a - 'a';
	else return (a - 'A') + 26;
}

int main() {
	int N, K;
	string enc, ans;
	cin >> N >> K >> enc >> ans;
	
	vvi A(K, vi(K, 0));

	for (int i = 0; i < N; i++) {
		A[idx(enc[i])][idx(ans[i])]--;
	}
	
	auto res = hungarian(A);
	
	cout << -res.first << endl;
	for (int k = 0; k < K; k++) {
		cout << c(res.second[k]);
	}
	cout << endl;
}
