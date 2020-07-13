/**
 * Unweighted Bipartite Matching - Hopcroft–Karp
 * @author Ryan Hechenberger
 * @date 2016-07-31
 * Reliability: 5
 * Tested On: uva10092, uva10080, uva11418, spoj_MATCHING, anzac2_2016
 *
 * Matches groups of nodes on the lest and the right to get a maximal matching.
 * This uses the Hopcroft–Karp approch. The left and right size are both
 * dynamically assigned, you do not explititly state them. The class itself
 * holds the matching of the left size.
 *
 * Complexity: O ( sqrt(R) * E ), where R = right no. of nodes, E = no. edges
 */
#include <bits/stdc++.h>
#define me (*this)

typedef vector<int> vi;

/**
 * The structure that handles the matching using the Ford-Fulkerson approch.
 * Principle from: https://en.wikipedia.org/wiki/Hopcroft%E2%80%93Karp_algorithm
 */
class BiMatchFord : public vi
{
	int p;
	vector<vi> R;
	vi D;
	queue<int> Q;
	list<int> umat;
	
public:
	/// adds a connection between l node and r node
	void add(int l, int r) {
		if (l >= (int)size()) resize(l+1, -1);
		if (++r >= (int)R.size()) R.resize(r+1);
		R[r].push_back(l);
	}

	int match() {
		int c = 0; D.assign(R.size(), 0);
		for (int i = 1; i < (int)R.size(); ++i) umat.push_back(i);
		for (p = 0; bfs(); p = D[0]+1)
		for (auto it = umat.begin(), ite = umat.end(); it != ite; )
			if (dfs(*it)) { c++; umat.erase(it++); } else ++it;
		return c;
	}

private:
	bool bfs() {
		for (int u : umat) { D[u] = p; Q.push(u); }
		while (!Q.empty()) {
			int u = Q.front(); Q.pop();
			if (D[u] != D[0])
				for (int v : R[u])
					if (D[me[v]+1] < p) { D[me[v]+1] = D[u]+1; Q.push(me[v]+1); }
		}
		return D[0] >= p;
	}
	
	bool dfs(int u) {
		if (u != 0) {
			for (int v : R[u])
				if (D[me[v]+1] == D[u]+1 && dfs(me[v]+1))
					{ me[v] = u-1; return true; }
			D[u] = D[0];
			return false;
		}
		return true;
	}

	// Optional functions from here
	/// Number of nodes on left side
	int left() const { return (int)size(); }

	/// Number of nodes on right side
	int right() const { return (int)R.size(); }

	/// Returns the matching pairs for the right side
	vi opposite() const {
		vi m(R.size(), -1);
		for (int i = 0; i < (int)size(); ++i) if (me[i] >= 0) m[me[i]] = i;
		return m;
	}
};
