/**
 * Unweighted Bipartite Matching - Ford–Fulkerson
 * @author Ryan Hechenberger
 * @date 2016-07-31
 * Reliability: 5
 * Tested On: uva10092, uva10080, uva11418, spoj_MATCHING, anzac2_2016
 * 
 * Matches groups of nodes on the lest and the right to get a maximal matching.
 * This uses the Ford-Fulkerson approch. The left and right size are both
 * dynamically assigned, you do not explititly state them. The class itself
 * holds the matching of the left size.
 * 
 * Complexity: O ( R * E ), where R = right no. of nodes, E = no. edges
 */
#include <bits/stdc++.h>
#define me (*this)

typedef vector<int> vi;

/**
 * The structure that handles the matching using the Ford-Fulkerson approch.
 * Principle from: http://stackoverflow.com/questions/563198/
 */
struct BiMatching : vi
{
private:
	int p;
	vector<vi> R;
	vi seen;
	
	bool dfs(int u) {
		seen[u] = p;
		for (int v : R[u])
			if (me[v] < 0 || ( seen[me[v]] != p && dfs(me[v]) ))
			{ me[v] = u; return true; }
		return false;
	}
	
public:
	/// adds a connection between l node and r node
	void add(int l, int r) {
		if (l >= (int)size()) resize(l+1, -1);
		if (r >= (int)R.size()) R.resize(r+1);
		R[r].push_back(l);
        }
	int match() {
		int c = 0; seen.assign(R.size(), -1);
		for (p = 0; p < (int)R.size(); ++p)
			c += dfs(p);
		return c;
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
