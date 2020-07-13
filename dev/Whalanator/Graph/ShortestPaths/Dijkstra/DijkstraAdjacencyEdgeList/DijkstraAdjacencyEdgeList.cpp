#include <bits/stdc++.h>

#define x first
#define y second

using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> ii;
typedef vector<ii> vii;
typedef vector<vii> vvii;

// Dijkstra
//
// Author      : Peter Whalan
// Date        : June 15, 2016
// Reliability : 3
// Tested On   : UVA 10986,762, CERC 12 F
//
// Single source shortest paths
//
// al is (0-based) adjacency list. al[i][j] is an int. el[al[i][j]] is an
// outgoing edge of vertex i. See edge structure below.
//
// Second version provides predecessor array with pred[s]==-1. Anything
// unreachable also has pred of -1. Pair of vectors is returned. dist first,
// pred second.
//
// path takes the result of running the second version of dijkstra and a
// destination vertex and produces a vector of vertices from the source to
// the destination along the shortest path. The source and destination are
// included in this path. Each element in the path is distinct. The vector
// returned contains only the destination vertex if there is no path. Note
// this can also happen if the source is the destination.
//
// Complexity O(ElogE + V)

struct Edge {
	int j,w; // Connecting vertex, weight
};

const int inf=INT_MAX/2;

vi dijkstra(const vvi& al, const vector<Edge>& el, int s) {
	priority_queue<ii,vector<ii>,greater<ii>> pq;
	vi dist(al.size(),inf);
	pq.emplace(0,s);
	int u,d;
	while (!pq.empty()) {
		tie(d,u)=pq.top();pq.pop();
		if (dist[u]<=d) continue;
		dist[u]=d;//if u==t return d (heuristic)
		for (int c:al[u])
			pq.emplace(min(inf,d+el[c].w),el[c].j);
	}
	return dist;
}

pair<vi,vi> dijkstra(const vvi& al, const vector<Edge>& el, int s) {
	priority_queue<ii,vector<ii>,greater<ii>> pq;
	vi dist(al.size(),inf),pred(al.size(),-1);
	dist[s]=0;
	pq.emplace(0,s);
	int u,d;
	while (!pq.empty()) {
		tie(d,u)=pq.top();pq.pop();
		if (dist[u]<d) continue;//if u==t return d (heuristic)
		for (int c:al[u]) {
			int j=el[c].j,w=el[c].w;
			if (d+w<dist[j]) {
				dist[j]=d+w;
				pred[j]=u;
				pq.emplace(d+w,j);
			}
		}
	}
	return {dist,pred};
}

vi path(const vi& pred,int t) {
	vi r;
	for (;t!=-1;t=pred[t]) r.push_back(t);
	return vi(r.rbegin(),r.rend());
}

int N,n,m,S,T;
vvii al;

int main() {
	
}
