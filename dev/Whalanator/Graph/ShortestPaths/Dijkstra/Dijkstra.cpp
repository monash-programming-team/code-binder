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
// al is (0-based) adjacency list. al[i][j] is a pair. The first element
// is the vertex and the second element is the weight.
// s is source.
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


vi dijkstra(const vvii& al, int s) {
	priority_queue<ii,vector<ii>,greater<ii> > pq;
	vi dist(al.size(),INT_MAX);
	pq.emplace(0,s);
	int u,d;
	while (!pq.empty()) {
		tie(d,u)=pq.top();pq.pop();
		if (dist[u]<=d) continue;
		dist[u]=d;//if u==t return d (heuristic)
		for (ii v:al[u])
			pq.emplace(d+v.y,v.x);
	}
	return dist;
}

pair<vi,vi> dijkstra(const vvii& al, int s) {
	priority_queue<ii,vector<ii>,greater<ii> > pq;
	vi dist(al.size(),INT_MAX),pred(al.size(),-1);
	dist[s]=0;
	pq.emplace(0,s);
	int u,d;
	while (!pq.empty()) {
		tie(d,u)=pq.top();pq.pop();
		if (dist[u]<d) continue;//if u==t return d (heuristic)
		for (ii v:al[u]) if (d+v.y<dist[v.x]) {
			dist[v.x]=d+v.y;
			pred[v.x]=u;
			pq.emplace(d+v.y,v.x);
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
	cin >> N;
	for (int cas=1;cas<=N;cas++) {
		printf("Case #%d: ",cas);
		cin >> n >> m >> S >> T;
		int i,j,k;
		al.assign(n,vii());
		for (int c=0;c<m;c++) {
			cin >> i >> j >> k;
			al[i].emplace_back(j,k);
			al[j].emplace_back(i,k);
		}

		int r=dijkstra(al,S)/*.x*/[T];
		if (r==INT_MAX) cout << "unreachable\n";
		else cout << r << endl;
	}
}
