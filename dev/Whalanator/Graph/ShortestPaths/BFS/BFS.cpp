#include <bits/stdc++.h>

#define x first
#define y second

using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int,int> ii;
typedef vector<ii> vii;
typedef vector<vii> vvii;

// BFS
//
// Author      : Peter Whalan
// Date        : June 29, 2016
// Reliability : 3
// Tested On   : UVA 627,336 CF 687E
//
// Single source shortest paths
//
// al is (0-based) adjacency list. al[i][j] is vertex connected to vertex i.
// s is source.
//
// Second version provides predecessor array with pred[s]==-1. Anything
// unreachable also has pred of -1. Pair of vectors is returned. dist first,
// pred second.
//
// Complexity O(E + V)

vi bfs(const vvi& al, int s) {
	queue<int> q;
	vi dist(al.size(),INT_MAX);
	dist[s]=0;
	q.push(s);
	while (!q.empty()) {
		int i=q.front();q.pop();
		for (int j:al[i]) if (dist[j]==INT_MAX) {
			dist[j]=dist[i]+1;
			q.push(j);
		}
	}
	return dist;
}

pair<vi,vi> bfs(const vvi& al, int s) {
	queue<int> q;
	vi dist(al.size(),INT_MAX),pred(al.size(),-1);
	dist[s]=0;
	q.push(s);
	while (!q.empty()) {
		int i=q.front();q.pop();
		for (int j:al[i]) if (dist[j]==INT_MAX) {
			dist[j]=dist[i]+1;
			pred[j]=i;
			q.push(j);
		}
	}
	return {dist,pred};
}

int main() {

}
