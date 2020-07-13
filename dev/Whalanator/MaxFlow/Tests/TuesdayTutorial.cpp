#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

//see CF 653D
// Ford Fulkerson (DFS) Max Flow
//
// Author      : Peter Whalan
// Date        : September 1, 2016
// Reliability : 2
// Tested On   : CF:653/D, Tuesday Tutorial
//
// Computes the max flow of a graph using Ford Fulkerson's method, implemented
// with DFS.
//
// This code works if the value of a residual edge cannot exceed an int. For a
// directed graph, this is the case if each edge in the input fits in an int.
// For an undirected graph this is the case if 2 * (max undirected capacity)
// fits in an int.
//
// Remember that the residual graph has twice as many edges as the input graph.
//
// There will be an infinite loop if s == t.
//
// Complexity: O((V+E)F) where F is the max flow

// Represents an edge going from some vertex i, to vertex j, with capacity C.
struct Edge {
	int j,C; // Incoming vertex (ie. el[ al[i][c] ] != i), Capacity
};

// Your code must set the variables in this block. al and el are for the
// residual graph.
int n,s,t;// Number of vertices, Source, Sink
vector<Edge> el; // Edge list. Ensure pairs of edges are together (xor trick).
vvi al; // Adjacency list. Stores indices of outgoing edges in edge list.

vector<bool> vis;

int dfs(int i,int cf) {
	if (vis[i]) return 0;
	vis[i]=1;
	if (i==t) return cf;
	for (int c:al[i]) {
		int ncf = min(cf,el[c].C),f;
		if (ncf && (f=dfs(el[c].j,ncf))) {
			el[c].C-=f;
			el[c^1].C+=f;
			return f;
		}
	}

	return 0;
}

ll maxflow() {
	ll Mf=0,f=1;
	while (f) {
		vis.assign(n,0);
		f=dfs(s,INT_MAX); // Change to LLONG_MAX if single flow can exceed it.
		Mf+=f;
	}
	return Mf;
}

int main() {
	int m;
	scanf("%d%d%d%d",&n,&m,&s,&t);
	s--;t--;
	al.assign(n,{});
	el.resize(2*m);
	for(int i=0;i<m;i++){
		int u,v,c;
		scanf("%d%d%d",&u,&v,&c);
		u--;v--;
		al[u].push_back(2*i);
		al[v].push_back(2*i+1);
		el[2*i]={v,c};
		el[2*i+1]={u,0};//change 0 to c for undirected egde
	}

	cout << maxflow() << endl;
}
