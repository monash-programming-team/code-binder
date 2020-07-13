#include <bits/stdc++.h>

using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;

// Find Euler Path/Tour
//
// Author      : Peter Whalan
// Date        : August 25, 2016
// Reliability : 1
// Tested On   : UVA 10054
//
// Finds the Euler path/tour of a graph (if one exists).
//
// A linked list is avoided by putting the vertex on the Euler tour after making
// the recursive call.
//
// et contains the vertices of the Euler path/tour after running the function.
// For an Euler tour the starting vertex is not put at the end of et. For an
// Euler path both starting and ending vertices are included.
//
// Complexity: O( N + M )

struct Eulerian {
	struct Edge {
		int j;
		bool used;
	};
	
	// Number of vertices, start vertex, end vertex
	int N,s,t;
	vvi al;
	vector<Edge> el;
	vi indeg,ct,et; // Vector of in degrees, counters and Euler path/tour

	Eulerian(int n): N(n), al(N), indeg(N), ct(N) {}

	// Add edge directed
	void addedir(int i,int j) {
		al[i].push_back(el.size());
		el.push_back({j,0});
		indeg[j]++;
	}

	// Add edge undirected
	void addeundir(int i,int j) {
		addedir(i,j);
		addedir(j,i);
	}

	// Find Euler path/tour treating graph as directed
	// Returns:
	//	1 if there was an Euler path
	//	2 if there was an Euler tour
	//	0 otherwise
	int finddir() {
		s=t=-1;
		for (int i=0;i<N;i++)
			if (al[i].size()==indeg[i]+1) {
				if (s==-1) s=i;
				else return 0;
			}
			else if (al[i].size()+1==indeg[i]) {
				if (t==-1) t=i;
				else return 0;
			}
			else if (al[i].size()!=indeg[i]) return 0;
		if (s==t) s=t=el[0].j;// This segfaults if there are no edges
		find(s);
		reverse(et.begin(),et.end());
		if (et.size()!=el.size()) return 0; //Didn't use every edge
		if (s!=t) {
			et.push_back(t);
			return 1;
		}
		return 2;
	}
	// Find Euler path/tour treating graph as undirected
	// Returns:
	//	1 if there was an Euler path
	//	2 if there was an Euler tour
	//	0 otherwise
	int findundir() {
		s=t=-1;
		for (int i=0;i<N;i++)
			if (al[i].size()%2) {
				if (s==-1) s=i;
				else if (t==-1) t=i;
				else return 0;
			}
		if (s==t) s=t=el[0].j;
		find(s);
		reverse(et.begin(),et.end());
		if (2*et.size()!=el.size()) return 0; //Didn't use every edge
		if (s!=t) {
			et.push_back(t);
			return 1;
		}
		return 2;
	}

	void find(int i) {
		while (ct[i]<al[i].size()) {
			int k=al[i][ct[i]++];
			if (el[k].used) continue;
			el[k].used = el[k^1].used = 1;//Only assign to el[k].used if graph
			find(el[k].j);                //is directed.
			et.push_back(i);
			// At this point we are returning from j to i over edge
			// el[al[i][ct[i]-1]. Can store the edge if you need it. Reverse the
			// vector of edges afterwards to put them in the right order.
		}
	}
};

int T,N;

int main() {
	cin >> T;
	int i,j;
	for (int cas=1;cas<=T;cas++) {
		if (cas>1) cout << endl;
		printf("Case #%d\n",cas);
		cin >> N;
		Eulerian eu(51);
		//vvi al(51);
		//vector<Edge> el(N<<1);
		for (int c=0;c<(N<<1);c+=2) {
			cin >> i >> j;
			eu.addeundir(i,j);
			/*el[c].j=j;
			el[c+1].j=i;
			al[i].push_back(c);
			al[j].push_back(c+1);*/
		}
		//ct.assign(51,0);
		//et.clear();

		int r =eu.findundir();
		//findEulertour(el[0].j,al,el);
		//reverse(et.begin(),et.end());

		//for (int c=0;c<51;c++) if (al[c].size()%2) goto fail;
		if (r!=2) goto fail;

		for (int i=0;i<N;i++) printf("%d %d\n",eu.et[i],eu.et[(i+1)%N]);
		continue;
fail:
		printf("some beads may be lost\n");
	}
}
