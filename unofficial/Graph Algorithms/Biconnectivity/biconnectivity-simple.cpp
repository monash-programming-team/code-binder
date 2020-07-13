// Biconnectivity Algorithms - Articulation points and bridges
//
// Author      : Daniel Anderson
// Date        : March 30, 2016
// Reliability : 5
// Tested On   : SPOJ - SUBMERGE (Articulation points),
//				 UVA - 315 (Articulation points),
//				 UVA - 10199 (Articulation points),
//				 UVA 796 (Bridges),
//				 Codeforces 652E (Bridges)
//				 UVA - 610 (Bridges & Biconnected Components),
// 				 UVA - 10765 (Biconnected Components)
//				 
// Finds biconnected components, articulation points and bridges of
// an undirected graph.
//
// NOTE: This code only works for simple graphs! i.e. no multiple edges
// between the same vertices or edges from a vertex to itself.
//
// Definitions:
//  Articulation point: 	A vertex whose removal would increase the number
//							of connected components in the graph.
//  Bridge:					An edge whose removal would increase the number
//							of connected components in the graph.
//  Biconnected Graph:		A graph is biconnected if it contains no
//							articulation points.
//	Biconnected Components: The biconnected components of a graph are the
//							maximal biconnected sub-graphs.
//  Edge-Biconnected Graph: A graph is edge-biconnected if it contains no
//							bridges
//  Edge-Biconnected Components: The edge-biconnected components of a graph
//							are the maximal edge-biconnected sub-graphs
//
// NOTE: If you want bridge-connected components (maximal sub-graphs
//		 containing no bridges), use this algorithm to find the bridges,
//		 remove them from the graph and then find the connected
//		 components of the remaining graph
//
// Complexity: O(V + E)

#include<bits/stdc++.h>
using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int, int> pii;

class biconnectivity {
	int n;		// Number of nodes in the graph
	vvi edges;	// Edge adjacency list
	
	int dfs_root, dfs_count, root_children;
	vi dfs_low, dfs_num, dfs_parent;
	vi cut_points, is_cut_point;		
	vector<pii> bridge;
	vector<vector<pii>> bccs;
	vector<pii> comp;		

	void dfs(int u) {
		dfs_low[u] = dfs_num[u] = dfs_count++;
		for (int v : edges[u]) {
			pii edge = minmax(u,v);
			if (dfs_num[v] == -1) {
				dfs_parent[v] = u;
				if (u == dfs_root) root_children++;
				comp.push_back(edge);
				dfs(v);
				if (dfs_low[v] >= dfs_num[u]) {
					is_cut_point[u] = true;
					bccs.emplace_back(comp.rbegin(), find(comp.rbegin(),
										comp.rend(),edge)+1);
					comp.resize(comp.size()-bccs.back().size());
				}
				if (dfs_low[v] > dfs_num[u]) bridge.push_back(edge);
				dfs_low[u] = min(dfs_low[u], dfs_low[v]);
			}
			else if (v != dfs_parent[u]) {
				dfs_low[u] = min(dfs_low[u], dfs_num[v]);
				if (dfs_num[v] < dfs_num[u]) comp.push_back(edge);
			}
		}
	}
	
	void search() {
		cut_points.clear(); bridge.clear(); bccs.clear(); comp.clear();
		is_cut_point.assign(n, 0); dfs_parent.assign(n, -1);
		dfs_num.assign(n, -1);  dfs_low.assign(n, 0); dfs_count = 0;
		for (int v = 0; v < n; v++) {
			if (dfs_num[v] == -1) {
				dfs_root = v; root_children = 0;
				dfs(v);
				is_cut_point[v] = (root_children > 1);
			}
		}
		for (int v = 0; v < n; v++) if (is_cut_point[v]) cut_points.push_back(v);
	}
	
public:
	// Create a graph
	//
	// Parameters:
	//  n	: The number of nodes in the graph
	biconnectivity(int n) : n(n), edges(n) { }
	
	// Add an undirected edge between vertex u and vertex v
	//
	// Parameters:
	//  u	: The (0-based) index of the first vertex
	//	v	: The (0-based) index of the second vertex
	void add_edge(int u, int v) { edges[u].push_back(v); edges[v].push_back(u); }
	
	// Finds all of the articulation points of the graph
	//
	// Return:
	//  The indices of the articulation points (vertices) of the graph
	vi articulation_points() { search(); return move(cut_points); }
	
	// Finds all of the bridges (cut edges) of the graph
	//
	// Return:
	//  The bridges of the graph as vertex pairs representing the
	//  start and end points of the edge with start < end
	vector<pii> bridges() {	search(); return move(bridge);	}
	
	// Finds all of the biconnected components of the graph
	//
	// Return:
	//  The biconnected components of the graph. The biconnected
	//  components are represented by the edges contained in them
	vector<vector<pii>> biconnected_components() { search(); return move(bccs);	}
	
	// Finds all three of the above
	//
	// Return:
	//  All of the above
	tuple<vi, vector<pii>, vector<vector<pii>>> get_all() {
		search();
		return make_tuple(move(cut_points), move(bridge), move(bccs));
	}
};


// ---------------------------------------------------
// Biconnectivity ends. Tests start.
// --------------------------------------------------

int main() {
	biconnectivity b(6);
	b.add_edge(1,0);
	b.add_edge(1,2);
	b.add_edge(1,3);
	b.add_edge(1,4);
	b.add_edge(1,5);
	b.add_edge(4,5);
	
	//If you want just one feature
	//Will recompute for each request even if inefficient
	vi articulation_points = b.articulation_points();
	auto bridges = b.bridges();
	auto biconnected_components = b.biconnected_components();
	
	cout << "Articulation Points: " << articulation_points.size() << endl;
	cout << "Bridges: " << bridges.size() << endl;
	cout << "Biconnected Components: " << biconnected_components.size() << endl;
	
	biconnectivity b2(6);
	b2.add_edge(0,1);
	b2.add_edge(1,2);
	b2.add_edge(0,2);
	b2.add_edge(3,5);
	b2.add_edge(4,5);
	b2.add_edge(3,4);
	b2.add_edge(2,4);
	
	//Get everything at once
	tie(articulation_points, bridges, biconnected_components) = b2.get_all();
	
	cout << "Articulation Points: " << articulation_points.size() << endl;
	cout << "Bridges: " << bridges.size() << endl;
	cout << "Biconnected Components: " << biconnected_components.size() << endl;
}
