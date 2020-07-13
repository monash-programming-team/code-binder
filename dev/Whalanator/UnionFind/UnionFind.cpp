#include <bits/stdc++.h>

using namespace std;

typedef vector<int> vi;

// UnionFind
//
// Author      : Peter Whalan
// Date        : March 1, 2016
// Reliability : 1
// Tested On   : CF 687D
//
// Union sets and query if elements are in the same set.
//
// Complexity: Amortized O( log N ) per query

//listings:union_find
// Union-Find with path compression. Complexity: O(log(N)) amortized per query.
struct UnionFind {
	vi A;
	UnionFind(int N): A(N) { iota(A.begin(),A.end(),0); }
	int find(int i) { return A[i]==i?i:A[i]=find(A[i]); }
	bool same(int i,int j) { return find(i)==find(j); }
	bool merge(int i,int j) { return same(i,j)?0:A[A[i]]=A[j],1; }
};
//listings:/union_find

int main() {
  UnionFind uf(10);
  uf.merge(0,9);
}
