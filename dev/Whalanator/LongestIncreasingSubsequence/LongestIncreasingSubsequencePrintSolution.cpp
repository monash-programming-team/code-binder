#include <bits/stdc++.h>

using namespace std;

typedef vector<int> vi;

// Longest Increasing Subsequence + Print Solution
//
// Author      : Peter Whalan
// Date        : March 1, 2016
// Reliability : 0
// Tested On   : 
//
// Finds the longest increasing subsequence and gives the sequence as a
// vector.
//
// Change lower_bound to upper_bound for longest non-decreasing subsequence.
//
// Complexity: O( N log N )

vi lis(const vi& A) {
	int k=0,N=A.size();
	int is[N],vs[N],pred[N];
	for (int i=0;i<N;i++) {
		int j=lower_bound(vs,vs+k,A[i])-vs;
		pred[i]=j?is[j-1]:-1;
		vs[j]=A[is[j]=i];
		k=max(k,1+j);
	}

	vi r(k);
	for (int i=is[k-1],j=k-1;i!=-1;i=pred[i],j--) r[j]=i;

	return r;
}

int main() {
}
