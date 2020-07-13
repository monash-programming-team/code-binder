#include <bits/stdc++.h>

using namespace std;

typedef vector<int> vi;

// Longest Increasing Subsequence
//
// Author      : Peter Whalan
// Date        : March 1, 2016
// Reliability : 1
// Tested On   : Codeforces group problem
//
// Finds the longest increasing subsequence
//
// Change lower_bound to upper_bound for longest non-decreasing subsequence.
// 
// B[i] is the smallest element that is at the end of an increasing sequence of
// length i+1 (for the part of A that has been processed). The line taking the
// lower_bound finds the first element not strictly less than A[i]. This means
// the element before it was less than A[i] so adding A[i] to this sequence we
// improve on the value at the iterator (or it stays the same).
//
// Although memory for B is preallocated, k keeps track of how many elements are
// being used. This corresponds to the current length of the longest increasing
// subsequence. lower_bound only ever uses elments that have been set so there
// is no need to initialise B to any particular value.
//
//
//
// Complexity: O( N log k )

int lis(const vi& A) {
	int k=0;
	vi B(A.size());
	for (int i=0;i<A.size();i++) {
		auto it=lower_bound(B.begin(),B.begin()+k,A[i]);
		*it=A[i];
		if (k==it-B.begin()) k++;
	}

	return k;
}

int s[1000000];

int main() {
	int n;
	long long a,b,c,d,e,f;
	cin >> n >> a >> b >> c >> d >> e >> f;
	s[0] = f;
	for(int i=1;i<n;i++)
		s[i] = (((a * s[i-1] + b) ^ c) + d) % e;

	cout << lis({s,s+n}) << ' ' << lis2({s,s+n}) << ' ';
	reverse(s,s+n);
	cout << lis({s,s+n}) << ' ' << lis2({s,s+n}) << endl;
}
