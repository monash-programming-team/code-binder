#include <bits/stdc++.h>

using namespace std;

typedef vector<int> vi;

// Fenwick Tree
//
// Author      : Peter Whalan (Modified by Daniel)
// Date        : February 29, 2016
// Reliability : 1
// Tested On   : ASC18:E
//
// Performs range updates and point queries on array
// 
// Treats arguments as 0-based
//
// Every element is contained within O(logn) intervals. The value for an
// element is determined by summing the values of these intervals. It is
// not guaranteed that the values for the intervals won't overflow even
// if you know the result doesn't.
// This probably doesn't matter for integers because the integers will
// just overflow back as the sum is taken.
//
// Complexity:
//		O( N ) to build
//		O( log N ) to update and query

//listings:ft
// Fenwick tree with ranged updates and point queries. Complexity: O(log(n))
template<typename T> struct FenwickTree {
	int N;  vector<T> A;
	FenwickTree(int n): N(n+1), A(N) {}                    // Create tree with n elements
	void adjust(int b, T v) { for (;b;b-=b&-b) A[b]+=v; }              // Add v to A[0,b)
	void adjust(int a,int b, T v) { adjust(b,v), adjust(a,-v); }       // Add v to A[a,b)
	T pq(int i) { T r=0; for (i++;i<N;i+=i&-i) r+=A[i]; return r;	}           // Get A[i]
};
//listings:/ft

int main() {}
