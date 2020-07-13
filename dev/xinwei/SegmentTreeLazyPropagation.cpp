// Segment Tree with Lazy Propagation
// 
// Time Complexity: O(n, log n)
// Memory: 4 * n
#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;

class SegmentTree {
private:
    int n;
    vi arr, st, lazy;
    int merge(int p1, int p2){
	return max(p1, p2);
    }
    void prop(int p, int L, int R){
	if (lazy[p]){
	    st[p] = lazy[p];
	    if (L != R){
		lazy[p * 2] = lazy[p * 2 + 1] = lazy[p];
	    }
	}
	lazy[p] = 0;
    }
    void build(int p, int L, int R){
	if (L == R) st[p] = arr[L];
	else {
	    int mid = (L + R) / 2;
	    build(p * 2, L, mid);
	    build(p * 2 + 1, mid + 1, R);
	    st[p] = merge(st[p * 2], st[p * 2 + 1]);
	}
    }
    int query(int p, int L, int R, int i, int j){
	prop(p, L, R);
	if (i > R || j < L) return 0;
	if (i <= L && j >= R) return st[p];
	int mid = (L + R) / 2;
	int p1 = query(p * 2, L, mid, i, j);
	int p2 = query(p * 2 + 1, mid + 1, R, i, j);
	return merge(p1, p2);
    }
    void update(int p, int L, int R, int i, int j, int val){
	prop(p, L, R);
	if (i > R || j < L) return;
	if (i <= L && j >= R) {
	    lazy[p] = val;
	    prop(p, L, R);
	    return;
	}
	int mid = (L + R) / 2;
	update(p * 2, L, mid, i, j, val);
	update(p * 2 + 1, mid + 1, R, i, j, val);
	st[p] = merge(st[p * 2], st[p * 2 + 1]);
    }

public:
    SegmentTree(int n, const vi &arr) : n(n), arr(arr) {
	st.assign(4 * n, 0);
	lazy.assign(4 * n, 0);
	build(1, 0, n-1);
    }
    int query(int i, int j){
	return query(1, 0, n-1, i, j);
    }
    void update(int i, int j, int val){
	update(1, 0, n-1, i, j, val);
    }
};

int main(){
  	int n = 5;
  	vi arr = {1, 2, 7, 4, 5};
  	SegmentTree segtree(n, arr);
  	cout << segtree.query(0, 4) << endl; // returns 7
  	cout << segtree.query(3, 4) << endl; // returns 5

  	return 0;
}