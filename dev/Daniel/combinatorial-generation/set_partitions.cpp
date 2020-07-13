// Set partition generation
//
// Author: Daniel Anderson
// Date: 20-01-2017
// Reliability: 1
// Tested on: Brute-force
//
// Reference: Knuth, The Art of Computer Programming, Volume 4A
//
#include<bits/stdc++.h>
using namespace std;

typedef long long int ll;

typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:set_partitions
// Generates set partitions for a set of size n in gray code order. A set partition
// is represented as a vector of size n where P[i] = the index of the set that
// element i belongs to. Safe to use for n <= 13, where B_n = 27 million.
struct set_partition_generator {
	int n;  vi a, b, d;  bool done = false;
	set_partition_generator(int n) : n(n), a(n), b(n, 1), d(n, 1) { }
	void fix(int j, int m) { fill(b.begin() + j + 1, b.end(), m); }
	bool has_next() { return !done; }
	vi next_partition() {
		vi ans = a;  int j = n - 1;
		while (a[j] == d[j]) d[j--] ^= 1;
		if (j == 0) done = true; 
		else if (d[j] != 0) {
			if (a[j] == 0) { a[j] = b[j];	fix(j, a[j] + 1);	}
			else if (a[j] == b[j]) { a[j] = b[j] - 1;	fix(j, b[j]);	}
			else a[j]--;
		}
		else {
			if (a[j] == b[j] - 1) {	a[j] = b[j]; fix(j, a[j] + 1); }
			else if (a[j] == b[j]) { a[j] = 0; fix(j, b[j]); }
			else a[j]++;
		}
		return ans;
	}
};
//listings:/set_partitions

namespace brute_force_testing {
  int MAXN = 14;
  int bell[] = {1, 1, 2, 5, 15, 52, 203, 877, 4140, 21147, 115975, 678570, 4213597, 27644437, 190899322};
  void test() {
    for (int n=1; n<=MAXN; n++) {
      cout << "Testing " << n << "/" << MAXN << "        \r" << flush;
      set_partition_generator sp(n);  int cnt = 0;
      vi items(n); iota(items.begin(), items.end(), 0); vi prev(n);
      do {
        vi partition = sp.next_partition();
        int hamming = 0;
        for (int i=0; i<n; i++) hamming += (prev[i] != partition[i]);
        assert(hamming == 1 || prev == vi(n,0));  prev = partition;
        sort(partition.begin(),partition.end());
        partition.erase(unique(partition.begin(),partition.end()),partition.end());
        for (int i=0;i<(int)partition.size();i++) assert(partition[i]==i);
        cnt++;
      } while (sp.has_next());
      assert(cnt == bell[n]);
    }
  }
}

int main() {
  brute_force_testing::test();
}
