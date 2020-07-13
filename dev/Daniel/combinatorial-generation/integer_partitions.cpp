// Integer partition generation
//
// Author: Daniel Anderson (based on code from indy256)
// Date: 20-01-2017
// Reliability: 1
// Tested on: Brute-force
//
#include<bits/stdc++.h>
using namespace std;

typedef long long int ll;

typedef vector<int> vi;
typedef vector<vi> vvi;

// Count the number of partitions of the integer n
ll count_partitions(int n) {
	vector<ll> p(n + 1);
	p[0] = 1;
	for (int i = 1; i <= n; i++)
		for (int j = i; j <= n; j++)
			p[j] += p[j - i];
	return p[n];
}
 
//listings:integer_partitions
// Generates the lexicographically next integer partition of sum(p). Start with
// p = [1,1,1,1,...] to generate all integer partitions in gray code order.
bool next_partition(vi& p) {
	int n = p.size(), i = n - 2;
	if (n <= 1) return false;
	int s = p.back() - 1; p.pop_back();
	while (i > 0 && p[i] == p[i - 1]) {	s += p[i--]; p.pop_back(); }
	p[i]++;
	while (s-- > 0) p.push_back(1);
	return true;
}
//listings:/integer_partitions

namespace brute_force_tests {
  int MAXN = 100;
  
  void test() {
    for (int n=1; n<=MAXN; n++) {
      cout << "Testing " << n << "/" << MAXN << "        \r" << flush;
      vi partition(n,1); int cnt = 0; vi prev;
      do {
        assert(partition > prev);  prev = partition;
        assert(accumulate(partition.begin(), partition.end(), 0) == n);
        cnt++;
      } while (next_partition(partition));
      assert(cnt == count_partitions(n));
    }
  }
}

int main() {
  brute_force_tests::test();
}