// Bitmask subset generation
//
// Author: Daniel Anderson
// Date: 20-01-2017
// Reliability: 1
// Tested on: Brute-force
//
#include<bits/stdc++.h>
using namespace std;

typedef long long int ll;

typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:subsets_of_a_subset
// Generate the lexicographically previous subset of mask. To generate all subsets,
// the initial subset should be sub = mask.
template<typename T> bool next_subset(T& sub, T mask) {
  if (sub == 0) return false; else return sub = (sub - 1) & mask, true;
}
//listings:/subsets_of_a_subset

//listings:combinations
// Generate the lexicographically next combination of n choose k. To generate all
// combinations, the initial combination is comb = (1 << k) - 1.
template<typename T> bool next_combination(T& comb, int n, int k) {
  T x = comb & -comb, y = comb + x;  comb = (((comb ^ y) >> 2) / x) | y;
  return comb < (T(1) << n);
}
//listings:/combinations

ll nCk(int n, int k) {
  ll ans = 1;
  for (int i=n-k+1; i<=n; i++) ans *= i;
  for (int i=1; i<=k; i++) ans /= i;
  return ans;
}  

// Generate all combinations of n choose k and check that
// 1. The correct number of subsets was generated
// 2. Each contained the correct number of elements
// 3. None were generated multiple times
// 4. All contained bits lower than position n
void test_n_choose_k(int n, int k) {
  ll comb = (1LL << k) - 1;
  set<ll> seen;
  int cnt = 0;
  do {
    cnt++;
    assert(__builtin_popcountll(comb) == k);
    assert(seen.insert(comb).second);
    assert(comb < (1LL << n));
  } while (next_combination(comb, n, k));
  assert(cnt == nCk(n,k));
}

// Generate all subsets of mask and check that
// 1. The correct number of subsets was generated
// 2. Each was generated once only
// 3. Each was a subset of the given mask
void test_all_subsets(ll mask) {
  ll sub = mask; int cnt = 0;
  set<ll> seen;
  do {
    cnt++;
    assert(seen.insert(sub).second);
    assert((sub & mask) == sub);
  } while (next_subset(sub, mask));
  assert(cnt == (1LL << __builtin_popcountll(mask)));
}

int main() {
  test_all_subsets((1 << 20) - 1);
  test_n_choose_k(60, 5);
}
