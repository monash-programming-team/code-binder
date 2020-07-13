/******************************************************************************
 *								Combinatorial Algorithms
 * ---------------------------------------------------------------------------
 * Algorithms for generating, listings and enumerating common combinatorial
 * patterns.
 * ---------------------------------------------------------------------------
 * References:
 * (1) The Art of Computer Programming, Volume 4A : Combinatorial Algorithms
 * 		- Donald Knuth
 * (2) https://github.com/indy256/codelibrary
 *****************************************************************************/
#include<bits/stdc++.h>
using namespace std;

typedef long long int ll;

typedef vector<int> vi;
typedef vector<vi> vvi;

/******************************************************************************
 *								Generate Arrangements
 * ---------------------------------------------------------------------------
 * Generates all of the vectors of length m = a.size() consisting of unique
 * elements from [0, n). The initial arrangement must be {a[i] = i}. The final
 * arrangement is {a[i] = n - 1 - i}.
 * ---------------------------------------------------------------------------
 *****************************************************************************/
bool next_arrangement(vi& a, int n) {
	vector<bool> used(n);
	for (int x : a) used[x] = true;
	int m = a.size();
	for (int i = m - 1; i >= 0; i--) {
		used[a[i]] = false;
		for (int j = a[i] + 1; j < n; j++) {
			if (!used[j]) {
				a[i++] = j;
				used[j] = true;
				for (int k = 0; i < m; k++) 
					if (!used[k]) a[i++] = k;
				return true;
			}
		}
	}
	return false;
}
 
/******************************************************************************
 *						Generate Arrangements with Repeats
 * ---------------------------------------------------------------------------
 * Generates all of the vectors of length m = a.size() consisting of elements 
 * from [0, n). The initial arrangement must be {a[i] = 0}. The final
 * arrangement is {a[i] = n - 1}.
 * ---------------------------------------------------------------------------
 *****************************************************************************/
 bool next_arrangement_with_repeats(vi& a, int n) {
	 int m = a.size();
	 for (int i = m - 1; i >= 0; i--) {
		 if (a[i] < n - 1) {
			++a[i];
			fill(a.begin() + i + 1, a.end(), 0);
			return true;
		 }
	 }
	 return false;
 }
 
 /******************************************************************************
 *							Generate Subsets
 * ---------------------------------------------------------------------------
 * Generates all of the subsets of a set.
 * ---------------------------------------------------------------------------
 *****************************************************************************/
 template<typename T>
 struct subset_enumerator {
	int s, n, m;
	vector<T> a;
	subset_enumerator(vector<T> a) : a(a) {
		s = 0;
		n = a.size();
		m = 1 << a.size();
	}
	bool has_next_subset() { return s < m; }
	vector<T> next_subset() {
		vector<T> subset;
		for (int i = 0; i < n; i++)
			if (s & (1 << i)) subset.push_back(a[i]);
		s++;
		return subset;
	}
 };
 
  /******************************************************************************
 *							Generate Combinations
 * ---------------------------------------------------------------------------
 * Generates all combinations of size k = p.size() elements from a set of
 * n elements.
 * ---------------------------------------------------------------------------
 *****************************************************************************/
// Enumerates combinations with no repeats. The initial combination must be
// p[i] = i.
bool next_combination(vi& p, int n) {
	int m = p.size();
	for (int i = m - 1; i >= 0; i--) {
		if (p[i] < n - m + i) {
			++p[i];
			while (++i < m) p[i] = p[i - 1] + 1;
			return true;
		}
	}
	return false;
}

// Enumerates combinations with repeats. The initial combination must be p[i] = 0.
bool next_combination_with_repeats(vi& p, int n) {
	int m = p.size();
	for (int i = m - 1; i >= 0; i--) {
		if (p[i] < n - 1) {
			++p[i];
			while (++i < m) {
				p[i] = p[i - 1];
			}
			return true;
		}
	}
	return false;
}
 
/*
//listings:subsets_of_a_subset
// Generate all non-empty subsets of a given subset (mask)
// Use ll or __int128 instead of int if you have > 32 items.
for (int sub = mask; sub > 0; sub = (sub - 1) & mask) { ... }
//listings:/subsets_of_a_subset
*/

/*
//listings:combinations
// Generate all subsets of size k from a set of size n
// Use ll or __int128 instead of int if you have > 32 items.
for (int comb = (1 << k) - 1; comb < (1 << n); ) {
  ...
  int x = comb & -comb, y = comb + x;
  comb = ((comb & -y) / x >> 1) | y;
}
//listings:/combinations
*/
 
/******************************************************************************
 *							Set Partitions
 * ---------------------------------------------------------------------------
 * Generates/enumerates all of the set partitions of a set of size n.
 * ---------------------------------------------------------------------------
 *****************************************************************************/
 
 //listings:set_partitions
// Generates set partitions for a set of size n in gray code order. A set partition
// is represented as a vector of size n where P[i] = the index of the set that
// element i belongs to. Safe to use for n <= 13, where B_n = 27 million.
struct set_partition_generator {
	int n;  vi a, b, d;  bool done = false;
	set_partition_generator(int n) : n(n), a(n), b(n, 1), d(n, 1) { }
	void fix(int j, int m) { fill(b.begin() + j + 1, b.end(), m); }
	bool has_next_partition() { return !done; }
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
 
// Generates all set partitions at once
vvi generate_set_partitions(int n) {
	if (n == 1) return {{0}};
	vvi p = generate_set_partitions(n - 1);
	int l = p.size();
	vi m(l);
	int len = 0;
	for (int i = 0; i < l; i++) {
		for (int v : p[i])
			m[i] = max(m[i], 1 + v);
		len += m[i] + 1;
	}
	vvi res(len, vi(n));
	for (int i = 0, pos = 0; i < l; i++) {
		for (int j = 0; j <= m[i]; j++, pos++) {
			copy(p[i].begin(), p[i].end(), res[pos].begin());
			res[pos][n - 1] = (i % 2 == 1 ? j + 1 : m[i] - j + 1) % (m[i] + 1);
		}
	}
	return res;
}	
 
/******************************************************************************
 *							Generate Nested Parentheses
 * ---------------------------------------------------------------------------
 * Generates all strings of matched nested parentheses
 * ---------------------------------------------------------------------------
 * - There are C_n sets of n matched parentheses (C_n = n'th Catalan number)
 *****************************************************************************/

//listings:parenthesis
// Generates all strings of n pairs (2n chars) of balanced, nested parentheses.
struct parentheses_enumerator {
	int n, m;
	bool done = false;
	string a;
	parentheses_enumerator(int n) : n(n), m(2 * n - 1) {
		for (int i = 0; i < m + 2; i++) a.push_back(")("[i % 2]);
	}
	bool has_next_string() { return !done; }
	string next_string() {
		string ans = string(a.begin()+1, a.end());
		a[m] = ')';
		if (a[m - 1] == ')') a[--m] = '(';
		else {
			int j = m - 1, k = 2 * n - 1;
			while (a[j] == '(') {	a[j--] = ')';	a[k] = '(';	k -= 2;	}
			if (j == 0) done = true;
			a[j] = '(',	m = 2 * n - 1;
		}
		return ans;
	}
};
//listings:/parenthesis
 
/******************************************************************************
 *							Generate Integer Partitions
 * ---------------------------------------------------------------------------
 * Generates all of the partitions of an integer.
 * ---------------------------------------------------------------------------
 *****************************************************************************/
 
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

// Aux function for generate_increasing_partitions(int n)
void generate_increasing_partitions(vi& p, int left, int last, int pos, vvi& partitions) {
	if (left == 0) {
		vi partition(pos);
		for (int i = 0; i < pos; i++)
			partition[i] = p[i];
		partitions.push_back(partition);
		return;
	}
	for (p[pos] = last + 1; p[pos] <= left; p[pos]++)
		generate_increasing_partitions(p, left - p[pos], p[pos], pos + 1, partitions);
}

// Generate all integer partitions of n that are strictly increasing sequences
vvi generate_increasing_partitions(int n) {
	vvi partitions;
	vi p1(n);
	generate_increasing_partitions(p1, n, 0, 0, partitions);
	return partitions;
}
 
/******************************************************************************
 *									Enumerators
 *****************************************************************************/
 
// Applies the given function to all subsets of the given vector
template<typename T, typename F>
void foreach_subset(vector<T> a, F f) {
	int n = a.size();
	for (int s = 0; s < (1 << n); s++) {
		vector<T> subset;
		for (int i = 0; i < n; i++)
			if (s & (1 << i)) subset.push_back(a[i]);
		f(subset);
	}
}

// Applies the given function to all permutations of the given vector
template<typename T, typename F>
void foreach_permutation(vector<T> a, F f) {
	vector<T> b = a;
	sort(b.begin(), b.end());
	do { f(b); } while (next_permutation(b.begin(), b.end()));
}

// Applies the given function to all permutations of all subsets of the given vector
template<typename T, typename F>
void foreach_permutation_subset(vector<T> a, F f) {
	foreach_subset(a, [&](vector<T>& s) {
		foreach_permutation(s, f);
	});
}

// Applies the given function to all subsets of size k of the given vector a
template<typename T, typename F>
void foreach_combination(vector<T> a, int k, F f) {
	int n = a.size();
	vi comb(k);
	iota(comb.begin(), comb.end(), 0);
	do {
		vi subset; subset.reserve(k);
		for (int i : comb) subset.push_back(a[i]);
		f(subset);
	} while (next_combination(comb, n));
}

// Applies the given function to all combinations (with repeats) of k elements of
// the given vector a
template<typename T, typename F>
void foreach_combination_with_repeats(vector<T> a, int k, F f) {
	int n = a.size();
	vi comb(k);
	do {
		vi subset; subset.reserve(k);
		for (int i : comb) subset.push_back(a[i]);
		f(subset);
	} while (next_combination_with_repeats(comb, n));
}

// Applies the given function to all arrangements of k elements of the
// given vector a
template<typename T, typename F>
void foreach_arrangement(vector<T> a, int k, F f) {
	int n = a.size();
	vi arr(k);
	iota(arr.begin(), arr.end(), 0);
	do {
		vi arrangement; arrangement.reserve(k);
		for (int i : arr) arrangement.push_back(a[i]);
		f(arrangement);
	} while (next_arrangement(arr, n));
}

// Applies the given function to all arrangements with repeats of k
// elements of the given vector a
template<typename T, typename F>
void foreach_arrangement_with_repeats(vector<T> a, int k, F f) {
	int n = a.size();
	vi arr(k);
	iota(arr.begin(), arr.end(), 0);
	do {
		vi arrangement; arrangement.reserve(k);
		for (int i : arr) arrangement.push_back(a[i]);
		f(arrangement);
	} while (next_arrangement_with_repeats(arr, n));
}

// Applies the given function to all set partitions of the given vector
template<typename T, typename F>
void foreach_set_partition(vector<T> a, F f) {
	int n = a.size();
	set_partition_generator sp(n);
	while (sp.has_next_partition()) {
		vi part = sp.next_partition();
		vector<vector<T>> partition(*max_element(part.begin(), part.end()) + 1);
		for (int i = 0; i < n; i++)
			partition[part[i]].push_back(a[i]);
		f(partition);
	}
}
 
 /******************************************
 *					Examples
 ******************************************/
 
void stuff(vector<vector<char>>& partition) {
	cout << partition.size() << ": ";
	for (auto& part : partition) cout << part.size() << ' ';
	cout << endl;
}
 
int main() {
	cout << "Arrangements of 0 to 3:" << endl;
	vi a = {0,1,2};
	do {
		for (int x : a) cout << x << ' '; cout << endl;
	} while (next_arrangement(a, 4));
	
	cout << "Arrangements of 0 to 2 with repeats:" << endl;
	a = {0,0,0};
	do {
		for (int x : a) cout << x << ' '; cout << endl;
	} while (next_arrangement_with_repeats(a, 3));
	
	cout << "Subsets:" << endl;
	a = {0,1,2,3,4,5};
	subset_enumerator<int> sub(a);
	while (sub.has_next_subset()) {
		auto s = sub.next_subset();
		for (int x : s) cout << x << ' '; cout << endl;
	}
	
	int outsider = 10;
	
	cout << "Sum subsets:" << endl;
	foreach_subset(a, [&](vi& s) {
		int sum = outsider;
		for (int x : s) sum += x;
		cout << sum << endl;
	});
	
	cout << "Permutations of subsets:" << endl;
	a = {0,1,2};
	foreach_permutation_subset(a, [](vi& p) {
		for (int x : p) cout << x << ' '; cout << endl;
	});
	
	cout << "Combinations:" << endl;
	a = {0,1,2,3};
	do {
		for (int x : a) cout << x << ' '; cout << endl;
	} while (next_combination(a, 5));
	
	cout << "Combinations with repeats:" << endl;
	a = {0,0,0};
	do {
		for (int x : a) cout << x << ' '; cout << endl;
	} while (next_combination_with_repeats(a, 3));
	
	cout << "Integer partitions:" << endl;
	vi p = {1,1,1,1,1,1};
	do {
		for (int x : p) cout << x << ' '; cout << endl;
	} while (next_partition(p));
	
	cout << "Increasing partitions:" << endl;
	vvi partitions = generate_increasing_partitions(8);
	for (auto& p : partitions) {
		for (int x : p) cout << x << ' '; cout << endl;
	}
	
	parentheses_enumerator pe(4);
	cout << "Nested parentheses: " << endl;
	while (pe.has_next_string()) {
		string p2 = pe.next_string();
		cout << p2 << endl;
	}
	
	cout << "Set partitions:" << endl;
	vvi parts = generate_set_partitions(4);
	set_partition_generator sp(4);
	for (auto& part : parts) {
		assert(sp.has_next_partition());
		vi part2 = sp.next_partition();
		assert(part == part2);
		for (auto x : part) cout << x << ' '; cout << endl;
	}
	assert(!sp.has_next_partition());
	
	cout << "Set Partitions:" << endl;
	vector<char> c = {'a', 'b', 'c', 'd'};
	foreach_set_partition(c, [](vector<vector<char>>& partition) {
		for(auto part : partition) {
			cout << '{';
			for (auto p : part) cout << p << ' ';
			cout << "}, ";
		}
		cout << endl << endl;
	});
	
	foreach_set_partition(c, stuff);
	
}