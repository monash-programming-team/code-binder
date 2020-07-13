#include <bits/stdc++.h>

using namespace std;

typedef vector<int> vi;

// Longest Palindrome Substring
//
// Author      : Peter Whalan
// Date        : March 1, 2016
// Reliability : 0
// Tested On   : 
//
// Longest palindrome substring centered at every point
// returns vector of "radii" (such that 2*r+1 = width) after padding with spaces (including edges)
// This happens to correspond to the width.
// Pass same string to both params for classical problem
// Second string is what you want the string to match after reversing.
//
// Complexity: O( N )

void is(string& s) {//insert spaces
	int N=s.size();
	s.resize(2*N+1,' ');
	for (int c=N-1;c>=0;c--) {s[2*c+1]=s[c];s[2*c]=' ';}
}

vi lps(string s,string t) {
	int N=s.size();
	is(s);is(t);//Comment this out if only want odd length palindromes
	//and interpret the answer differently.
	N=2*N+1;
	vi r(N);
	
	int Mi=0;
	for (int c=1;c<N;c++) {
		int Mw=r[Mi]+Mi-c;
		int w=r[c]=(Mw>=0?min(Mw,r[2*Mi-c]):0);
		if (w>=Mw) {
			int i=c-w,j=c+w;
			while (++j<N && --i>=0 && s[i]==t[j] && s[j]==t[i]);//Only place where t used.
			r[Mi=c]=j-1-c;
		}
	}

	return r;
}

int main() {
	int N;
	string s;
	cin >> N >> s;
	vi r=lps(s,s);
	cout << *max_element(r.begin(),r.end()) << endl;
	return 0;
}

