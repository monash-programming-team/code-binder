#include <bits/stdc++.h>

using namespace std;

typedef vector<int> vi;

// Suffix Array + LCP
//
// Author      : Peter Whalan
// Date        : July 16, 2016
// Reliability : 2
// Tested On   : tests below and Monash Group Problem
//
// SA Sorts the suffixes of a string lexicographically. LCP determines the
// longest common prefix between adjacent suffixes.
//
// The offset to the start of the suffix in the string is used to refer to that
// suffix. For example the whole string is suffix 0. sa[0] is always the empty
// suffix.
//
// SA takes a vi. A string, s, can be cast to a vi using vi(s.begin(),s.end).
// The vi must have only positive values.
// 
// lcp[i] is the longest common prefix of sa[i] and sa[i-1]. lcp[0]==0.
//
// Complexity:
//		O( N log N ) - SA
//		O( N ) - LCP

// If s.size() is 1 greater than a power of 2 then k (in bs) can be one less
// than a power of 2. For this value of k empty suffix r[2*s.size()-1] will be
// accessed. Furthermore, for any s.size(), r[s.size()+1] will also be accessed.
// This is greater for s.size()<2. This is why there is a +2 in the
// initialisation of r and tr.

//max(max length of string, max alphabet value)
const int maxr = 10000000;//max rank
vi r/*(2*maxr+2)*/,tr/*(2*maxr+2)*/,sa,tsa,ct;//rank, temp rank, suffix array, temp suffix array, count
int N,Mr;//size of suffix array (includes empty suffix), current max rank

void bs(int k) {//bucket sort
	ct.assign(maxr+2,0);
	for (int c=0;c<N;c++) ct[r[c+k]+1]++;
	for (int c=1;c<Mr;c++) ct[c+1]+=ct[c];
	for (int c=0;c<N;c++) tsa[ct[r[sa[c]+k]]++]=sa[c];//Do in suffix array order
	swap(sa,tsa);                                     //so the sort is stable.
}


void SA(const vi& s) {
	N=s.size()+1;
	sa.resize(N);tsa.resize(N);
	iota(sa.begin(),sa.end(),0);

	//Modify this to change the rank of characters.
	for (int c=0;c<N-1;c++) r[c]=s[c];//char cast to int for char rank
	fill(r.begin()+N-1,r.end(),0);//0 used as rank of terminating char

	Mr=maxr;//This can be changed to max in s
	for (int k=1;k==1 || Mr+1!=N;k<<=1) {
		bs(k);
		bs(0);
		Mr=-1;
		for (int c=0;c<N;c++) tr[sa[c]]=
			c && r[sa[c]]==r[sa[c-1]] && r[sa[c]+k]==r[sa[c-1]+k]?Mr:++Mr;
		swap(r,tr);
	}
}

vi LCP(const vi& s) {
	vi lcp(N);

	//Uncomment this if r wasn't just computed by running SA.
	//for (int c=0;c<N;c++) r[sa[c]]=c;
	
	for (int i=0,w=0;i<N;i++) {
		int k=r[i];
		if (k) {
			for (int j=sa[k-1];max(i,j)+w+1<N && s[i+w]==s[j+w];w++);
			lcp[k]=w;
		}
		if (w) w--;
	}

	return lcp;
}



int main() {
	string s;

	cin >> s;
	SA(vi(s.begin(),s.end()));

	for (int c=0;c<N;c++) {
		if (c!=1) printf(" ");
		printf("%d",sa[c]);
	}
	cout << endl;

	vi lcp=LCP(vi(s.begin(),s.end()));
	for (int c=0;c<N;c++) {
		if (c!=0) printf(" ");
		printf("%d",lcp[c]);
	}
	cout << endl;
	
	
	/*for (int len=0;len<=10000;len++) {
		maxr=len;
		r.assign(2*maxr+2,0);tr.assign(2*maxr+2,0);
		SA(vi(s.begin(),s.end()));

		for (int i=0;i<N;i++) assert(sa[i]==N-1-i);

		vi lcp=LCP(vi(s.begin(),s.end()));

		for (int i=0;i<N;i++) assert((i<2 && lcp[i]==0) || (i>=2 && lcp[i]==i-1));

		s.push_back(1);
	}*/
	
}
