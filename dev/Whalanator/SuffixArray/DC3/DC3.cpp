#include <bits/stdc++.h>

using namespace std;

typedef vector<int> vi;

struct SA : vi {
	vi s;

	SA(string s) s(s.begin(),s.end()) {
		sa(0);
	}

	void sa(int i) {
		int N=s.size()-i;
		int A=i?256:N+1;
		s.insert(s.end(),vi(3));
		//for (int c=0;c<3;c++) s.push_back(-1);
		vi dc,b,tem;
		for (int c=1;c<=N;c++) if (c%3) dc.push_back(i+c);
		tem.resize(dc.size());
		for (int o=2;o>=0;o--) {
			b.assign(A,0);
			for (int c=0;c<dc.size();c++) b[s[dc[c]+o]]++;
			for (int c=0,sum=0;c<A;c++) {
				int t=b[c];b[c]=sum;sum+=t;
			}

			for (int c=0;c<dc.size();c++) tem[b[s[dc[c]+o]]++]=dc[c];
			dc=tem;
		}
	}

	bool eq(int i,int j) {
		int o=0;
		for (;o<3 && s[i+o]==s[j+o];o++);
		return o==3;
	}
