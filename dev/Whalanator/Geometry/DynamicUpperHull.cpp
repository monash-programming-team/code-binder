#include <bits/stdc++.h>

#define x first
#define y second

using namespace std;

typedef pair<int,int> ii;
typedef vector<ii> vii;

// Dynamic Upper Hull
//
// Author      : Peter Whalan
// Date        : June 22, 2016
// Reliability : 
// Tested On   : 
//
// Performs range sum queries on array. Allows individual values to be
// modified.
// 
// Treats arguments as 0-based
//
// Complexity:
//		O( N ) to build
//		O( log N ) to update and query


struct UpperEnvelope {
	set<ii> U;

	bool midrel(ii a,ii b,ii c) {//Watch out for overflow
		return (a.y-b.y)*(c.x-b.x) < (b.y-c.y)*(b.x-a.x);
		//Change < to <= to include lines with single point on envelope
		//and identical lines
	}

	void add(ii line) {
		if (U.empty()) {U.insert(line);return;}
		auto r=U.lower_bound(line),l=r;
		if (r==U.begin()) {
			if (line.x==r->x && line.y<=r->y) return;
		}
		else {
			l--;
			if (r!=U.end() && !midrel(*l,line,*r)) return;
		}

		while (l!=U.begin()) {
			auto ol=l--;
			if (midrel(*l,*ol,line)) break;
			U.erase(ol);
		}
		if (r!=U.end()) {
			auto olr=r++;
			while (r!=U.end()) {
				if (midrel(line,*olr,*r)) break;
				U.erase(olr);
				olr=r++;
			}
		}
		
		U.insert(line);
		printf("lsize %d\n",U.size());
	}

	int envval(int x) {
		if (U.empty()) return INT_MIN;
		int l=U.begin()->x,u=U.rbegin()->x+1;//half open interval
		set<ii>::iterator it;
		while (l+1!=u) {
			int m=(l+u)/2;//This m isn't necessarily the gradient of a line in U
			it=U.lower_bound({m,INT_MIN});
			if (it==U.begin() || prev(it)->y-it->y<=x*(it->x-prev(it)->x)) l=m;
			else u=m;
		}
		it=U.lower_bound({l,INT_MIN});
		return it->x*x+it->y;
	}
};

int main() {
	UpperEnvelope ue;
	while (1) {
		string s;
		if (!getline(cin,s)) break;
		stringstream ss(s);
		int i,j;
		if (ss >> i >> j) ue.add({i,j});
		else cout << ue.envval(i) << endl;
	}
}
