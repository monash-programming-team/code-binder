#include <bits/stdc++.h>

#define x first
#define y second

using namespace std;

typedef long long ll;
typedef pair<ll,ll> pll;
typedef pair<int,int> ii;

ll operator* (pll a,pll b) {
	return a.x*b.x+a.y*b.y;
}

ll abs2(pll a) {return a*a;}

pll operator-(pll a) {
	return {-a.x,-a.y};
}

pll operator+(pll a,pll b) {
	return {a.x+b.x,a.y+b.y};
}

pll operator-(pll a,pll b) {
	return a+-b;
}

ll sq(ll a) {return a*a;}

vector<pll> byy,byyl,byyr;
pll P1,P2;
ll m;

void cpop(int,int);

void cpop(const vector<pll>& pts) {
	byy=pts;
	sort(byy.begin(),byy.end());
	m=LLONG_MAX;
	cpop(0,(int)byy.size());
}

void cpop(int i,int j) {
	if (j-i<2) return;
	int mid=(i+j)/2;
	ll div=byy[mid].x;
	cpop(i,mid);
	cpop(mid,j);

	byyl.clear();
	byyr.clear();
	for (int c=i;c<j;c++) if (sq(div-byy[c].x)<=m) {
		if (c<mid) byyl.push_back(byy[c]);
		else byyr.push_back(byy[c]);
	}
	
	for (int c=0,l=0,u=0;c<byyl.size();c++) {
		while (u<byyr.size() && (byyr[u].y<byyl[c].y || sq(byyr[u].y-byyl[c].y)<=m)) u++;
		while (l<u && sq(byyr[l].y-byyl[c].y)>=m) l++;
		assert(u-l<=6);
		ll cur;
		for (int d=l;d<u;d++) if ((cur = abs2(byyr[d]-byyl[c]))<m) {
			m=cur;
			P1=byyl[c];
			P2=byyr[d];
		}
	}

	inplace_merge(byy.begin()+i,byy.begin()+mid,byy.begin()+j,
			[](const pll a,const pll b) {return a.y<b.y || (a.y==b.y && a.x<b.x);});
}

int main() {
	int N;
	while (cin >> N && N) {
		byy.resize(N);
		for (int c=0;c<N;c++) scanf("%I64d%I64d",&byy[c].x,&byy[c].y);
		cpop(byy);
		cout << fixed << setprecision(10);
		/*if (m<=100000000ll)*/ cout << sqrt(m) << endl;
		//else cout << "INFINITY\n";//
	}
}
