#include <bits/stdc++.h>

using namespace std;

typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;

struct Edge {
	int j,C;
};

int n,s,t;
vector<Edge> el;
vvi al;

vector<bool> vis;

int dfs(int i,int cf) {
	if (vis[i]) return 0;
	vis[i]=1;
	if (i==t) return cf;
	for (int c:al[i]) {
		int ncf = min(cf,el[c].C),f;
		if (ncf && (f=dfs(el[c].j,ncf))) {
			el[c].C-=f;
			el[c^1].C+=f;
			return f;
		}
	}

	return 0;
}
int x;
ll maxflow() {
	ll Mf=0,f=1;
	while (f && Mf<x) {
		vis.assign(n,0);
		f=dfs(s,INT_MAX);
		Mf+=f;
	}
	return Mf;
}

int oc[1000];

int main() {
	int m;
	scanf("%d%d%d",&n,&m,&x);
	s=0;t=n-1;
	al.assign(n,{});
	el.resize(2*m);
	for(int i=0;i<m;i++){
		int u,v,c;
		scanf("%d%d%d",&u,&v,&c);
		u--;v--;
		al[u].push_back(2*i);
		al[v].push_back(2*i+1);
		oc[2*i]=c;
		oc[2*i+1]=0;
		el[2*i]={v,c};
		el[2*i+1]={u,0};//change 0 to c for undirected egde
	}

	double l=0,u=5e10;
	double mid=(l+u)/2;
	//while (u-l>1e-9) {
	while ((mid-l)*x/max(1.0,l*x)>1e-7) {
	//for (int c=0;c<60;c++) {
		//double mid=(l+u)/2;
		for (int i=0;i<2*m;i+=2) el[i].C=min((double)INT_MAX,oc[i]/mid);
		if (maxflow() >= x) l=mid;
		else u=mid;
		mid=(l+u)/2;
	}

	cout << fixed << setprecision(20) << mid*x << endl;
}
