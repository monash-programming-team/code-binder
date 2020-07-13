#include <bits/stdc++.h>

using namespace std;

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef long long ll;

int s,t;
vector<ll> OF;vi H; //overflow,height

int n,m;
vector<Edge> es;
vvi al;

struct Edge {
	int v[2];//vertices connected
	int f[2];//f[0] capacity from v[0] to v[1]

	bool send(int i) {//returns capacity hit or height too low
		int j=i==v[1];
		if (H[i]<=H[v[!j]]) return 1;
		int a=min((ll)f[j],(i!=s)*OF[i]);
		f[j]-=a;
		f[!j]+=a;
		OF[v[!j]]+=a;
		if (i) OF[i]-=a;
		return !f[j];
	}
};


ll maxflow() {
	OF.assign(n,0);
	H.assign(n,0);

	for (int c:al[s]) es[c].send(s);

	for (int done=0,i=0;done<n-2

int main() {
	scanf("%d%d%d%d",&n,&m,&s,&t);s--;t--;
	es.resize(m);
	al.resize(n);
	int j=0;
	for(int i=0;i<m;i++){
  		int u,v,c;
		scanf("%d%d%d",&u,&v,&c);
		
		al[u-1].push_back(j);
		al[v-1].push_back(j);
		es[j++]=Edge(u-1,v-1,c,c);
	}
