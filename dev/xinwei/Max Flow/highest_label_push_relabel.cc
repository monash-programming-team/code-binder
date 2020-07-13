#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
const ll INF = LLONG_MAX;
struct edge {
    int a, b, cap, flow;
};

vector<vector<int> > g;
vector<edge> edges;
vector<ll> excess;
vector<int> h, inqueue, ct;
vector<vector<int> > bkt;
int maxBkt;
int n, s, t;

void add_edge(int u, int v, int cap){
    edge e1 = {u, v, cap, 0};
    edge e2 = {v, u, 0, 0};
    g[u].push_back( (int) edges.size() );
    edges.push_back(e1);
    g[v].push_back( (int) edges.size() );
    edges.push_back(e2);
}

void gap_heuristic(int k){
    for (int u = 0; u < n; u++){
	if (u == s || h[u] < k) continue;
	ct[h[u]]--;
	if (inqueue[u] && h[u] <= n) {
	    bkt[h[u]].clear();
	    bkt[n+1].push_back(u);
	}
	h[u] = max(h[u], n + 1);
	if (h[u] > maxBkt) maxBkt = h[u];
	ct[h[u]]++;
    }
}

void push(int u, int v, int id){
    ll tmp = min(excess[u], (ll) edges[id].cap - edges[id].flow); // int because it is limited by capacity
    excess[u] -= tmp;
    excess[v] += tmp;
    edges[id].flow += tmp;
    edges[id^1].flow -= tmp;
}

void discharge(int u){
    while (excess[u] > 0){
	int minH = 2 * n;
	for (int id : g[u]){
	    int v = edges[id].b;
	    if (edges[id].cap == edges[id].flow) continue;
	    if (h[u] == h[v] + 1){
		push(u, v, id);
		if (inqueue[v] == 0 && v != s && v != t){
		    bkt[h[v]].push_back(v);
		    if (h[v] > maxBkt) maxBkt = h[v];
		    inqueue[v] = 1;
		}
	    }
	    else 
		minH = min(minH, h[v]);
	}
	if (excess[u] > 0){
	    if (ct[h[u]] == 1){
		gap_heuristic(h[u]);
	    } else {
		ct[h[u]]--;
		h[u] = 1 + minH; // relabel
		ct[h[u]]++;
	    }
	}
    }
    inqueue[u] = 0;
}

ll push_relabel(int _n, int _s, int _t){
    n = _n; s = _s; t = _t;
    // initialize
    excess.assign(n, 0);
    h.assign(n, 0); h[s] = n;
    ct.assign(2 * n, 0); ct[0] = n - 1; ct[n] = 1;
    bkt.assign(2 * n, vector<int>());
    inqueue.assign(n, 0);
    maxBkt = 0;

    for (int id : g[s]){
	if (edges[id].cap == 0) continue;
	int u = edges[id].b;
	excess[u] += edges[id].cap;
	if (inqueue[u] == 0 && u != s && u != t){
	    bkt[0].push_back(u);
	    inqueue[u] = 1;
	}
	edges[id].flow += edges[id].cap;
	edges[id^1].flow -= edges[id].cap;
    }
    while (maxBkt >= 0){
	if (maxBkt < n && !bkt[maxBkt].empty()){
	    int u = bkt[maxBkt].back(); bkt[maxBkt].pop_back();
	    discharge(u);
	} else 
	    maxBkt--;
    }
    return excess[t];
}


 
int main(){
    int n,m,s,t;
    scanf("%d%d%d%d",&n,&m,&s,&t);
    g.assign(n, vector<int>());

    for(int i=0;i<m;i++){
	int u,v,cap;
	scanf("%d%d%d",&u,&v,&cap); u--; v--;
	if (u == v) continue;
	add_edge(u, v, cap);
    }


    ll mf = push_relabel(n, s-1, t-1);
    cout << mf << endl;
    return 0;
}
