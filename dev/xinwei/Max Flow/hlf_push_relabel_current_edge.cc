// Tested on: CF - 653D
//            Tuesday Tutorial: Network Flow
// 	          UVA - 820 (Internet Bandwidth)
//            CF - 269C
// 	          SPOJ - FASTFLOW
#include <bits/stdc++.h>

using namespace std;
typedef long long ll;
typedef vector<int> vi;
typedef vector<vi> vvi;
const ll INF = LLONG_MAX;
struct edge {
    int a, b, cap, flow;
};

vvi g, bucket;
vector<edge> edges;
vector<ll> excess;
vi h, inqueue, num_height, current_edge;
int max_bucket, n, s, t;

void add_edge(int u, int v, int cap){
    if (u == v) return;
    edge e1 = {u, v, cap, 0};
    edge e2 = {v, u, 0, 0};
    g[u].push_back( (int) edges.size() );
    edges.push_back(e1);
    g[v].push_back( (int) edges.size() );
    edges.push_back(e2);
}

void gap_heuristic(int k){
    for (int u = 0; u < n; u++){
	if (u == s) continue; 
	if (h[u] > k && h[u] <= n){  
	    num_height[h[u]]--;
	    if (inqueue[u]){
		bucket[h[u]].clear();
		bucket[n+1].push_back(u);
	    } 
	    h[u] = n + 1;
	    num_height[h[u]]++;
	    if (h[u] > max_bucket) max_bucket = h[u];
	    current_edge[u] = 0;
	}
    }
}

void push(int u, int v, int id){
    ll tmp = min(excess[u], (ll) edges[id].cap - edges[id].flow); 
    excess[u] -= tmp;
    excess[v] += tmp;
    edges[id].flow += tmp;
    edges[id^1].flow -= tmp;
}

int relabel(int u){
    int minH = 2 * n;
    for (int id : g[u]){
	if (edges[id].cap != edges[id].flow){
	    minH = min(minH, h[edges[id].b]);
	}
    }
    return 1 + minH;
}

void discharge(int u){
    inqueue[u] = 0;
    while (excess[u] > 0){
	for (; current_edge[u] < (int) g[u].size(); current_edge[u]++){
	    int id = g[u][current_edge[u]];
	    int v = edges[id].b;
	    if (edges[id].cap == edges[id].flow) continue;
	    if (h[u] == h[v] + 1){
		push(u, v, id);
		if (inqueue[v] == 0 && v != s && v != t){
		    bucket[h[v]].push_back(v);
		    if (h[v] > max_bucket) max_bucket = h[v];
		    inqueue[v] = 1;
		}
	    } 
	    if (excess[u] == 0) break; // need to remain at current_edge since still admissable
	}

	if (excess[u] > 0){
	    int prev_h = h[u];
	    // relabel
	    num_height[h[u]]--;
	    h[u] = relabel(u);
	    num_height[h[u]]++;
	    current_edge[u] = 0;

	    // apply gap heuristic
	    if (num_height[prev_h] == 0 && prev_h <= n - 1){
		gap_heuristic(prev_h);
	    } 
	} 
    }
}

ll push_relabel(int _n, int _s, int _t){
    n = _n; s = _s; t = _t;
    // initialize
    excess.assign(n, 0);
    current_edge.assign(n, 0);
    h.assign(n, 0); h[s] = n;
    num_height.assign(2 * n, 0); num_height[0] = n - 1; num_height[n] = 1;
    bucket.assign(2 * n, vector<int>());
    inqueue.assign(n, 0);
    max_bucket = 0;

    for (int id : g[s]){
	if (edges[id].cap == 0) continue;
	int u = edges[id].b;
	excess[u] += edges[id].cap;
	if (inqueue[u] == 0 && u != s && u != t){
	    bucket[0].push_back(u);
	    inqueue[u] = 1;
	}
	edges[id].flow += edges[id].cap;
	edges[id^1].flow -= edges[id].cap;
    }
    while (max_bucket >= 0){
	if (!bucket[max_bucket].empty()){ // optional (max_bucket < n)
	    int u = bucket[max_bucket].back(); bucket[max_bucket].pop_back();
	    discharge(u);
	} else {
	    max_bucket--;
	}
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

	add_edge(u, v, cap);
    }

    ll mf = push_relabel(n, s-1, t-1);
    cout << mf << endl;
    return 0;
}
