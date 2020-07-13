include <bits/stdc++.h>

using namespace std;
typedef long long ll;
typedef pair<int, int> pii;
typedef vector<int> vi;
typedef vector<vi> vvi;

template<typename T>
class rmq_update {
  typedef pair<T, int> pti;
private:
  int n;
  vector<T> v;
  vector<pti> st, lazy;

  pti merge(pti p1, pti p2) {
    if (p1.second == -1) return p2;
    if (p2.second == -1) return p1;
    return min(p1, p2);
  }
  void build(int p, int L, int R){
    if (L == R) st[p] = {v[L], L};
    else {
      int mid = (L + R) / 2;
      build(p * 2, L, mid);
      build(p * 2 + 1, mid + 1, R);
      st[p] = merge(st[p * 2], st[p * 2 + 1]);
    }
  }
  pti query(int p, int L, int R, int i, int j){
    if (i > R || j < L) return {-1, -1};
    if (i <= L && j >= R) return st[p];
    int mid = (L + R) / 2;
    pti p1 = query(p * 2, L, mid, i, j);
    pti p2 = query(p * 2 + 1, mid + 1, R, i, j);
    return merge(p1, p2);
  }
  void update(int p, int L, int R, int pos, T val){
    if (pos > R || pos < L) return;
    if (pos == L && pos == R) st[p].first = val;
    else {
      int mid = (L + R) / 2;
      update(p * 2, L, mid, pos, val);
      update(p * 2 + 1, mid + 1, R, pos, val);
      st[p] = merge(st[p * 2], st[p * 2 + 1]);
    }
  }

public:
  rmq_update(vector<T> &v) : v(v) {
    n = v.size();
    st.resize(4 * n);
    build(1, 0, n-1);
  }
  pti query(int i, int j){
    return query(1, 0, n-1, i, j);
  }
  void update(int pos, T val){
    update(1, 0, n-1, pos, val);
  }
};

#define ff first
#define ss second

#ifndef ONLINE_JUDGE
#define dbg(x) cerr << __LINE__ << ": " << #x << " = " << (x) << ' ';
#else
#define dbg(x)
#endif

const int INF = int(1e9);
int count(const vi &v){
    int ct = 0;
    for (int x : v) ct += (x > 0);
    return ct;
}

int diff(const vi &v1, const vi &v2){
    vi v(26);
    for (int i = 0; i < 26; i++)
        v[i] = v1[i] - v2[i];
    return count(v);
}

int main(){
    std::ios_base::sync_with_stdio(false);
    int k, n;
    cin >> k;
    string str; cin >> str;
    n = (int) str.size();
    vi res(n, INF);
    vi ans(n, -1);
    rmq_update<int> cc(res);
    vvi LT(n);
    vi letters(26, 0);
    for (int i = 0; i < n; i++){
        letters[str[i] - 'a']++;
        LT[i] = letters;
    }
    int L = 0, R = 0;
    for (int i = 0; i < n; i++){
        int ct = count(LT[i]);
        if (ct < k) continue;
        if (ct == k) ;
        else {
            while (diff(LT[i], LT[L]) > k) L++;
        }
        while (diff(LT[i], LT[R]) >= k) R++;

        if (ct == k){
            ans[i] = 1;
            cc.update(i, 1);
        } else {
            int Mi = cc.query(L, R - 1).first;
            if (Mi < INF){
                ans[i] = Mi + 1;
                cc.update(i, Mi + 1);
            }
        }
    }
    for (int x : ans)
        cout << x << ' ';
    cout << endl;

    return 0;
}
