#include <bits/stdc++.h>

using namespace std;
typedef vector<int> vi;

vi prefix(const string &pat){
    int m = (int) pat.size(), j = -1;
    vi P(m + 1, -1);
    for (int i = 0; i < m; i++){
        while (j >= 0 && pat[i] != pat[j]) j = P[j];
        j++;
        P[i + 1] = j;
    }
    return P;
}

vi match(const string &txt, const string &pat, const vi &P){
    int n = (int) txt.size(), m = (int) pat.size(), j = 0;
    vi idx;
    for (int i = 0; i < n; i++){
        while (j >= 0 && txt[i] != pat[j]) j = P[j];
        j++;
        if (j == m){
            idx.push_back(i - j + 1);
            // OPTIONAL: return idx if only one match required
            j = P[j];
        }
    }
    return idx;
}


int main(){
    string pat = "ababa";
    string txt = "acabaababdababa";
    vi P = prefix(pat);
    vi ans = match(txt, pat, P);
    for (int x : ans){
        cout << x << ' ';
    }
    cout << endl;

    return 0;
}
