namespace ECNA07_A {
  void solve() {
    int n, caseno=1;
    while (cin >> n && n) {
      int ans = 0;
      while (n--) {
        string s;
        cin >> s;
        ans += to_arabic(s);
      }
      cout << "Case " << to_roman(caseno++) << ": " << to_roman(ans) << endl;
    }
  }
}

namespace SPOJ_ROMAN008 {
  void solve() {
    string a, b;
    char op;
    while (cin >> a) {
      cin >> b >> op;
      int x = to_arabic(a), y = to_arabic(b);
      if (op == '+')
        cout << to_roman(x+y) << endl;
      else if (op == '-')
        cout << to_roman(x-y) << endl;
      else if (op == '*')
        cout << to_roman(x*y) << endl;
      else if (op == '/')
        cout << to_roman(x/y) << endl;
      else
        cout << to_roman(x%y) << endl;
    }
  }
}

namespace SPOJ_ROMANN {
  void solve() {
    int K;
    cin >> K;
    for (int k=1; k<=K; k++) {
      string s;
      cin >> s;
      cout << "Case #" << k << ": " << to_arabic(s) << '\n';
    }
  }
}

namespace UVA_759 {
  void solve() {
    set<string> valid;
    for (int i=1; i<4000; i++)
      valid.insert(to_roman(i));
    string s;
    while (cin >> s) {
      if (valid.count(s) == 0)
        cout << "This is not a valid number\n";
      else
        cout << to_arabic(s) << '\n';
    }
  }
}

namespace UVA_11616 {
  void solve() {
    char c;
    while (cin >> c) {
      cin.putback(c);
      if (c >= '0' && c <= '9') {
        int x;
        cin >> x;
        cout << to_roman(x) << '\n';
      } else {
        string s;
        cin >> s;
        cout << to_arabic(s) << '\n';
      }
    }
  }
}

namespace UVA_12397 {
  void solve() {
    map<char,int> val = {{'I',1},{'V',2},{'X',2},{'L',2},{'C',2},{'D',3},{'M',4}};
    int N;
    while (cin >> N) {
      int ans = 0;
      for (const char c : to_roman(N))
        ans += val[c];
      cout << ans << '\n';
    }
  }
}
