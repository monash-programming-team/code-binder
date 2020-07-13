namespace CF_284A {
  void solve() {
    ll p;
    cin >> p;
    cout << phi(phi(p)) << endl;
  }
}

namespace SPOJ_DISCRT1 {
  void solve() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    ll n, k, a;
    cin >> n >> k >> a;
    auto res = discrete_root(n, k, a);
    sort(res.begin(), res.end());
    cout << res.size() << '\n';
    for (auto x : res)
      cout << x << '\n';
  }
}

namespace SPOJ_ETF {
  void solve() {
    int T;
    cin >> T;
    while (T--) {
      int n;
      cin >> n;
      cout << phi(n) << endl;
    }
  }
}

namespace SPOJ_MOD {
  void solve() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    ll x, z, k;
    while (cin >> x >> z >> k, x > 0) {
      auto ans = discrete_log(x, k, z);
      if (ans == -1)
        cout << "No Solution" << '\n';
      else
        cout << ans << '\n';
    }
  }
}

namespace SPOJ_PROOT {
  void solve() {
    ll p, n;
    while (cin >> p >> n, p > 0) {
      ll tot = p - 1;
      auto fact = slow_factors(tot);
      sort(fact.begin(),fact.end());
      fact.resize(unique(fact.begin(),fact.end(),fact.begin())-fact.begin());
      for (int i = 0; i < n; i++) {
        ll r;
        cin >> r;
        if (is_proot(p, r, tot, fact))
          cout << "YES" << '\n';
        else
          cout << "NO" << '\n';
      }
    }
  }
}

namespace UVA_374 {
  void solve() {
    int B, P, M;
    while (cin >> B >> P >> M)
      cout << expmod(B, P, M) << endl;
  }
}

namespace UVA_543 {
  void solve() {
    // Pre-compute
    const int maxn = 20000000;
    sieve(maxn);
    fast_sieve(maxn);
    // Answer queries
    int n;
    while (cin >> n, n != 0) {
      for (int p : pr) {
        if (isprime[n - p]) {
          assert(is_prime(n - p));
          assert(binary_search(pr.begin(), pr.end(), n - p));
          assert(trial_division(n - p));
          printf("%d = %d + %d\n", n, p, n - p);
          break;
        }
      }
    }
  }
}

namespace UVA_583 {
  void solve() {
    ll n;
    while (cin >> n, n != 0) {
      bool neg = false;
      if (n < 0) {
        n *= -1;
        neg = true;
      }
      auto F2 = slow_factors(n);
      for (int f : F2)
        assert(is_prime(f) && n % f == 0);
      sort(F2.begin(), F2.end());
      if (neg)
        cout << '-';
      cout << n << " =";
      if (neg)
        cout << " -1 x";
      for (int i=0; i<(int)F2.size(); i++) {
        cout << ' ' << F2[i];
        if (i < (int)F2.size() - 1)
          cout << " x";
      }
      cout << endl;
    }
  }
}

namespace UVA_756 {
  void solve() {
    const char* msg = "Case %d: the next triple peak occurs in %d days.\n";
    const int g = __gcd(23, __gcd(28, 33)), lcm = 23 * 28 * 33 / g;
    int p, e, i, d, caseno = 1;
    while (cin >> p >> e >> i >> d, p != -1) {
      vi m = {23, 28, 33}, a = {p, e, i};
      int x = cra(a, m);
      while (x <= d)
        x += lcm;
      printf(msg, caseno++, x - d);
    }
  }
}

namespace UVA_1230 {
  void solve() {
    int c;
    cin >> c;
    while (c--) {
      ll x, y, n;
      cin >> x >> y >> n;
      cout << expmod(x, y, n) << endl;
    }
  }
}

namespace UVA_10090 {
  void solve() {
    ll n, c1, c2, n1, n2, x0, y0;
    while (cin >> n, n > 0) {
      cin >> c1 >> n1 >> c2 >> n2;
      ll G = gcd(n1, n2, x0, y0);
      if (n % G != 0) {
        cout << "failed" << endl;
        continue;
      } else {
        ll x = x0 * n / G, y = y0 * n / G;
        ll d1 = n2 / G, d2 = n1 / G;
        // Make x non-negative
        if (x < 0) {
          ll q = -(x-d1+1)/d1;
          x += q * d1, y -= q * d2;
          if (y < 0) {
            cout << "failed" << endl;
            continue;
          }
        }
        // Make y non-negative
        if (y < 0) {
          ll q = -(y-d2+1)/d2;
          y += q * d2, x -= q * d1;
          if (x < 0) {
            cout << "failed" << endl;
            continue;
          }
        }
        // Minimise c1 x + c2 y
        if (c1 * d1 > c2 * d2) {
          ll q = x / d1;
          x -= q * d1, y += q * d2;
        } else if (c1 * d1 < c2 * d2) {
          ll q = y / d2;
          y -= q * d2, x += q * d1;
        }
        // Output
        cout << x << ' ' << y << endl;
      }
    }
  }
}

namespace UVA_10104 {
  void solve() {
    int A, B;
    ll X, Y, D;
    while (cin >> A) {
      cin >> B;
      D = gcd(A, B, X, Y);
      cout << X << ' ' << Y << ' ' << D << endl;
    }
  }
}

namespace UVA_10179 {
  void solve() {
    ll n;
    while (cin >> n, n > 0)
      cout << phi(n) << endl;
  }
}

namespace UVA_10392 {
  void solve() {
    // NOTE: Only works with __int128. WA if you don't have it. You'll
    // also need to choose the biggest test-numbers for is_prime.

    ll n;
    while (cin >> n, n != -1) {
      auto res = factor_huge(n);
      sort(res.begin(), res.end());
      for (auto x : res) {
        cout << "    " << x << '\n';
      }
      cout << '\n';
    }
  }
}

namespace UVA_10394 {
  void solve() {
    // Pre-compute
    const int maxn = 20000000;
    sieve(maxn);
    fast_sieve(maxn);
    vector<pii> twins;
    vector<bool> seen(maxn);
    for (int p : pr)
      seen[p] = true;
    for (int x = 1; x <= maxn; x++)
      if (isprime[x] && isprime[x+2]) {
        assert(is_prime(x) && is_prime(x+2));
        assert(seen[x] && seen[x+2]);
        assert(trial_division(x) && trial_division(x+2));
        twins.emplace_back(x, x+2);
      }
    // Answer input
    int S;
    while (cin >> S)
      printf("(%d, %d)\n", twins[S-1].X, twins[S-1].Y);
  }
}

namespace UVA_11904 {
  const int MOD = 1e9 + 7;
  const int MAXN = 1000100;
  ll fact[MAXN];

  void solve() {
    fact[0] = fact[1] = 1;
    for (int i=2; i<MAXN; i++)
      fact[i] = (fact[i-1] * i) % MOD;
    int T, c = 1;
    cin >> T;
    while (T--) {
      int n;
      cin >> n;
      vi k(n);
      for (int i = 0; i < n; i++)
        cin >> k[i];
      ll sum = 0, ans = 1;
      for (auto x : k) {
        ans = (ans * fact[sum + x - 1]) % MOD;
        ans = (ans * inv((fact[x - 1] * fact[sum]) % MOD, MOD)) % MOD;
        sum += x;
      }
      cout << "Case " << c++ << ": " << ans << endl;
    }
  }
}
