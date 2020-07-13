// Date manipulation
//
// Author: Daniel Anderson
// Date: 28-01-2017
// Reliability: 5
// Tested on: UVA11356, SPOJ-CODEIT03, UVA11947, UVA12019, Random tests
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef pair<int,int> pii;
typedef vector<int> vi;
typedef vector<vi> vvi;

//listings:dates
// Date manipulation -- Conversion from Gregorian dates to Julian days. The Julian
// day is the number of days since November 24th 4714 BC. Note that there is no year
// zero in the AD calendar, so 4714 BC corresponds to year -4713. Gregorian dates
// are expressed as {year,month,day}.

// Determine the day of the week for the given Julian date. 0 = Monday ... 6 = Sunday
int day_of_week(int jd) { return jd % 7; }

// Converts the given Gregorian date into the corresponding Julian day
int to_julian(int y, int m, int d) {
    return 1461 * (y + 4800 + (m - 14) / 12) / 4 +
        367 * (m - 2 - (m - 14) / 12 * 12) / 12 -
        3 * ((y + 4900 + (m - 14) / 12) / 100) / 4 + d - 32075;
}

// Converts the given Julian day into the corresponding Gregorian date
tuple<int,int,int> to_gregorian(int jd) {
    int x, n, i, j, y, m, d;
    x = jd + 68569, n = 4 * x / 146097, x -= (146097 * n + 3) / 4;
    i = (4000 * (x + 1)) / 1461001, x -= 1461 * i / 4 - 31;
    j = 80 * x / 2447, d = x - 2447 * j / 80, x = j / 11;
    m = j + 2 - 12 * x, y = 100 * (n - 49) + i + x;
    return make_tuple(y,m,d);
}

// Returns true if the given year is a leap year in the Gregorian calendar
bool leap_year(int y) { return (y % 400 == 0 || (y % 4 == 0 && y % 100 != 0)); }

// Returns the number of days in the given month/year in the Gregorian calendar.
int days_in(int y, int m) { return m == 2 ? 28 + leap_year(y) : 31 - (m-1) % 7 % 2; }
//listings:/dates

namespace problems {
  // Verdict: AC
  namespace UVA11356 {
  string months[] = {"January","February","March","April","May","June","July","August",
                      "September","October","November","December"};
    void solve() {
      int T; cin >> T;
      for (int t=1; t<=T; t++) {
        int yyyy, dd, mm, K; string m, temp;
        getline(cin,temp,'-'); yyyy = stoi(temp);
        getline(cin,temp,'-'); m = temp, mm = find(begin(months), end(months),m) - begin(months) + 1;
        getline(cin,temp); dd = stoi(temp);
        cin >> K;
        int jd = to_julian(yyyy,mm,dd); jd += K;
        tie(yyyy,mm,dd) = to_gregorian(jd);
        printf("Case %d: %04d-%s-%.2d\n", t, yyyy, months[mm-1].c_str(), dd);
      }
    }
  }
  // Verdict: AC
  namespace SPOJ_CODEIT03 {
    string day[] = {"Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"};
    void solve() {
      int T; cin >> T;
      while (T--) {
        int d, m, y; cin >> d >> m >> y;
        int jd = to_julian(y,m,d);
        cout << day[day_of_week(jd)] << '\n';
      }
    }
  }
  // Verdict: AC
  namespace UVA11947 {
    string sign(int m, int d) {
      if (pii(m,d) <= pii(1,20)) return "capricorn";
      else if (pii(m,d) <= pii(2,19)) return "aquarius";
      else if (pii(m,d) <= pii(3,20)) return "pisces";
      else if (pii(m,d) <= pii(4,20)) return "aries";
      else if (pii(m,d) <= pii(5,21)) return "taurus";
      else if (pii(m,d) <= pii(6,21)) return "gemini";
      else if (pii(m,d) <= pii(7,22)) return "cancer";
      else if (pii(m,d) <= pii(8,21)) return "leo";
      else if (pii(m,d) <= pii(9,23)) return "virgo";
      else if (pii(m,d) <= pii(10,23)) return "libra";
      else if (pii(m,d) <= pii(11,22)) return "scorpio";
      else if (pii(m,d) <= pii(12,22)) return "sagittarius";
      else return "capricorn";
    }
    void solve() {
      int N; cin >> N;
      for (int t=1; t<=N; t++) {
        int mm,dd,yyyy;
        scanf("%02d %02d %04d", &mm, &dd, &yyyy);
        int jd = to_julian(yyyy,mm,dd); jd += 40 * 7;
        tie(yyyy,mm,dd) = to_gregorian(jd);
        printf("%d %02d/%02d/%04d %s\n", t, mm, dd, yyyy, sign(mm,dd).c_str());
      }
    }
  }
  // Verdict: AC
  namespace UVA12019 {
    string day[] = {"Monday","Tuesday","Wednesday","Thursday","Friday","Saturday","Sunday"};
    void solve() {
      int T, M, D;
      cin >> T;
      while (T--) {
        cin >> M >> D;
        cout << day[day_of_week(to_julian(2011,M,D))] << '\n';
      }
    }
  }
}

namespace random_tests {
  int T = 1000000;
  void test_conversion() {
    srand(time(0));
    for (int t=1; t<=T; t++) {
      cout << "Running test " << t << "/" << T << "          \r" << flush;
      int Y = rand();
      int M = rand() % 12 + 1;
      int D = rand() % days_in(Y,M) + 1;
      int jd = to_julian(Y,M,D);
      assert(tie(Y,M,D) == to_gregorian(jd));
    }
  }
  int num_days[] = {31,28,31,30,31,30,31,31,30,31,30,31};
  void test_num_days() {
    srand(time(0));
    for (int t=1; t<=T; t++) {
      cout << "Running test " << t << "/" << T << "          \r" << flush;
      int Y = rand();
      int M = rand() % 12 + 1;
      if (leap_year(Y) && M == 2) assert(days_in(Y,M) == 29);
      else assert(days_in(Y,M) == num_days[M-1]);
    }
  }
}

int main() {
  //problems::UVA11356::solve();
  //problems::SPOJ_CODEIT03::solve();
  //problems::UVA11947::solve();
  //problems::UVA12019::solve();
  //brute_force_tests::test_conversion();
  random_tests::test_num_days();
}
