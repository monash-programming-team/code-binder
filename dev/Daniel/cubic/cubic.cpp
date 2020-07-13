// Hashing custom types example
//
// Author: Daniel. Based on Darcy's code binder
// Date: 19-01-2017
//
// Assumed reliable from Darcy's code binder. If you find any problems that
// require this, you should verify it.
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;

//listings:cubic
// Cubic equation solver. Solves ax^3 + bx^2 + cx + d = 0.
// a must be non-zero, does NOT work well when a is NEAR 0.
const double EPSILON = 1e-8, PI = acos(-1);

template<typename T> vector<T> cubic(T a, T b, T c, T d) { 
  b /= a, c /= a, d /= a;  // Make sure T is non-integral (double or long double)!
  T q = (b*b - 3*c)/9, r = (2*b*b*b - 9*b*c + 27*d) / 54, z = r*r - q*q*q;
  if (z <= EPSILON) {
    vector<T> sol;  T theta = acos(r/pow(q,1.5));
    for (int i=0; i<3; i++) sol.push_back(-2*sqrt(q)*cos((theta+i*2*PI)/3) - b/3);
    return sol;
  }
  T s = cbrt(sqrt(z)+abs(r)); s = (s + q/s) * (r < 0 ? 1 : -1) - b/3;
  return {s};
}
//listings:/cubic

int main() {
  auto sol = cubic(1.0,1.0,1.0,1.0);
  for (auto x : sol) cout << x << endl;
}
