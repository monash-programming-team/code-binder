// Numerical integration and differentiation
//
// Author: Darcy's code binder
// Date: 19-01-2017
//
// Assumed tested and reliable from Darcy's code binder. (If anyone can actually
// find a problem that uses these then feel free to verify).
//
#include<bits/stdc++.h>
using namespace std;

#define X first
#define Y second

typedef long long ll;
typedef vector<int> vi;

//listings:numerical
// Numerical integration. Integrate f(x) for x in [a,b]. n is the number
// of intervals (it must be even). If K is an upper bound on the 4th derivative
// of f for all x in [a,b], then the error is bounded by (K h^5) / (180 n^4)
template<typename F> double integrate(F f, double a, double b, int n) {
  double ans = f(a) + f(b), h = (b-a)/n;
  for (int i=1; i<n; i++) ans += f(a+i*h) * (i%2 ? 4 : 2);
  return ans * h / 3;
}

// Numerical differentiation. h is the step size. Error is O(h^4)
template<typename F> double differentiate(F f, double x, double h) {
  return (-f(x+2*h) + 8*(f(x+h) - f(x-h)) + f(x-2*h)) / (12*h);
}
//listings:/numerical

const double EPSILON = 1e-6;

int main() {
  auto f = [](double x) { return x*x; };
  
  for (double x = 0; x <= 1000; x += 0.01) {
    double df = differentiate(f, x, 0.001);
    double exact = 2*x;
    assert(abs(df-exact) < EPSILON);
  }
  
  for (double a = -100, b = 0; b <= 100; a += 0.1, b += 0.1) {
    double area = integrate(f, a, b, 1000);
    double exact = 1.0/3.0 * (b*b*b - a*a*a);
    assert(abs(area-exact) < EPSILON);
  }
}