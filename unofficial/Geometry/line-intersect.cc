/**
 * Line Intersection for 2D
 * @author Ryan Hechenberger
 * @date 2016-03-29
 * Reliability: 5
 * Tested On: uva191, uva378, uva866, uva10902, uva11343
 * 
 * Computes the shortest path between two points given a list of points
 * (duplicates permitted)
 * 
 * Complexity: O( 1 )
 */
#define X real()
#define Y imag()

using namespace std;
typedef long long ll;
typedef complex<double/ll> pt;
/// Required operands, same for both versions, choose ll for integer
double/ll cross(pt a, pt b) { return (conj(a) * b).Y; }
double/ll dot(pt a, pt b) { return (conj(a) * b).X; }
constexpr double EPS = 1e-9; /// Floating-point only
/// Integer fraction class
struct frac { ll num, den; };

/**
 * @param a Line 1 Point 1
 * @param b Line 1 Point 2
 * @param c Line 2 Point 1
 * @param d Line 2 Point 2
 * @param r Intersection point
 * @param seg Is a segment line (true) or infinate line (false)
 * @returns  0 on no intersection between the lines, 1 on an intersection,
 *          -1 on an infinate number of intersections
 * 
 * Will determain if lines intersect, and returns the intersection point in
 * r if there is only one intersection.
 * Principle from: http://stackoverflow.com/questions/563198/
 */
int line_intersect(pt a, pt b, pt c, pt d, pt& r, bool seg)
{
	pt u = b - a, v = d - c;
	double x1 = cross(u, v), x2 = cross(c - a, u);
	if (abs(x1) > EPS) {
		double s1 = cross(c - a, v) / x1, s2 = x2 / x1;
		if (seg && (s1 < -EPS || s1 > 1+EPS || s2 < -EPS || s2 > 1+EPS)) return 0;
		r = pt(a.X + s1 * u.X, a.Y + s1 * u.Y);
		return 1;
	}
	if (abs(x2) > EPS) return 0; // parallel
	// coliniear onwards
	if (!seg) return -1;
	double d1 = dot(u, v), d2 = dot(u, u);
	if (d1 < 0) { c = d; d1 = -d1; }
	double t0 = dot(c - a, u) / d2, t1 = t0 + d1 / d2;
	if (t1 < -EPS || t0 > 1+EPS) return 0;
	if (t1 > EPS && t0 < 1-EPS) return -1;
	r = t1 <= EPS ? a : c;
	return 1;
}

/**
 * @param a Line 1 Point 1
 * @param b Line 1 Point 2
 * @param c Line 2 Point 1
 * @param d Line 2 Point 2
 * @param x The x intersection point, is not normalised
 * @param y The y intersection point, is not normalised
 * @param seg Is a segment line (true) or infinate line (false)
 * @returns  0 on no intersection between the lines, 1 on an intersection,
 *          -1 on an infinate number of intersections
 * 
 * The integer version of the line intersect. Will set x and y to the
 * intersection if there is only a single intersect point.
 * Principle from: http://stackoverflow.com/questions/563198/
 */
int line_intersect(pt a, pt b, pt c, pt d, frac& x, frac& y, bool seg)
{
	pt u = b - a, v = d - c;
	ll x1 = cross(u, v);
	if (x1 < 0) { x1 = -x1; a = b; u = -u; }
	ll x2 = cross(c - a, u);
	if (x1 != 0) {
		ll s = cross(c - a, v);
		if (seg && (s < 0 || s > x1 || x2 < 0 || x2 > x1)) return 0;
		x = { a.X * x1 + s * u.X, x1 }; y = { a.Y * x1 + s * u.Y, x1 };
		return 1;
	}
	if (x2 != 0) return 0; // parallel
	// coliniear onwards
	if (!seg) return -1;
	ll d1 = dot(u, v), d2 = dot(u, u); // both will be made positive non-zeros
	if (d1 < 0) { d1 = -d1; c = d; }
	ll t0 = dot(c - a, u), t1 = t0 + d1;
	if (t1 < 0 || t0 > d2) return 0;
	if (t1 > 0 && t0 < d2) return -1;
	x = { t1 <= 0 ? a.X : c.X, 1 }; y = { t1 <= 0 ? a.Y : c.Y, 1 };
	return 1;
}

int main()
{
	int casen;
	cin >> casen;
	printf("INTERSECTING LINES OUTPUT\n");
	while (casen-- > 0) {
		point a, b, c, d, r;
		cin >> a.x >> a.y >> b.x >> b.y >> c.x >> c.y >> d.x >> d.y;
		int t = line_intersect(a, b, c, d, r, false);
		if (t == 0)
			printf("NONE\n");
		else if (t == 1)
			printf("POINT %.2f %.2f\n", r.x, r.y);
		else
			printf("LINE\n");
	}
	printf("END OF OUTPUT\n");
	
	return 0;
}

bool in_poly(vector<pt> P, pt v, bool boundary)
{
	bool in = false;
	for (int i = 0, j = (int)P.size() - 1; i < (int)P.size(); j = i++) {
		pt a = P[i], b = P[j];
		if (abs(ccw(v, a, b))<EPS && max(norm(v-a), norm(v-b)) <= norm(a-b)+EPS)
			return boundary;
		if ( min(a.Y, b.Y) <= v.Y+EPS && v.Y+EPS < max(a.Y, b.Y) &&
				(v.Y - a.Y) / (b.Y - a.Y) * (b.X - a.X) < v.X )
			in = !in;
	}
	return in;
}