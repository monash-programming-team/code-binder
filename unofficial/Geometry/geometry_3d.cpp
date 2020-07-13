#include <bits/stdc++.h>
using namespace std;

//listings:point
// --------------------------- 3D Computational Geometry -------------------------------
using namespace rel_ops;

#define x first
#define y second

#define cpt const Pt&
#define cpt2 const Pt2&

const double EPS=1e-8;
const double pi=acos(-1);
bool deq(double a,double b) {return abs(a-b)<EPS;}

struct Pt {
	double x=0,y=0,z=0;
	bool operator==(cpt b) const { return deq(x,b.x) && deq(y,b.y) && deq(z,b.z); }
	bool operator<(cpt b) const {
		return x<b.x || (x==b.x && (y<b.y || (y==b.y && z<b.z)));
	}
	double& operator[](int i) {return i==0?x:i==1?y:z;}
	Pt operator+=(cpt b) {return {x+=b.x,y+=b.y,z+=b.z};}
	Pt operator-=(cpt b) {return {x-=b.x,y-=b.y,z-=b.z};}
	Pt operator*=(double c) {return {x*=c,y*=c,z*=c};}
	Pt operator/=(double c) {return {x/=c,y/=c,z/=c};}
};
Pt operator+(cpt a,cpt b) {return {a.x+b.x,a.y+b.y,a.z+b.z};}
Pt operator-(cpt a) {return {-a.x,-a.y,-a.z};}
Pt operator-(cpt a,cpt b) {return {a.x-b.x,a.y-b.y,a.z-b.z};}
Pt operator*(double c,cpt a) {return {c*a.x,c*a.y,c*a.z};}
Pt operator*(cpt a,double c) {return {c*a.x,c*a.y,c*a.z};}
Pt operator/(cpt a,double c) {return {a.x/c,a.y/c,a.z/c};}
double operator*(cpt a,cpt b) {return a.x*b.x+a.y*b.y+a.z*b.z;}

Pt cross(cpt a,cpt b) {return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x};}
double det(cpt a,cpt b,cpt c) {return a*cross(b,c);}
double norm(cpt a) {return a*a;}
double abs(cpt a) {return sqrt(norm(a));}

bool areperp(cpt a,cpt b,cpt p,cpt q) {	return deq((b-a)*(q-p),0); }
bool arepara(cpt a,cpt b,cpt p,cpt q) {	return cross(b-a,q-p)==Pt{0,0,0}; }
//listings:/point

//listings:great_circle
double greatcircledist(cpt a,cpt b) {
	return abs(a)*acos((a*b)/(abs(a)*abs(b)));
}
//listings:/great_circle

//listings:geometry
typedef pair<double,double> Pt2;

double det(cpt2 a,cpt2 b) {return a.x*b.y-a.y*b.x;}

// Finds a line that is perpendicular to two lines. Divides by 0 if they are
// parallel. The if statments below ensure the resulting line intersects
// the lines taken as segments. The first point returned lies on ab and the
// second on pq. If the given lines intersect then the points returned are the same.
pair<Pt,Pt> perpline(Pt a,Pt b,Pt p,Pt q) {
	Pt ab=b-a,qp=p-q,ap=p-a;
	Pt2 c{ab*ab,ab*qp},d{ab*qp,qp*qp},e{ab*ap,qp*ap};//[c,d] is Gram matrix
	double s=det(e,d)/det(c,d),t=det(c,e)/det(c,d);
	if (-EPS<s && s<1+EPS //Answer intersects ab
		&& -EPS<t && t<1+EPS) //Answer intersects pq
		return {a+s*ab,p-t*qp};
	return {Pt{NAN,NAN,NAN},{}};
}

//Distance between line and point (Infinite line and line segment respectively)
double distlinept(cpt a,cpt b,cpt p) { return abs(cross(b-a,p-a))/abs(b-a); }
double distsegpt(Pt a,Pt b,Pt p) {
	b-=a;p-=a;	double sp=b*p/norm(b);	Pt closest;
	if (sp>=0) {
		if (sp>1) closest=b;
		else closest=sp*b;
	}
	return abs(closest-p); // Note that actual closest Pt on line is closest + a
}

//Project p onto the plane through the origin spanned by a and b. Coordinates
//are given with respect to the basis {a,b}. Divides by 0 if a and b are
//parallel.
Pt2 projectplanept(cpt a,cpt b,cpt p) {
	Pt2 c{a*a,a*b},d{a*b,b*b},e{a*p,b*p};
	return {det(e,d)/det(c,d),det(c,e)/det(c,d)};
}

//Divides by 0 if a, b and c are collinear.
double disttrianglept(Pt a,Pt b,Pt c,Pt p) {
	b-=a;c-=a;p-=a;
	double s,t; tie(s,t)=projectplanept(b,c,p);
	if (0<s && 0<t && s+t<1) return abs(s*b+t*c-p); // Projection within tri.
	return min({distsegpt(a,b,p),distsegpt(a,c,p),distsegpt(b,c,p)});
}

//Distance between two finite lines. Modify perpline to get infinite lines
double distfinitelineline(cpt a,cpt b,cpt p,cpt q) {
	if (!arepara(a,b,p,q)) {
		Pt u,v;	tie(u,v)=perpline(a,b,p,q);
		if (!std::isnan(u.x)) return abs(v-u);
	}
	return min({distsegpt(a,b,p),distsegpt(a,b,q),distsegpt(p,q,a),distsegpt(p,q,b)});
}

//Rotate a point around a line by theta radians. Anticlockwise when looking from b to a.
Pt rotatelinept(Pt a,Pt b,double theta,Pt p) {
	b-=a;p-=a;	b/=abs(b);	double C=cos(theta);
	return C*p+(1-C)*(b*p)*b+sin(theta)*cross(b,p)+a;
}

//Use quaternions when composition of 3D rotations is required. Note that both a
//quaternion and its negative represent the same rotation.
typedef pair<double,Pt> Quaternion;
#define cq const Quaternion&

// Gives the rotation equivalent to doing the b rotation then the a rotation.
Quaternion operator*(cq a,cq b) {
	return {a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x+cross(a.y,b.y)};
}

double norm(cq a) {return norm(a.x)+norm(a.y);}
double abs(cq a) {return sqrt(norm(a));}
Quaternion operator/(cq a,double c) {return {a.x/c,a.y/c};}

// Careful of divide by zero if you invert this
Quaternion quaternionforrotation(Pt a,double theta) {
	return {cos(theta/2),sin(theta/2)/abs(a)*a};
}

Pt rotatept(Quaternion q,cpt p) {
	q=q/abs(q); // Need this only if quaternion not already normalized.
	return p+cross(2*q.y,cross(q.y,p)+q.x*p);
}

// 3D Convex Hull O(n^2)
// faces is an array of triangles covering the convex hull. f is the number of
// faces. Edges and Tris store indices of p. For any Tri of the hull,
// (p[b]-p[a]) X (p[c]-p[a]) points outward.
// Fun fact: Any triangulation of a (non-degenerate) polyhedron with n vertices
// has 3*(n-2) edges and 2*(n-2) faces.
// If f is two after running convexhull then all points lie in the plane
// described by the two faces but they are not necessarily touching the
// triangle. If f is 0 then all points are collinear.
namespace Hull {  // Set maxn to max number of points
	const int maxn=1000; int f,inh[2*maxn],in[maxn],out[maxn],modif[maxn];
	struct Tri {int a,b,c;} faces[2*maxn];
	struct Edge { int u,v,f[2];	int& operator[](int i) {return i?v:u;} } edges[3*maxn];
	void convexhull(const vector<Pt>& p) {
		int n=p.size(),m=0,i,j,k; f=0; fill(modif,modif+n,-1);
		for (i=1;i<n;i++) {
			if (m==0 && p[i]!=p[0]) edges[m++]={0,i,0,0};
			bool use=m==1 && cross(p[edges[0][1]]-p[0],p[i]-p[0])!=Pt{0,0,0};
			for (j=0;j<f;j++) {
				Tri &t=faces[j];
				if (inh[j] && det(p[t.a]-p[i],p[t.b]-p[i],p[t.c]-p[i])<-EPS)
					inh[j]=0,use=1;
			}
			if (!use) continue;
			for (j=0,k=0;j<m;j++) {
				int nk=1; Edge &e=edges[j];
				if (m==1 || (nk=(int)inh[e.f[0]]+inh[e.f[1]])==1)
					for (int c=0;c<2;c++) if (m==1 || !inh[e.f[c]]) {
						for (;inh[k] && k<f;k++);
						Tri &t=faces[k]={e[c],e[1-c],i};
						e.f[c]=k, in[t.b]=out[t.a]=k, modif[t.b]=i,	k++;
					}
				if (nk==0) e=edges[--m], j--;
			}
			bool reset=f==2; f=max(f,k);
			for (j=0;j<n;j++) if (modif[j]==i)
				edges[m++]={i,j,out[j],in[j]},inh[in[j]]=1;
			if (reset) i=0;
		}
		for (i=0,j=0;j<f;j++) if (inh[j]) faces[i++]=faces[j];
		f=i;
	}
}

struct Tri {Pt a,b,c;};

//Signed Volume of polyhedron. Positive when (b-a) X (c-a) points outward for each Tri.
double volume(const vector<Tri>& poly) {
	double r=0;
	for (const Tri &t:poly) r+=det(t.a,t.b,t.c);
	return r/6;
}
//Surface area of polyhedron
double surfacearea(const vector<Tri>& poly) {
	double r=0;
	for (const Tri &t:poly) r+=abs(cross(t.b-t.a,t.c-t.a));
	return r/2;
}
//listings:/geometry

//listings:dt
// Delauney Triangulation O(n^2)
// Triangulation of a set of points so that no point p is inside the
// circumcircle of any triangle. Maximizes the minimum angle of all angles of
// the triangles in the triangulation. Each Tri in the result holds 3 indices of
// p. The indices are such that det(p[b]-p[a],p[c]-p[a]) is positive. If all
// points are collinear, then the triangulation will be empty.
vector<Hull::Tri> delauneytriangulation(const vector<Pt2>& p) {
	using namespace Hull;
	vector<Pt> q(p.size());
	for (int i=0;i<p.size();i++) q[i]={p[i].x,p[i].y,-norm(p[i].x)-norm(p[i].y)};
	convexhull(q);
	for (int i=0;i<f;i++) {
		Hull::Tri &t=faces[i];
		if (cross(q[t.b]-q[t.a],q[t.c]-q[t.a]).z<EPS) faces[i--]=faces[--f];
	}
	return {faces,faces+f};
}
//listings:/dt


int main() {
	ios::sync_with_stdio(0); cin.tie(0);

}
