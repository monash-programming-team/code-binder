#include <bits/stdc++.h>

//listings:geometry
#define x first
#define y second

#define cpt const Pt&
#define cpt2 const Pt2&

using namespace std;
using namespace rel_ops;

const double EPS=1e-8;
const double pi=acos(-1);
bool deq(double a,double b) {return abs(a-b)<EPS;}

struct Pt {
	double x=0,y=0,z=0;
	bool operator==(cpt b) const {
		return deq(x,b.x) && deq(y,b.y) && deq(z,b.z);
	}
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
Pt cross(cpt a,cpt b) { return {a.y*b.z-a.z*b.y,a.z*b.x-a.x*b.z,a.x*b.y-a.y*b.x}; }
double det(cpt a,cpt b,cpt c) {return a*cross(b,c);}
double norm(cpt a) {return a*a;}
double abs(cpt a) {return sqrt(norm(a));}

struct Pt2 { double x,y; };
double det(cpt2 a,cpt2 b) {return a.x*b.y-a.y*b.x;}
//listings:/geometry

// Finds a line that is perpendicular to two lines. Divides by 0 if they are
// parallel. The if statments below ensure the resulting line intersects
// the lines taken as segments. The first point returned lies on ab and the
// second on pq. If the given lines intersect then the points returned are the
// same.
pair<Pt,Pt> perpline(Pt a,Pt b,Pt p,Pt q) {
	Pt ab=b-a,qp=p-q,ap=p-a;
	Pt2 c{ab*ab,ab*qp},d{ab*qp,qp*qp},e{ab*ap,qp*ap};//[c,d] is Gram matrix
	double s=det(e,d)/det(c,d),t=det(c,e)/det(c,d);

	if (-EPS<s && s<1+EPS //Answer intersects ab
		&& -EPS<t && t<1+EPS) //Answer intersects pq
		return {a+s*ab,p-t*qp};

	return {Pt{NAN,NAN,NAN},{}};
}

//Distance between infinite line and point
double distlinept(cpt a,cpt b,cpt p) {
	return abs(cross(b-a,p-a))/abs(b-a);
}

double distfinitelinept(Pt a,Pt b,Pt p) {
	b-=a;p-=a;
	double sp=b*p/norm(b);
	Pt closest;
	if (sp>=0) {
		if (sp>1) closest=b;
		else closest=sp*b;
	}
	return abs(closest-p); // Note that actual closest Pt on line is closest + a
}

//Are lines perpendicular
bool areperp(cpt a,cpt b,cpt p,cpt q) {
	return deq((b-a)*(q-p),0);
}

//Are lines parallel
bool arepara(cpt a,cpt b,cpt p,cpt q) {
	return cross(b-a,q-p)==Pt{0,0,0};
}

// TODO: Point to plane distance

// TODO: Point to triangle distance

//Distance between two finite lines
//Modify perpline to get infinite lines
double distfinitelineline(cpt a,cpt b,cpt p,cpt q) {
	if (!arepara(a,b,p,q)) {
		Pt u,v;
		tie(u,v)=perpline(a,b,p,q);
		if (!std::isnan(u.x)) return abs(v-u);
	}

	return min({
			distfinitelinept(a,b,p),
			distfinitelinept(a,b,q),
			distfinitelinept(p,q,a),
			distfinitelinept(p,q,b)
	});
}

//Rotate a point around a line by theta radians. Rotation is anticlockwise when
//looking from b to a.
Pt rotatelinept(Pt a,Pt b,double theta,Pt p) {
	b-=a;p-=a;
	b/=abs(b);
	double C=cos(theta);
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
	            // This could be the case after many multiplications.
	return p+cross(2*q.y,cross(q.y,p)+q.x*p);
}

// TODO: Volume of convex polyhedron

// 3D convex hull O(n^2)
// For any Tri of the hull, (b-a) X (c-a) points outward
namespace Hull {
	struct Tri {
		Pt a,b,c;
		Pt& operator[](int i) {return i==0?a:i==1?b:c;}
	};
	const int maxn=1000; // Set this to max number of points
	// Fun fact: Any triangulation of a (non-degenerate) polyhedron with n vertices
	// has 3*(n-2) edges and 2*(n-2) faces.
	int n,f; // Number of faces
	Tri faces[2*maxn];
	// 3 triangles a triangle shares a side with
	//int adj[2*maxn][3]{{2,3,1},{0,3,2},{1,3,0},{0,2,1}};
	int adj[2*maxn][3],prep[4][3]{{2,3,1},{0,3,2},{1,3,0},{0,2,1}};
	bool inhull[2*maxn]{1,1,1,1};

	int nc(int i) {
		return (int)inhull[adj[i][0]]+inhull[adj[i][1]]+inhull[adj[i][2]];
		int r=0;for (int j=0;j<3;j++) r+=inhull[adj[i][j]];return r;
	}

	void convexhull(const vector<Pt>& p) {
		fill(inhull,inhull+4,1);
		fill(inhull+4,inhull+2*p.size(),0);
		copy(prep[0],prep[0]+12,adj[0]);


		n=p.size();
		Pt tetra[4];
		int i=0,j=0;
		for (;i<n && j<4;i++) if (
				j==0 ||
				(j==1 && p[i]!=tetra[0]) ||
				(j==2 && cross(tetra[1]-p[0],p[i]-p[0])!=Pt{0,0,0}) ||
				(j==3 && !deq(0,det(tetra[1]-p[0],tetra[2]-p[0],tetra[3]-p[i]))))
			tetra[j++]=p[i];

		//cerr << j << endl;
		if (j<4) {
			f=0;
			return;
		}
		assert(j==4);
		if (det(tetra[1]-tetra[0],tetra[2]-tetra[0],tetra[3]-tetra[0])<0)
			swap(tetra[2],tetra[3]);

		faces[0]={tetra[0],tetra[2],tetra[1]};
		faces[1]={tetra[0],tetra[1],tetra[3]};
		faces[2]={tetra[0],tetra[3],tetra[2]};
		faces[3]={tetra[1],tetra[2],tetra[3]};
		f=4;

		
		for (i=0;i<n;i++) {
		/*for (int i=0;i<f;i++) {
			//cerr << "triangle: " << i << endl;
			//cerr << "adj: " << adj[i][0] << ' ' << adj[i][1] << ' ' << adj[i][2] << endl;
			for (int j=0;j<3;j++) cerr << faces[i][j].x << ' ' << faces[i][j].y << ' ' << faces[i][j].z << endl;
			cerr << endl;
		}*/

		//cerr << "break--\n";

			bool use=0;
			for (j=0;j<f;j++) {
				Tri &t=faces[j];
				if (inhull[j] && det(t.a-p[i],t.b-p[i],t.c-p[i])<-EPS) {
					use=1;
					inhull[j]=0;
				}
			}

			/*cerr << (use?"in":"out") << endl;
			for (int j=0;j<f;j++) cerr << inhull[j] << ' ';
			cerr << endl;*/
			if (!use) continue;
			//cerr << "yo\n";
			
			int s=0,e=0;
			for (;1;s++) for (j=0;j<3;j++) if (inhull[s] && !inhull[adj[s][j]]) {
					e=j;
					goto brk;
			}
brk:
			//cerr << s << ' ' << e << endl;
			//for (int i=0;i<4;i++) cerr << inhull[i] << ' ';
			//cerr << endl;
			int l=-1,t=s,ce=e,firs;
			for (j=0;l==-1 || t!=s || ce!=e;j++) {
				//cerr << "hey: " << t << ' '<< ce << endl;
				// add tri
				while (inhull[j]) j++;
				//inhull[j]=1;
				//cerr << "j: " << j << endl;
				/*if (t==3 && ce==2) {
					cerr << "ayo\n";
					   for (int i=0;i<3;i++) {
					   cerr << faces[t][i].x << ' ' << faces[t][i].y << ' ' << faces[t][i].z << endl;
					   //cerr << endl;
					   }
					cerr << "hol up\n";

				}*/
				faces[j].a=faces[t][(ce+1)%3];
				faces[j].b=faces[t][ce];
				faces[j].c=p[i];
				adj[j][0]=t;
				adj[t][ce]=j;
				// match to last it l!=-1
				if (l!=-1) {
					adj[j][1]=l;
					adj[l][2]=j;
				}
				else firs=j;
				// get next edge
				l=j;
				for (ce=(ce+1)%3;inhull[adj[t][ce]];) {
					for (int x=0;x<3;x++) if (adj[adj[t][ce]][x]==t) {
						t=adj[t][ce];
						ce=(x+1)%3;
						break;
					}
				}
			}
			// join first and last edge
			adj[firs][1]=l;
			adj[l][2]=firs;
			for (int k=0;k<j;k++) inhull[k]=1;
			f=max(f,j);
		}
		for (i=0,j=0;j<f;j++) if (inhull[j]) {
			faces[i++]=faces[j];
			//copy(adj[j],adj[j]+3,adj[i++]);
		}
		f=i;
	}
}

// TODO: Delaunay triangulation

// TODO: Voronori diagram

namespace Test {
	int rn() {
		int range=100;
		return rand()%(2*range+1)-range;
	}

	void testperpline() {
		int ct[9];
		int range=5,found=0;
		for (ct[0]=-range;ct[0]<range;ct[0]++)
		for (ct[1]=-range;ct[1]<range;ct[1]++)
		for (ct[2]=-range;ct[2]<range;ct[2]++)
		for (ct[3]=-range;ct[3]<range;ct[3]++)
		for (ct[4]=-range;ct[4]<range;ct[4]++)
		for (ct[5]=-range;ct[5]<range;ct[5]++)
		for (ct[6]=-range;ct[6]<range;ct[6]++)
		for (ct[7]=-range;ct[7]<range;ct[7]++)
		for (ct[8]=-range;ct[8]<range;ct[8]++)
		{
			Pt base{rn(),rn(),rn()};
			Pt pts[4]{
				base,
				base+Pt{ct[0],ct[1],ct[2]},
				base+Pt{ct[3],ct[4],ct[5]},
				base+Pt{ct[6],ct[7],ct[8]}
			};
			pair<Pt,Pt> ans=perpline(pts[0],pts[1],pts[2],pts[3]);
			assert(std::isnan(ans.x.x) || (
				deq(distfinitelinept(pts[0],pts[1],ans.x),0)
				&& deq(distfinitelinept(pts[2],pts[3],ans.y),0)
				&& areperp(pts[0],pts[1],ans.x,ans.y)
				&& areperp(pts[2],pts[3],ans.x,ans.y)
				));
			if (!std::isnan(ans.x.x)) found++;
		}
		assert(found==351001068);
	}

	void testrotations() {
		Pt ans={4,5,10};
		assert(rotatelinept({1,1,0},{1,1,1},pi,{-2,-3,10})==ans);
		ans-={1,1,0};
		assert(rotatept(quaternionforrotation({0,0,1},pi),{-3,-4,10})==ans);
		for (int i=0;i<1000000;i++) {
			Pt a{rand(),rand(),rand()},b{rand(),rand(),rand()};
			double ang1=(double)rand()/RAND_MAX*4*pi;
			double ang2=(double)rand()/RAND_MAX*4*pi;
			Pt s={rand(),rand(),rand()};
			Pt t1=rotatelinept({0,0,0},b,ang2,rotatelinept({0,0,0},a,ang1,s));
			Quaternion combined=
				quaternionforrotation(b,ang2)*quaternionforrotation(a,ang1);
			Pt t2=rotatept(combined,s);
			assert(abs(t2-t1)/(abs(t2)+abs(t1))<EPS);
		}
	}

	void testconvexhull() {
		cerr << "start test hull\n";
		vector<Pt> inp{{0,0,-5},{-5,-7,0},{5,-7,0},{0,20,0},{0,0,8}};
		Hull::convexhull(inp);
		for (int i=0;i<Hull::f;i++) {
			for (int j=0;j<3;j++) cerr << Hull::faces[i][j].x << ' ' << Hull::faces[i][j].y << ' ' << Hull::faces[i][j].z << endl;
			cout << endl;
		}
	}

	void solve() {
		Pt a{1,2,3},b{4,5,9};
		a=a*5.0;
		Pt ans{5,10,15};
		assert(a==ans);
		assert(!arepara(Pt(),a,Pt(),b));
		assert(deq(a*b,205));
		assert(cross(a,b)==-cross(b,a));
		ans=Pt{15,15,-15};
		assert(cross(a,b)==ans);
		assert(deq(abs(a),5*sqrt(14)));
		assert(deq(det({1,0,0},{0,1,0},{0,0,1}),1));
		assert(deq(det({-1,5,2},{8,2,3},{-1,-3,7}),-362));

		//testperpline();

		assert(deq(distlinept({1,2,3},{-1,0,1},{-1,-6,-47}),6*sqrt(38)));
		
		assert(deq(distfinitelinept({1,2,3},{-1,0,1},{-1,-6,-47}),6*sqrt(65)));
		assert(deq(distfinitelinept({1,2,3},{-1,0,1},{18,13,-28}),6*sqrt(38)));
		assert(deq(distfinitelinept({1,2,3},{-1,0,1},{39,34,-7}),2*sqrt(642)));

		testrotations();

		testconvexhull();
	}
}

namespace SPOJ_CH3D {
	void solve() {
		//srand(time(0));
		int N,T;
		scanf("%d",&T);
		cout << fixed << setprecision(2);
		vector<Pt> p;
		while (T--) {
			scanf("%d",&N);
			p.resize(N);
			for (Pt &i:p) scanf("%lf %lf %lf",&i.x,&i.y,&i.z);
			//random_shuffle(p.begin(),p.end());
			Hull::convexhull(p);
			double a=0,v=0;
			for (int i=0;i<Hull::f;i++) {
				using namespace Hull;
				v+=det(faces[i].a,faces[i].b,faces[i].c);
				a+=abs(cross(faces[i].b-faces[i].a,faces[i].c-faces[i].a));

				/*for (int j=0;j<3;j++)
					cerr << faces[i][j].x << ' ' << faces[i][j].y << ' ' << faces[i][j].z << endl;
				cerr << endl;*/
			}
			v/=6;
			a/=2;
			printf("%.4f %.4f\n",a,v);
			//cout << a << ' ' << v << '\n';
		}
	}
}



int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);

	//Test::solve();
	SPOJ_CH3D::solve();
}
