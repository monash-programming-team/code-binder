#include <bits/stdc++.h>

#define x first
#define y second

using namespace std;

typedef long long ll;
typedef pair<int,int> ii;
typedef pair<ll,ll> pll;

// Basic Maths Functions
//
// Author      : Peter Whalan
// Date        : June 6, 2016
// Reliability : 1
// Tested On   :
//	CF:666/C (pmod, gcd, inv)

typedef __int128 big;
ll mult(ll a,ll b,ll mod) {
	return big(a)*b%mod;
}


//b^e mod mod. O(log e)
ll pmod(ll b,ll e,ll mod) {
	if (e==0) return 1%mod;
	ll r=pmod(b,e/2,mod);
	r=mult(r,r,mod);//(ll)r*r%mod;
	if (e%2) r=mult(r,b,mod);//(ll)r*b%mod;
	return r;
}

//Extended Euclidean Algorithm. O(log(min(m,n)))
//All operations don't overflow
//abs(g) is gcd (Note that it could be negative)
//am+bn=g with abs(a)<abs(n)/g and and abs(b)<abs(m)/g except when m or n is 0
//in which case a and b are 0 and 1 (in some order).
//
//The proof follows easily from
// |b| <= |a'|+|m/n||b'| <= |r|+|m/n||n| (inductive step when m,n coprime)
// == |r+m/n*n| (sign r == sign m for truncated division)
// == |m|

int a,b,g;
void gcd(int m,int n) {
    if (n==0) {a=1;b=0;g=m;return;}
    gcd(n,m%n);
    int t=a-b*(m/n);
    a=b;b=t;
}

//Multiplicative inverse of a mod mod. O(log(min(a,mod)))
//The remarks above show that |inv(a,mod)|<|mod|/|gcd(a,mod)| but inv(a,mod)
//can be positive or negative even when both a and mod are positive.
//If a and mod aren't coprime then the inverse doesn't exist and in general we
//have a*inv(a,mod) = abs(gcd(a,mod)) (mod mod)

int inv(int a,int mod) {
	gcd(mod,a);
	return g>0?b:-b;
}

//Discrete Logarithm. O(sqrt(M)log(M))
//Solves a^x == b (mod mod) for integer x. The smallest nonnegative x is chosen.
//Returns -1 if there is no such x. x is assumed to be strictly less than M;  
//To optimise set M to
//	phi(mod/gcd(mod,lcm(a^lg(mod),b^lg(mod)))) + lg(mod)
//	lg(mod) is the maximum multiplicity of a prime factor. This is length of
//	path before entering cycle.
//Can also derive extra conditions to determine when solution exists before
//running algorithm.

pll seen[5000000]; // Must be at least ceil(sqrt(M))
ll discretelog(ll a,ll b,ll mod) {
	ll M=mod;
	ll s=0,as=1,bas; // step size, a^s, ba^s
	for (;s*s<M;s++) {
		as=mult(as,a,mod);
		bas=mult(b,as,mod);//(ll)b*as%mod;
		seen[s]={bas,s+1};
	}
	sort(seen,seen+s);

	for (ll i=1,ap=1,ct=0,p;i<=s && ct<=s;i++) {
		ap=mult(ap,as,mod);//(ll)ap*as%mod;
		int j=lower_bound(seen,seen+s,pll{ap+1,0})-seen;
		for (;--j>=0 && seen[j].x==ap && ct<=s;ct++)
			if (pmod(a,p=(ll)i*s-seen[j].y,mod)==b) return p;
	}
	return -1;
}

namespace Test {
	void test_discretelog() {
		int maxmod=250;
		for (int mod=1;mod<=maxmod;mod++) for (int a=0;a<mod;a++) for (int b=0;b<mod;b++) {
			bool f=0;
			int p=1%mod;
			for (int ans=0;ans<mod;ans++,p=p*a%mod) if (p==b) {
				f=1;
				if (discretelog(a,b,mod)!=ans) {
					cerr << a << ' ' << b << ' ' << mod << ' ' << discretelog(a,b,mod) << ' ' << ans << endl;
					return;
				}
				break;
			}
			if (!f && discretelog(a,b,mod)!=-1) {
					cerr << a << ' ' << b << ' ' << mod << ' ' << discretelog(a,b,mod) << endl;
					//return;
			}
		}
		cerr << "success\n";
	}

	void solve() {
		cout << pmod(3,2,5) << endl;

		gcd(9,6);
		printf("%d %d %d\n",a,b,g);

		cout << inv(2,5) << endl;

		//Very general test of inv
		for (int i=-10000;i<=10000;i++) for (int mod=-10000;mod<=10000;mod++) {
			if (mod!=0) {
				int p=i*inv(i,mod)%mod;
				int gmod=abs(g)%mod;
				if (p<0) p+=abs(mod); //mod could be negative
				if (p!=gmod) {
					printf("%d %d %d\n",i,mod,inv(i,mod));
					assert(p==gmod);
				}
			}
		}
	}
};

namespace SPOJ_MOD {
	typedef long long ll;

	void solve() {
		int x,z,k;
		while (cin >> x >> z >> k && x) {
			int r=discretelog(x%z,k%z,z);
			if (r==-1) cout << "No Solution\n";
			else cout << r << '\n';
		}
	}
}

//https://www.codechef.com/problems/DISLOG
//one test case and server doesn't have __int128
namespace Codechef_DISLOG {
	void solve() {
		assert(305811030771 == discretelog(21309,696969,999998999999));
	}
}

// Don't seem to be able to submit to this one
namespace SPOJ_KRYPT_67_DLOG {
	void solve() {
		int p,g,y;
		cin >> p >> g >> y;
		cout << discretelog(g,y,p) << endl;
	}
}

namespace UVA_10225 {
	void solve() {
		int p,b,n;
		while (cin >> p >> b >> n) {
			ll ans=discretelog(b,n,p);
			if (ans==-1) cout << "no solution\n";
			else cout << ans << '\n';
		}
	}
}

namespace ACM_ICPC_LIVE_ARCHIVE_7457 {
	void solve() {
		int p,a,b;
		cin >> p;
		while (cin >> a >> b && a) {
			int ans=discretelog(a,b,p);
			if (ans==-1) ans++;
			cout << ans << '\n';
		}
	}
}


int main() {
	ios::sync_with_stdio(0);
	cin.tie(0);
	//Test:solve();
	//Test::test_discretelog();
	//SPOJ_MOD::solve();
	//Codechef_DISLOG::solve();
	//SPOJ_KRYPT_67_DLOG::solve();
	//UVA_10225::solve();
	ACM_ICPC_LIVE_ARCHIVE_7457::solve();
}
