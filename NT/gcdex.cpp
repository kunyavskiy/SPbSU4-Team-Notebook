#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <string>
#include <vector>
#include <deque>
#include <queue>
#include <list>
#include <set>
#include <map>
#include <sstream>
#include <iostream>
#include <testlib.h>

#define pb push_back
#define mp make_pair
#define TASKNAME ""

#ifdef LOCAL
#define eprintf(...) fprintf(stderr,__VA_ARGS__)
#else
#define eprintf(...)
#endif

#define TIMESTAMP(x) eprintf("[" #x "] Time = %.3lfs\n",clock()*1.0/CLOCKS_PER_SEC)

#ifdef linux
#define LLD "%lld"
#else
#define LLD "%I64d"
#endif

#define sz(x) ((int)(x).size())
#define strstr strstr_wregthrtu

using namespace std;

typedef long double ld;
typedef long long ll;
typedef vector<ll> vll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<bool> vb;
typedef vector<vb> vvb;
typedef pair<int, int> pii;
typedef pair<ll,ll> pll;

//const int inf = 1e9;
const double eps = 1e-9;
//const double INF = inf;
const double EPS = eps;

const long long checkItCnt = 1000000;
const long long maxn = (ll)1e18;

// BEGIN ALGO

//a*x+b*y=__gcd(a,b)  a>0 b>0

void gcdex2(ll a,ll b,ll *x,ll *y) {
  if (a<b)
    swap(a,b), swap(x,y);
  ll a1=1, a2=0, b1=0, b2=1, tmp;
  while (b)
  {
    tmp=a/b, a%=b, a1-=tmp*b1, a2-=tmp*b2;
    swap(a,b), swap(a1,b1), swap(a2,b2);
  }
  *x=a1, *y=a2;
}

void gcdex(ll a,ll b,ll &x,ll &y) {
  if (!a){
    x = 0, y = 1;
    return;
  }
  ll x1, y1;
  gcdex(b % a, a, x1, y1);
  x = y1 - (b / a) * x1;
  y = x1;
}

// END ALGO

int main(){

  /*#ifdef LOCAL
  freopen(TASKNAME".in","r",stdin);
  freopen(TASKNAME".out","w",stdout);
  #endif*/

  rnd.setSeed(239);

  for (ll i=0; i<checkItCnt; i++)
  {
    ll a = rnd.next(1ll,maxn);
    ll b = rnd.next(1ll,maxn);
    ll x, y, x2, y2;
    gcdex(a,b,x,y);
    gcdex2(a,b,&x2,&y2);

    assert(abs(x2)<b && abs(y2)<a);
    assert(a*x+b*y==__gcd(a,b));
    assert(x==x2 && y==y2);
  }

  TIMESTAMP(end);
  return 0;
}

