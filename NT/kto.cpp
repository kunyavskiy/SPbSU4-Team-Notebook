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

const long long checkItCnt = 2000000;
const long long maxMod = (ll)2e9;
const long long maxT = 1000;

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

// BEGIN ALGO
//return pair(nmod,nr)
//nr%mod1=r1, nr%mod2=r2
//nmod=mod1*mod2/gcd(mod1,mod2)
//if input incosistent return mp(-1,-1)

pll kto (ll mod1, ll r1, ll mod2, ll r2)
{
  ll d=__gcd(mod1,mod2);
  if (r1%d!=r2%d)
    return mp(-1,-1);
  ll rd=r1%d;
  mod1/=d, mod2/=d, r1/=d, r2/=d;
  
  if (mod1<mod2)
    swap(mod1,mod2), swap(r1,r2);
  
  ll k=(r2-r1)%mod2;
  if (k<0)
    k+=mod2;
  
  ll x, y;
  gcdex(mod1,mod2,x,y);
  x%=mod2;
  if (x<0)
    x+=mod2;
  k*=x, k%=mod2;
  return mp(mod1*mod2*d,(k*mod1+r1)*d+rd);
}
// END ALGO

int main(){

  /*#ifdef LOCAL
  freopen(TASKNAME".in","r",stdin);
  freopen(TASKNAME".out","w",stdout);
  #endif*/

  rnd.setSeed(239);

  int cnt1=0, cnt2=0;

  for (ll i=0; i<checkItCnt; i++)
  {
    ll mod1 = rnd.next(1ll,maxMod);
    ll mod2 = rnd.next(1ll,maxMod);
    ll r1 = rnd.next(0ll,mod1-1);
    ll r2 = rnd.next(0ll,mod2-1);
    pll ret = kto(mod1,r1,mod2,r2);

    ll d=__gcd(mod1,mod2);
    ll dr1=r1%d;
    ll dr2=r2%d;

    if (ret.first==-1)
    {
      assert(dr1!=dr2);
      cnt1++;
      continue;
    }

    assert(ret.first==(mod1/d)*mod2);
    assert(ret.second%mod1==r1 && ret.second%mod2==r2);
    cnt2++;

  }

  for (ll i=0; i<checkItCnt; i++)
  {

    ll t1=rnd.next(1ll,maxT);

    ll mod1 = rnd.next(1ll,maxMod/t1);
    ll mod2 = rnd.next(1ll,(maxMod*maxMod/mod1)/t1);

    ll t=rnd.next(1ll,t1);
    mod1*=t, mod2*=t;

    ll r1 = rnd.next(0ll,mod1-1);
    ll r2 = rnd.next(0ll,mod2-1);
    pll ret = kto(mod1,r1,mod2,r2);

    ll d=__gcd(mod1,mod2);
    ll dr1=r1%d;
    ll dr2=r2%d;

    if (ret.first==-1)
    {
      assert(dr1!=dr2);
      cnt1++;
      continue;
    }

    assert(ret.first==(mod1/d)*mod2);
    assert(ret.second%mod1==r1 && ret.second%mod2==r2);
    cnt2++;

  }

  cerr<<cnt1<<" "<<cnt2<<endl;
  TIMESTAMP(end);
  return 0;
}
