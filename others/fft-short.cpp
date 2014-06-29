//#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <cstring>
#include <string>
#include <cmath>
#include <cassert>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <list>
#include <queue>
#include <deque>
#include <stack>
#include <cstdlib>
#include <cstdio>
#include <iterator>
#include <functional>
#include <bitset>
#define mp make_pair
#define pb push_back

//#define DEBUG

#ifdef LOCAL
#define eprintf(...) fprintf(stderr,__VA_ARGS__)
#else
#define eprintf(...)
#endif

#define TIMESTAMP(x) eprintf("["#x"] Time : %.3lf s.\n", clock()*1.0/CLOCKS_PER_SEC)
#define TIMESTAMPf(x,...) eprintf("[" x "] Time : %.3lf s.\n", __VA_ARGS__, clock()*1.0/CLOCKS_PER_SEC)

#if ( ( _WIN32 || __WIN32__ ) && __cplusplus < 201103L)
  #define LLD "%I64d"
#else
  #define LLD "%lld"
#endif

using namespace std;

#define TASKNAME "F"

typedef vector<int> vi;
#define sz(a) ((int)(a).size())

#ifdef LOCAL
static struct __timestamper {
  string what;
  __timestamper(const char* what) : what(what){};
  __timestamper(const string& what) : what(what){};
  ~__timestamper(){
    TIMESTAMPf("%s", what.data());
  }
} __TIMESTAMPER("end");
#else
struct __timestamper {};
#endif

typedef long long ll;
typedef long double ld;


const int MOD = 998244353;
const int ROOT = 3;

int mmul(int a, int b) {
	return (a * 1LL * b) % MOD;
}
void madd(int& a, int b) {
    if ((a += b) >= MOD) a -= MOD;
}

int mpow(int a, int b){
  if (!b) return 1;
  if (b & 1) return (mpow(a, b-1) * 1LL * a) % MOD;
  int temp = mpow(a, b/2);
  return (temp * 1LL * temp) % MOD;
}

bool checkRoot(int x){
  return mpow(x, (MOD - 1) / 7) != 1 && mpow(x, (MOD - 1) / 17) != 1 && mpow(x, (MOD - 1) / 2) != 1;
}

// BEGIN ALGO
const int MROOT = 19; 
const int MROOTP = 1<<MROOT;
int rts[MROOTP + 10], brev[MROOTP + 10];
// \emph{Don't forget to call before}
void PreCalcRoots(){
  rts[0] = 1; // ROOT is primary root for MOD
  rts[1] = mpow(ROOT, (MOD-1) / MROOTP);
  for (int i = 2; i < MROOTP; i++)
    rts[i] = mmul(rts[i-1], rts[1]);
  for (int i = 0; i < MROOTP; i++) /*BOXNEXT*/
    brev[i] = (brev[i>>1]>>1) | ((i&1) << (MROOT-1));
}
inline void butterfly(int &a, int &b, int x){
  int temp = mmul(x, b); b = a;
  madd(a, temp); madd(b, MOD - temp);
}
void fft(vi &a, bool inv){
  int n = __builtin_ctz(sz(a));
  for (int i = 0; i < (1<<n); i++){
    int temp = brev[i] >> (MROOT - n);
    if (temp > i) swap(a[i], a[temp]);
  }
  for (int step = 0; step < n; step++){
    int pos = 0; /*BOXNEXT*/
    int dlt = (inv ? -1 : 1) * (1 << (MROOT - step - 1));
    for (int i = 0; i < (1<<n); i++){
      if ((i ^ (1<<step)) > i) /*BOXNEXT*/
        butterfly(a[i], a[i ^ (1<<step)], rts[pos]);
      pos = (pos + dlt) & (MROOTP-1);
    }
  }
}
vi multiply(vi a, vi b){
  int rsz = sz(a) + sz(b), rsz2 = 1;
  while (rsz2 < rsz) rsz2 *= 2;
  a.resize(rsz2); b.resize(rsz2);
  fft(a, false); fft(b, false);
  for (int i = 0; i < sz(a); i++)
    a[i] = mmul(a[i], b[i]);
  fft(a, true);
  int in = mpow(sz(a), MOD - 2);
  for (int i = 0; i < sz(a); i++)
    a[i] = mmul(a[i], in);
  return a;
}
vi inverse(vi a){
  assert(a[0] != 0);
  vi x(1, mpow(a[0], MOD - 2));
  while (sz(x) < sz(a)) { /*BOXNEXT*/
    vi temp(a.begin(), a.begin() + min(sz(a), 2*sz(x)));
    vi nx = multiply(multiply(x, x), temp);
    x.resize(2*sz(x));
    for (int i = 0; i < sz(x); i++) /*BOXNEXT*/
      madd(x[i], x[i]), madd(x[i], MOD - nx[i]);
  }
  return x;
}
// END ALGO

const int MAXN = 257000;

int facs[MAXN];
int ifacs[MAXN];

void PreCalcFacs(){
  facs[0] = ifacs[0] = 1;
  for (int i = 1; i < MAXN; i++){
    facs[i] = (facs[i-1] * 1LL * i) % MOD;
    ifacs[i] = mpow(facs[i], MOD - 2);
  }
}

int cnk(int n, int k){
  if (n < 0 || k < 0 || k > n) return false;
  int res = facs[n];
  res = (res * 1LL * ifacs[k]) % MOD;
  res = (res * 1LL * ifacs[n-k]) % MOD;
  return res;
}

ll get2(int x1, int x2, int t){
  int dlt = t - abs(x1 - x2);
  if (dlt < 0 || dlt % 2) return 0;
  return cnk(t, dlt / 2);
}

ll get1(int x1, int x2, int y1, int y2, int t){
  return (get2(x1 + y1, x2 + y2, t) * 1LL * get2(x1 - y1, x2 - y2, t)) % MOD;
}

int get(int x1, int x2, int y1, int y2, int t){
  return (get1(x1, x2, y1, y2, t) -
      get1(x1,-x2, y1, y2, t) -
      get1(x1, x2, y1,-y2, t) +
      get1(x1,-x2, y1,-y2, t) +
      3LL * MOD) % MOD;
}

int main(){
  assert(checkRoot(ROOT));
  PreCalcRoots();
  PreCalcFacs();
  #ifdef LOCAL
  assert(freopen(TASKNAME".in","r",stdin));
  assert(freopen(TASKNAME".out","w",stdout));
  #endif

  int x1, y1, x2, y2, t;
  scanf("%d %d %d %d %d",&x1,&y1,&x2,&y2, &t);

  vi res(t+1);
  vi res0(t+1);

  for (int i = 0; i <= t; i++){
    res[i] = get(x1, x2, y1, y2, i);
    res0[i] = get(x2, x2, y2, y2, i);
  }

  res = multiply(res, inverse(res0));

  printf("%d\n", res[t]);

  return 0;
}