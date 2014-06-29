#include <cstdio>
#include <cstring>
#include <cassert>
#include <ctime>
#include <algorithm>
#include <vector>

using namespace std;

#define pb push_back
#define mp make_pair
#define sz(x) ((int)(x).size())
#define eprintf(...) fprintf(stderr, __VA_ARGS__)
#define TIMESTAMP(x) eprintf("[" x "] Time=%.3lfs\n", clock() * 1.0 / CLOCKS_PER_SEC)

typedef long long ll;
typedef vector<int> vi;
typedef pair<ll, int> pli;

#ifdef _WIN32
#define LLD "%I64d"
#elif linux
#define LLD "%lld"
#endif

// BEGIN ALGO
// mod should be not greater than 4e18
inline ll mmul(ll a, ll b, ll mod) {
  assert(0 <= a && a < mod);
  assert(0 <= b && b < mod);
  if (mod < int(1e9)) return a * b % mod;
  ll k = (ll)((long double)a * b / mod);
  ll res = a * b - k * mod;
  res %= mod;
  if (res < 0) res += mod;
  return res;
}
inline ll mpow(ll a, ll b, ll mod) {
  if (mod==1)
    return 0;
  ll res = 1;
  for (; b; b >>= 1, a = mmul(a, a, mod))
    if (b & 1) res = mmul(res, a, mod);
  return res;
}
inline ll gcd(ll a, ll b) {
  return b ? gcd(b, a % b) : a;
}

inline int randint(int l, int r) {
  static unsigned int seed = 239017u;
  seed=seed*1664525u+1013904223u;
  return l+seed%(r-l+1);
}

vector<pli> res;

//run \emph{only} if n%as[i]!=0 for all i
inline bool run_miller_rubbin(ll n) {
  assert(n > 27 && (n & 1));
  ll t = n - 1; int s = 0;
  while (!(t & 1)) {
    t >>= 1;
    s++;
  }
  const int as[] = { 2, 3, 5, 7, 11, 13,
                     17, 19, 23, 27, -1 };
  for (int i = 0; as[i] >= 0; i++) {
    ll cur = mpow(as[i], t, n);
    if (cur==1)
      continue;
    bool good = false;
    for (int i2 = 0; i2 < s; i2++) {
      if (cur == n - 1) {
        good = true;
        break;
      }
      assert(0 <= cur && cur < n);
      cur = mmul(cur, cur, n);
    }
    if (!good) return false;
  }
  return true;
}

const int MAXD = 110;

//\emph{only} for n>MAXD
inline void pollard(ll n) {
  if (n == 1) return;
  if (run_miller_rubbin(n)) {
    res.pb(mp(n, 1));
    return;
  }
  assert(n > MAXD);
  
  const ll MAXX = int(1e9);
  
  int sumIt = 0;
  for (int i = 0; i<3 || sumIt<1000; i++) {
    int q = randint(2, MAXD);
    ll x = randint(1, min(MAXX, n - 1));
    ll y = x, g = 0;
    int maxJ = 1 << (i + 18);
    for (int j = 1; j < maxJ; j++) {
      sumIt++;
      x = mmul(x, x, n);
      if ((x += q) >= n) x -= n;
      ll z = abs(x - y);
      g = gcd(z, n);
      if (g != 1) break;
      if ((j & (j - 1)) == 0) {
        y = x;
      }
    }
    if (1 < g && g < n) {
      pollard(g);
      pollard(n / g);
      return;
    }
  }
  eprintf("Failed at n="LLD"\n", n);
  assert(false);
}

vi ps; // all primes <= MAXD

//big primes can be not grouped
void factorize(ll n) {
  res.clear();
  for (int i = 0; i < sz(ps); i++) {
    int deg = 0;
    while (n % ps[i] == 0)
      n /= ps[i], deg++;
    if (deg)
      res.pb(mp(ll(ps[i]), deg));
  }
  if (n > 1) {
    assert(n > MAXD);
    pollard(n);
  }
}
// END ALGO

bool isp[MAXD + 1];
int main() {
  #ifdef DEBUG
  freopen(".in", "r", stdin);
  freopen(".out", "w", stdout);
  #endif

  //eprintf("%d\n", run_miller_rubbin(5394826801ll));
  memset(isp, 1, sizeof isp);
  for (int i = 2; i <= MAXD; i++) if (isp[i]) {
    ps.pb(i);
    for (int i2 = i * i; i2 <= MAXD; i2 += i)
      isp[i2] = false;
  }

  ll n;
  while (scanf(LLD, &n) == 1) {
    if (!n) break;

    factorize(n);
    sort(res.begin(), res.end());
    int ptr = 0;
    for (int i = 0; i < sz(res); i++)
      if (ptr == 0 || res[ptr - 1].first != res[i].first) {
        res[ptr++] = res[i];
      } else {
        res[ptr - 1].second += res[i].second;
      }
    res.resize(ptr);
    for (int i = 0; i < sz(res); i++)
      printf(LLD"^%d%c", res[i].first, res[i].second, "\n "[i + 1 < sz(res)]);
  }
  TIMESTAMP("end");
  return 0;
}
