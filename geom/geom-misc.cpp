// Problem: ACM ICPC 2008 WF, problem H (Painter)
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

#define pb push_back
#define mp make_pair
#define TASKNAME "painter"

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

using namespace std;

typedef long double ld;
typedef long long ll;
typedef vector<ll> vll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<bool> vb;
typedef vector<vb> vvb;
typedef pair<int, int> pii;

const int inf = 1e9;
const double eps = 1e-9;
const double INF = inf;
const double EPS = eps;

int sgn(ll x) { return x < 0 ? -1 : x > 0; }
struct pt {
  int x, y;
  pt() : x(0), y(0) {}
  pt(int x, int y) : x(x), y(y) {}
  inline pt operator-(const pt &p2) const { return pt(x - p2.x, y - p2.y); }
  inline int operator*(const pt &p2) const {
    return sgn(ll(x) * p2.y - ll(y) * p2.x);
  }
  inline bool operator<(const pt &p2) const {
    if (x != p2.x) return x < p2.x;
    return y < p2.y;
  }
  inline bool operator==(const pt &p2) const {
    return x == p2.x && y == p2.y;
  }
};

struct segm {
  pt a, b;
  int id;

  segm() : id(-1) {}
  segm(const pt &p1, const pt &p2, int _id) : a(min(p1, p2)), b(max(p1, p2)), id(_id) {}
  inline bool operator==(const segm &s2) const {
    return a == s2.a && b == s2.b && id == s2.id;
  }
  void print() const {
    eprintf("(%d %d)-->(%d %d)\n", a.x, a.y, b.x, b.y);
  }
};

//BEGIN ALGO
//abs(coordinates)<=1e6

// next two checks intersection of [a1,b1], [a2,b2]
inline bool if_cross(int a1, int b1,
                     int a2, int b2) {
  if (a1 > b1) swap(a1, b1);
  if (a2 > b2) swap(a2, b2);
  return max(a1, a2) <= min(b1, b2);
}

inline bool if_cross(pt a1, pt b1,
                     pt a2, pt b2) {
  if (!if_cross(a1.x, b1.x, a2.x, b2.x))
     return false;
  if (!if_cross(a1.y, b1.y, a2.y, b2.y))
     return false;

  if ( ((a2 - a1) * (b1 - a1)) *
       ((b2 - a1) * (b1 - a1)) > 0)
     return false;
  swap(a1, a2);
  swap(b1, b2);
  if ( ((a2 - a1) * (b1 - a1)) *
       ((b2 - a1) * (b1 - a1)) > 0)
     return false;
  return true;
}

//return pair<num,den> : line(a,b) cross x=const with y=num/den;
inline pair<ll, ll> get_intersec(const pt &a,
                                 const pt &b,
                                 int x) {
  if (a.x == b.x)
    return pair<ll, ll>(ll(a.y), 1LL);
  return pair<ll, ll>(
      ll(a.y) * ll(b.x - a.x) +
      ll(b.y - a.y) * ll(x - a.x),
    ll(b.x - a.x)
  );
}

//struct segm { pt a, b; int id; }
//id = id of polygon, NOT SEGMENT
//assumes that only segments from the same simple polygon
//can intersect each other
inline bool operator<(const segm &s1,
                      const segm &s2) {
  if (s1 == s2) return false;
  if (s1.id != s2.id &&
      if_cross(s1.a, s1.b, s2.a, s2.b))
    throw mp(s1.id, s2.id);

  int xl = max(s1.a.x, s2.a.x);
  int xr = min(s1.b.x, s2.b.x);
  int xmid = (xl + xr) / 2;
  /*BOXNEXT*/
  pair<ll, ll> y1 = get_intersec(s1.a, s1.b, xmid); /*BOXNEXT*/
  pair<ll, ll> y2 = get_intersec(s2.a, s2.b, xmid);
  int cmp = sgn(y1.first * y2.second -
                y1.second * y2.first);
  if (cmp != 0) return cmp < 0;

  return s1.a.x + s1.b.x  < s2.a.x + s2.b.x;
}

// END ALGO

struct Ev {
  int x, ty;
  segm s;

  bool operator<(const Ev &e2) const {
    if (x != e2.x) return x < e2.x;
    return ty > e2.ty;
  }
};

vi pars, deep;
int calc(int v) {
  if (deep[v] >= 0) return deep[v];
  if (pars[v] == v) {
    return deep[v] = 1;
  }
  assert(pars[v] >= 0);
  return deep[v] = calc(pars[v]) + 1;
}

int main(){
  freopen(TASKNAME".in","r",stdin);
  #ifdef LOCAL
  freopen(TASKNAME".out","w",stdout);
  #endif

  int n;
  int TN = 0;
  while (scanf("%d", &n) == 1) {
    if (n < 0) break;

    vector<segm> ss;
    for (int i = 0; i < n; i++) {
      pt a, b, c;
      scanf("%d%d%d%d%d%d", &a.x, &a.y, &b.x, &b.y, &c.x, &c.y);
      a.x *= 2;
      a.y *= 2;
      b.x *= 2;
      b.y *= 2;
      c.x *= 2;
      c.y *= 2;
      ss.pb(segm(a, b, i));
      ss.pb(segm(a, c, i));
      ss.pb(segm(b, c, i));
    }

    vector<Ev> evs;
    for (int i = 0; i < sz(ss); i++) {
      Ev e;
      e.x = ss[i].a.x; e.ty = 1; e.s = ss[i];
      evs.pb(e);
      e.x = ss[i].b.x; e.ty = -1;
      evs.pb(e);
    }
    sort(evs.begin(), evs.end());

    TN++;
    printf("Case %d: ", TN);
    try {
      const int INFCOORD = 1e6;
      set<segm> curs;
      curs.insert(segm(pt(-INFCOORD, -INFCOORD), pt(INFCOORD, -INFCOORD), n));
      curs.insert(segm(pt(-INFCOORD,  INFCOORD), pt(INFCOORD,  INFCOORD), n));
      pars = vi(n + 1, -1);
      pars[n] = n;

      for (int i = 0; i < sz(evs); i++) {
        set<segm>::iterator it = curs.lower_bound(evs[i].s);
        assert(it != curs.begin());
        assert(it != curs.end());
        int a = it->id;

        it--;
        int b = it->id;
        int me = evs[i].s.id;

        assert(pars[a] >= 0);
        assert(pars[b] >= 0);

        if (pars[me] == -1) {
          assert(a != me && b != me);
          if (a == b) {
            pars[me] = a;
          } else if (pars[a] == pars[b]) {
            pars[me] = pars[a];
          } else if (a == pars[b]) {
            pars[me] = a;
          } else if (b == pars[a]) {
            pars[me] = b;
          } else {
            assert(false);
          }
        }


        if (evs[i].ty == 1) {
          curs.insert(evs[i].s);
        } else {
          assert(curs.erase(evs[i].s));
        }
      }
      assert(sz(curs) == 2);

      deep = vi(sz(pars), -1);
      int ans = 0;
      for (int i = 0; i < sz(pars); i++) {
        ans = max(ans, calc(i));
      }

      printf("%d shades\n", ans);
      eprintf("nya %d\n", ans);
    } catch (const pii &e) {
      printf("ERROR\n");
      eprintf("botva: %d %d\n", e.first, e.second);
    }
  }
  TIMESTAMP(end);
  return 0;
}

