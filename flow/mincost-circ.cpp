// Problem: ACM ICPC World Finals 2011, Problem D 'Chips'
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
#define TASKNAME "chips"

#ifdef LOCAL
#define eprintf(...) fprintf(stderr,__VA_ARGS__)
#else
#define eprintf(...)
#endif

#define TIMESTAMP(msg) eprintf("[" msg "] Time = %.3lfs\n",\
                             clock() * 1.0 /CLOCKS_PER_SEC)
#define TIMESTAMPF(msg, ...) eprintf("[" msg "] Time = %.3lfs\n",\
                      __VA_ARGS__, clock() * 1.0/ CLOCKS_PER_SEC)

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
typedef pair<ll, ll> pll;
typedef pair<ld, ld> pld;

const int inf = 1e9;
const double eps = 1e-9;
const double INF = inf;
const double EPS = eps;

// BEGIN ALGO
class Solver2 { // Min-cost circulation
  struct Edge {
    int to, ne, w, c;
  };
  vector<Edge> es;
  vi firs;
  int curRes;

  public:
  Solver2(int n) : es(), firs(n, -1), curRes(0) {}
  // from, to, capacity (max.flow), cost
  int adde(int a, int b, int w, int c) {
    Edge e;
    e.to = b; e.ne = firs[a];
    e.w = w; e.c = c;
    es.pb(e);
    firs[a] = sz(es) - 1;
    
    e.to = a; e.ne = firs[b];
    e.w = 0; e.c = -c;
    es.pb(e);
    firs[b] = sz(es) - 1;
    return sz(es) - 2;
  }
  // increase capacity of edge 'id' by 'w'
  void ince(int id, int w) {
    es[id].w += w;
  }
  int solve() {
    const int n = sz(firs);
    
    for (;;) {
      vi d(n, 0), fre(n, -1);
      vi chd(n, -1);
      
      int base = -1;
      for (int step = 0; step < n; step++) { /*BOXNEXT*/
        for (int i = 0; i < sz(es); i++) if (es[i].w > 0) {
          int b = es[i].to;
          int a = es[i ^ 1].to; /*BOXNEXT*/
          if (d[b] <= d[a] + es[i].c) continue;
          d[b] = d[a] + es[i].c;
          fre[b] = i;
          if (step == n - 1)
            base = b;
        }
      }
      if (base < 0) break;
      
      vi seq;
      vb was(n, false);/*BOXNEXT*/
      for (int x = base;; x = es[fre[x] ^ 1].to) {
        if (!was[x]) {
          seq.pb(x);
          was[x] = true;
        } else {
          seq.erase(
              seq.begin(),
              find(seq.begin(), seq.end(),
              x
          ));
          break;
        }
      }
      for (int i = 0; i < sz(seq); i++) {
        int v = seq[i];
        int eid = fre[v];
        assert(es[eid].w > 0);
        es[eid].w--;
        es[eid ^ 1].w++;
        curRes += es[eid].c;
      }
    }
    return curRes;
  }
};
// END ALGO

class Solver1 { // Min-cost LR circulation
  static const int BASEC = 10000;

  Solver2 flow;
  int st, en;
  int minFlow;

  public:
  Solver1(int n) : flow(n + 2), st(n), en(n + 1), minFlow(0) {
    flow.adde(en, st, inf, -BASEC);
  }
  int adde(int a, int b, int minf, int maxf, int c) {
    assert(minf <= maxf);
    assert(c >= -1 && c <= 1);
    if (minf > 0) {
      flow.adde(st, b, minf, c);
      flow.adde(a, en, minf, 0);
      minFlow += minf;
    }
    return flow.adde(a, b, maxf - minf, c);
  }
  void ince(int id, int toadd) {
    flow.ince(id, toadd);
  }
  int solve() {
    int res = flow.solve();
//    eprintf("res=%d\n", res);
    assert(res >= (minFlow + 1) * -BASEC);
    if (res > minFlow * -BASEC) return inf;
    return res + minFlow * BASEC;
  }
};

const int MAXN = 40;
char f[MAXN][MAXN + 1];

int main(){
  #ifdef LOCAL
  freopen(TASKNAME".in","r",stdin);
  freopen(TASKNAME".out","w",stdout);
  #endif

  int n, a, b;
  int TN = 0;
  while (scanf("%d%d%d", &n, &a, &b) == 3) {
    if (!n && !a && !b) break;

    int minCnt = 0, maxCnt = 0;
    for (int i = 0; i < n; i++) {
      scanf("%s", f[i]);
      for (int j = 0; j < n; j++) {
        if (f[i][j] == 'C')
          minCnt++;
        if (f[i][j] != '/')
          maxCnt++;
      }
    }

    TN++;
    printf("Case %d: ", TN);

    Solver1 s(2 * n);
    for (int y = 0; y < n; y++)
    for (int x = 0; x < n; x++) if (f[y][x] != '/') {
      s.adde(2 * y + 1, 2 * x, (f[y][x] == 'C') ? 1 : 0, 1, 0);
    }
    vi es(n);
    for (int i = 0; i < n; i++)
      es[i] = s.adde(2 * i, 2 * i + 1, 0, 0, -1);

    int pcnt = 0;
    int best = -1;
    for (int ans = minCnt; ans <= maxCnt; ans++) {
      int maxInRow = a * ans / b;

      if (pcnt < maxInRow) {
        for (int i = 0; i < n; i++) {
          s.ince(es[i], maxInRow - pcnt);
//          s.adde(2 * i, 2 * i + 1, 0, maxInRow - pcnt, -1);
        }
        pcnt = maxInRow;
      }
      int res = -s.solve();
      if (res == ans)
        best = ans;
    }
    if (best < 0) printf("impossible\n");
    else printf("%d\n", best - minCnt);
  }

  TIMESTAMP("end");
  return 0;
}
