// Problem: http://acm.timus.ru/problem.aspx?space=1&num=1016
#include <cstdio>
#include <cassert>
#include <cstring>
#include <algorithm>
#include <queue>
#include <vector>

using namespace std;

#define pb push_back
#define mp make_pair
#define sz(x) ((int)(x).size())

typedef vector<int> vi;
typedef vector<vi> vvi;
typedef pair<int, int> pii;

#ifdef DEBUG
#define eprintf(...) fprintf(stderr, __VA_ARGS__)
#else
#define eprintf
#endif

class Solver {
  vvi es, esw;
  vi d, fr;

  public:
  Solver(int n) : es(n), esw(n) {}
  void adde(int a, int b, int w) {
    es[a].pb(b);
    esw[a].pb(w);
  }
  void solve(int st) {
    const int n = sz(es);
    d = vi(n, 1e9);
    fr = vi(n, -1);

    priority_queue<pii> q;
    d[st] = 0;
    q.push(mp(-d[st], st));
    while (!q.empty()) {
      int cd = -q.top().first;
      int v = q.top().second;
      q.pop();
      if (d[v] != cd) continue;

      for (int i = 0; i < sz(es[v]); i++) {
        int b = es[v][i];
        int nd = cd + esw[v][i];
        if (nd >= d[b]) continue;
        d[b] = nd;
        q.push(mp(-d[b], b));
        fr[b] = v;
      }
    }
  }
  int get(int v) { return d[v]; }
  vi getP(int v) {
    vi res;
    while (v >= 0) {
      res.pb(v);
      v = fr[v];
    }
    reverse(res.begin(), res.end());
    return res;
  }
};

// BEGIN ALGO
enum CubeSide {
  FAR = 0, RIG = 1, NEA = 2, LEF = 3,
  UP = 4, DN = 5
};
// OX grows right
// OY grows up
// if cube rolls dx[d], dy[d], it turns(d)
const int dx[] = { 0, 1, 0, -1 };
const int dy[] = { 1, 0, -1, 0 };
int invturn[6] = { 2, 3, 0, 1, 5, 4 };
class Cube {
  int p[6]; // p[x] = label on side 'x'

  public:
  Cube() {
    for (int i = 0; i < 6; i++)
      p[i] = i;
  }
  Cube(int *ptr) {
    for (int i = 0; i < 6; i++)
      p[i] = ptr[i];
  }
  int& operator[](int x) {
    return p[x];
  }
  const int& operator[](int x) const {
    return p[x];
  }
  int getId() const { // nea * 4 + (rig - ...)
    int x = p[RIG];
    if (x >= max(p[NEA], p[FAR])) x--;
    if (x >= min(p[NEA], p[FAR])) x--;
    assert(0 <= x && x < 4);
    return p[FAR] * 4 + x;
  }
  /*
    dir in [0..3] => roll to direction 'dir'
      0 - roll forward; FAR --> DN
    dir = 4 => turn counterclockwise
    dir = 5 => turn clockwise
  */
  void turn(int dir) {
    const int fix[6][5] = {
      { 0, 5, 2, 4, 0 }, // a --> b, b --> c
      { 1, 5, 3, 4, 1 },
      { 0, 4, 2, 5, 0 }, // turn(0) x3
      { 1, 4, 3, 5, 1 }, // turn(1) x3
      { 0, 3, 2, 1, 0 }, // turn(0), 1, 2
      { 0, 1, 2, 3, 0 }  // turn(0), 3, 2
    };
    int np[6];
    memcpy(np, p, sizeof np);
    for (int i = 0; i < 4; i++)
      np[fix[dir][i + 1]] = p[fix[dir][i]];
    memcpy(p, np, sizeof np);
  }
  bool operator==(const Cube &c) const { return memcmp(p, c.p, sizeof p) == 0; } // NOT ALGO
  // c1.turn(a) * c2.turn(b) ==
  // == c.turn(a), c.turn(b)
  Cube operator*(const Cube &c2) {
    Cube res;
    for (int i = 0; i < 6; i++)
      res[i] = p[c2[i]];
    return res;
  }
};
int cubs[24][6]; // cube by its id
int cubgo[24][6]; // a.turn(d) ==> cubgo[a][d]
void cube_init() {
  memset(cubs, -1, sizeof cubs);
  Cube st;
  for (int a = 0; a < 4; a++, st.turn(0))
  for (int b = 0; b < 4; b++, st.turn(1))
  for (int c = 0; c < 4; c++, st.turn(4)) {
    int id = st.getId();
    assert(0 <= id && id < 24);
    if (cubs[id][0] < 0) {
      for (int i = 0; i < 6; i++)
        cubs[id][i] = st[i];
    }
    for (int i = 0; i < 6; i++)
      assert(cubs[id][i] == st[i]);
  }
  for (int i = 0; i < 24; i++)
    assert(cubs[i][0] >= 0);
  {
    Cube c;
    c.turn(0);
    assert(c[FAR] == UP);
  }

  // Precalc
  for (int a = 0; a < 24; a++)
  for (int d = 0; d < 6; d++) {
    Cube nst = cubs[a];
    nst.turn(d);
    cubgo[a][d] = nst.getId();
  }
}
void cube_check() {
  for (int d = 0; d < 6; d++) {
    Cube st;
    assert(st.getId() == 0);
    Cube c = st;
    for (int i = 0; i < 4; i++) {
      if (i) assert(!(c == st));
      c.turn(d);
    }
    assert(c == st);

    c.turn(d);
    assert(!(c == st));
    c.turn(invturn[d]);
    assert(c == st);
  }
}
// END ALGO

const int h = 8;
const int w = 8;

int main() {
  #ifdef DEBUG
  freopen("cube.in", "r", stdin);
  freopen("cube.out", "w", stdout);
  #endif

  cube_init();
  cube_check();

  char pos1[3], pos2[3];
  while (scanf("%s%s", pos1, pos2) == 2) {
    int cs[6];
    scanf("%d%d%d%d%d%d", &cs[NEA], &cs[FAR], &cs[UP], &cs[RIG], &cs[DN], &cs[LEF]);

//    assert(cubs[cubgo[0][1]

    int sx = pos1[0] - 'a', sy = pos1[1] - '1';
    int ex = pos2[0] - 'a', ey = pos2[1] - '1';

    Solver s(w * h * 24);
    for (int y = 0; y < h; y++)
    for (int x = 0; x < w; x++)
    for (int ty = 0; ty < 24; ty++) {
      int a = (y * w + x) * 24 + ty;
      for (int d = 0; d < 4; d++) {
        int nx = x + dx[d], ny = y + dy[d];
        if (nx < 0 || ny < 0 || nx >= w || ny >= h) continue;

        int nty = cubgo[ty][d];
        int b = (ny * w + nx) * 24 + nty;
        s.adde(a, b, cs[cubs[nty][DN]]);
      }
    }
    s.solve((sy * w + sx) * 24 + 0);
    int ans = 1e9;
    vi ap;
    for (int ety = 0; ety < 24; ety++) {
      int b = (ey * w + ex) * 24 + ety;
      int cans = s.get(b);
      if (cans < ans) {
        ans = cans;
        ap = s.getP(b);
      }
    }
    ans += cs[DN];
    printf("%d", ans);
    for (int i = 0; i < sz(ap); i++) {
      int pos = ap[i] / 24;
      int x = pos % w, y = pos / w;
      printf(" %c%c", 'a' + x, '1' + y);
    }
    printf("\n");
  }
  return 0;
}
