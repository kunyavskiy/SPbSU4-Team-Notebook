// OpenCup XII Open Cup named after E. V. Pankratiev
// Stage 4: Grand Prix of Eastern Europe, Division 1, Sunday, November 4, 2012
// Problem G. General
#include <cstdio>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

#define mp make_pair
#define pb push_back
#define sz(x) ((int)(x).size())
#define eprintf(...) fprintf(stderr, __VA_ARGS__)
#ifdef _WIN32
#define LLD "%I64d"
#else
#define LLD "%lld"
#endif

typedef long long ll;
typedef pair<int, int> pii;

int sgn(ll x) { return x < 0 ? -1 : x > 0; }
struct pt {
  int x, y;
  pt() : x(0), y(0) {}
  pt(int x, int y) : x(x), y(y) {}
  pt operator+(const pt &p2) const { return pt(x + p2.x, y + p2.y); }
  pt operator-(const pt &p2) const { return pt(x - p2.x, y - p2.y); }
  pt operator*(const int &b) const { return pt(x * b, y * b); }
  int operator*(const pt &p2) const { return sgn(ll(x) * p2.y - ll(y) * p2.x); }
  int operator^(const pt &p2) const { return sgn(ll(x) * p2.x + ll(y) * p2.y); }
  bool operator==(const pt &p2) const { return x == p2.x && y == p2.y; }
  ll dist2() const { return ll(x) * x + ll(y) * y; }
};

bool operator<(const pt &a, const pt &b) {
  if (a.y != b.y) return a.y < b.y;
  return a.x < b.x;
}

bool good(const pt &a, const pt &b, const pt &c) {
  return (b - a) * (c - b) > 0;
}
pt _root;
bool rcmp(pt a, pt b) {
  a = a - _root; b = b - _root;
  int s = a * b;
  if (s) return s > 0;
  return a.dist2() < b.dist2();
}
void convex(vector<pt> &pts) {
  swap(pts[0], *min_element(pts.begin(), pts.end()));
  _root = pts[0];
  sort(pts.begin() + 1, pts.end(), rcmp);
  int ste = 0;
  for (int i = 0; i < sz(pts); i++) {
    while (ste >= 2 && !good(pts[ste - 2], pts[ste - 1], pts[i]))
      ste--;
    pts[ste++] = pts[i];
  }
  pts.resize(ste);
}

// BEGIN ALGO
// 'root/3' is strictly inside of convex hull
pt root;
// compares 'a' and 'b/3'
bool acmp3(pt a, pt b) {
  a = a * 3 - root; b = b - root; /*BOXNEXT*/
  bool x1 = pt() < a; // is in upper half-plane
  bool x2 = pt() < b;
  if (x1 != x2) return x1 > x2;
  return a * b > 0;
}

void rotate_for_tangent(vector<pt> &pts) {
  assert(sz(pts) >= 3);
  root = pts[0] + pts[1] + pts[2];
  int bi = 0;
  for (int i = 0; i < sz(pts); i++)
    if (acmp3(pts[i], pts[bi] * 3))
      bi = i; /*BOXNEXT*/
  rotate(pts.begin(), pts.begin() + bi, pts.end());
}

/*BOXNEXT*/
// No three points in hull are collinear, area is >0
// Counter-clockwise, first point has the least polar angle
// Let f(a)=vp(a - p, b - a), where (a,b) is a side of hull
// result.first = first x | f(x) < 0
// result.first = first x | f(x) >= 0 /*BOXNEXT*/
// literally, two nearest points on tangents, in order
// from left to right.
// Returns (-1, -1) if point lies inside or on a border
// !!!WARNING!!! Call rotate_for_tangent first, it will set 'root'
pii get_tangent_to_conv(const vector<pt> &pts,
                        const pt &p) { /*BOXNEXT*/
  #define modn(x) (((x) >= n) ? ((x) - n) : (x))
  #define f(x) ((pts[x] - p) * (pts[modn(x + 1)] - pts[x]))
  const int n = sz(pts);
  assert(n >= 3);
  
  if (p * 3 == root) return mp(-1, -1);
  
  int a = lower_bound(pts.begin(), pts.end(),
            p * 3, acmp3) - pts.begin();
  a = modn(a + n - 1);
  assert(0 <= a && a < n);
  // check for 'lies inside'
  if (f(a) >= 0) return mp(-1, -1);
  assert(f(a) < 0);
  
  int b = lower_bound(pts.begin(), pts.end(), /*BOXNEXT*/
                      root * 2 - p * 3, acmp3) -
          pts.begin();
  b = modn(b + n - 1);
  assert(0 <= b && b < n);
  assert(f(b) > 0);
  
  // now f(a)==-1, f(b)==1
  int L, R;
  L = a, R = b;
  while (modn(L + 1) != R) {
    int M = L;
    M += modn(R - L + n) / 2;
    M = modn(M);
    if (f(M) < 0) L = M;
    else R = M;
  }
  int res2 = R;
  
  L = b, R = a;
  while (modn(L + 1) != R) {
    int M = L;
    M += modn(R - L + n) / 2;
    M = modn(M);
    if (f(M) < 0) R = M;
    else L = M;
  }
  int res1 = R;
  
  assert(f(res1) < 0);
  assert(f(modn(res1 + n - 1)) >= 0);
  assert(f(res2) >= 0);
  assert(f(modn(res2 + n - 1)) < 0);
  return mp(res1, res2);
  #undef modn
  #undef f
}
// END ALGO

ll vp(const pt &a, const pt &b) {
  return ll(a.x) * b.y - ll(a.y) * b.x;
}

int main() {
  #ifdef DEBUG
  freopen("general.in", "r", stdin);
  freopen("general.out", "w", stdout);
  #endif

  int n, m;
  while (scanf("%d%d", &n, &m) == 2) {
    vector<pt> pts(n);
    for (int i = 0; i < n; i++)
      scanf("%d%d", &pts[i].x, &pts[i].y);
    convex(pts);
    rotate_for_tangent(pts);
    n = sz(pts);

    vector<ll> sums(n + 1);
    for (int i = 0; i < n; i++) {
      pt a = pts[i], b = pts[(i + 1) % n];
      sums[i + 1] = sums[i];
      sums[i + 1] += vp(a, b);
    }

    while (m --> 0) {
      pt p;
      scanf("%d%d", &p.x, &p.y);
      ll ans = sums[n];

      pii tg = get_tangent_to_conv(pts, p);
      if (tg.first >= 0) {
        assert(tg.first != tg.second);
        ans -= sums[tg.second];
        ans += sums[tg.first];
        if (tg.first >= tg.second)
          ans -= sums[n];
        ans += vp(pts[tg.first], p);
        ans += vp(p, pts[tg.second]);
      }
      ans = abs(ans);
      printf(LLD"."LLD"\n", ans / 2, ans % 2 * 5);
    }
  }
  return 0;
}
