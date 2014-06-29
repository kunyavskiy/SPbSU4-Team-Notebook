#include <cstdio>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <vector>

using namespace std;

#define pb push_back
#define sz(x) ((int)(x).size())
#define eprintf(...) fprintf(stderr, __VA_ARGS__)

typedef long long ll;

int sgn(ll x) { return x < 0 ? -1 : x > 0; }
struct pt {
  int x, y;
  pt() : x(0), y(0) {}
  pt(int x, int y) : x(x), y(y) {}
  pt operator+(const pt &p2) const { return pt(x + p2.x, y + p2.y); }
  pt operator-(const pt &p2) const { return pt(x - p2.x, y - p2.y); }
  int operator*(const pt &p2) const { return sgn(ll(x) * p2.y - ll(y) * p2.x); }
  int operator^(const pt &p2) const { return sgn(ll(x) * p2.x + ll(y) * p2.y); }
  bool operator==(const pt &p2) const { return x == p2.x && y == p2.y; }
};

// BEGIN ALGO
// ALL coordinates should be |x| <= 1e9, 2e9 is bad!

// x in [min(l, r), max(l, r)]
bool is_in(int l, int r, int x) {
  if (l > r) swap(l, r);
  return l <= x && x <= r;
}

// [a, b] contains p
bool contains(const pt &a, const pt &b,
              const pt &p) {
  if ((p - a) * (p - b) != 0) return false;
  return is_in(a.x, b.x, p.x) &&
         is_in(a.y, b.y, p.y);
}

// gon should not have self-intersections
// (self-touches are ok, even S=0 is ok)
// 0 - p is strictly outside
// 1 - p lies on border
// 2 - p is strictly inside
int inside(const vector<pt> &gon, const pt &p) {
  int cnt = 0;
  for (int i = 0; i < sz(gon); i++) {
    pt a = gon[i];
    pt b = gon[i + 1 >= sz(gon) ? 0 : i + 1];
    if (contains(a, b, p)) return 1;
    
    if (a.y > b.y) swap(a, b);
    if (a.y == b.y || p.y <= a.y || b.y < p.y)
       continue;

    if (ll(b.x - a.x) * ll(p.y - a.y) -
        ll(b.y - a.y) * ll(p.x - a.x) >= 0)
      cnt++;
  }
  return (cnt & 1) ? 2 : 0;
}

bool operator<(const pt &a, const pt &b) {
  if (a.y != b.y) return a.y < b.y;
  return a.x < b.x;
}

// gon should not have self-intersections
// self-touches are ok, consecutive points /*BOXNEXT*/
// on the same line are ok (not unique points are ok) /*BOXNEXT*/
// works ok with degenerate gons (sz(gon) >= 1)
// returns gon without consecutive points on the same line
void normalize(vector<pt> &gon) {
  gon.erase(unique(gon.begin(), gon.end()),
            gon.end());
  while (sz(gon) > 1 && gon[0] == gon.back())
    gon.pop_back();
  rotate(gon.begin(),
         min_element(gon.begin(), gon.end()),
         gon.end());
  
  int ptr = 1;
  for (int i = 1; i < sz(gon); i++) {
    int pr = ptr - 1;
    int ne = (i + 1 == sz(gon) ? 0 : i + 1);

    pt a = gon[pr], b = gon[i], c = gon[ne];
    if (((b - a) * (c - b)) != 0 ||
        ((b - a) ^ (c - b)) < 0) {
      gon[ptr++] = gon[i];
      continue;
    }
  }
  gon.resize(ptr);
}
/*BOXNEXT*/
// gon should not have self-intersections and should be
// convex (counter-clockwise order, i.e. oriented area >0)
// self-touches are not ok /*BOXNEXT*/
// consecutive points on the same line are not ok,
// {0,1,2} = {"OUTSIDE", "BORDER", "INSIDE"} /*BOXNEXT*/
int inside_convex(const vector<pt> &gon, const pt &p) {
  assert(!gon.empty());
  
  pt root = gon[0];
  if (root == p) return 1;
  if (sz(gon) <= 1) return 0;
  
  pt a = gon[1], b = gon.back();
  assert((a - root) * (b - root) > 0);
  if ((p - root) * (a - root) > 0) return 0;
  if ((p - root) * (b - root) < 0) return 0;
  
  if ((p - root) * (b - root) == 0)
    return contains(root, b, p) ? 1 : 0;
  
  // (p - root) * (a - root) <= 0
  // (p - root) * (b - root) > 0
  int L = 0, R = sz(gon) - 1;
  while (L + 1 < R) {
    int M = (L + R) / 2; /*BOXNEXT*/
    if ((p - root) * (gon[M] - root) > 0) R = M;
    else L = M;
  }
  pt x = gon[L], y = gon[R];
  
  if (contains(x, y, p)) return 1;
  if ((x - root) * (a - root) == 0 &&
      contains(root, x, p))
    return 1;
  return (p - x) * (y - x) > 0 ? 0 : 2;
}
// END ALGO

const char* res[] = {
  "OUTSIDE", "BORDER", "INSIDE"
};

int main() {
  #if 0
  freopen("norm.in", "r", stdin);
  freopen("norm.out", "w", stdout);
  int n;
  while (scanf("%d", &n) == 1) {
    vector<pt> pts(n);
    for (int i = 0; i < n; i++)
      scanf("%d%d", &pts[i].x, &pts[i].y);
    normalize(pts);
    printf("%d\n", sz(pts));
    for (int i = 0; i < sz(pts); i++)
      printf("%d %d\n", pts[i].x, pts[i].y);
  }
  #else
  freopen("inside.in", "r", stdin);
  freopen("inside.out", "w", stdout);

  int n;
  while (scanf("%d", &n) == 1) {
    vector<pt> pts(n);
    for (int i = 0; i < n; i++)
      scanf("%d%d", &pts[i].x, &pts[i].y);

    {
      ll sum = 0;
      for (int i = 0; i < n; i++) {
        pt a = pts[i], b = pts[(i + 1) % n];
        sum += ll(a.x) * b.y - ll(a.y) * b.x;
      }
      if (sum < 0)
        reverse(pts.begin(), pts.end());
    }

    int m;
    scanf("%d", &m);
    for (int i = 0; i < m; i++) {
      pt p;
      scanf("%d%d", &p.x, &p.y);
      printf("%s\n", res[inside_convex(pts, p)]);
    }
  }
  #endif
  return 0;
}
