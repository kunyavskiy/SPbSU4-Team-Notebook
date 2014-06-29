#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cmath>
#include <algorithm>
#include <vector>
#include <string>
#include <deque>
#include <map>
#include <set>

using namespace std;

#define pb push_back
#define mp make_pair
#define eprintf(...) fprintf(stderr, __VA_ARGS__)
#define sz(x) ((int)(x).size())
#define TASKNAME "forest"

typedef long long ll;
typedef vector<ll> vll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<bool> vb;
typedef vector<vb> vvb;
typedef pair<int, int> pii;
typedef long double ld;
//typedef double ld;

#define EPS 2e-8
int sgn(ld x) { return x < -EPS ? -1 : x > EPS; }

struct pt {
  ld x, y;
  pt(ld _x = 0, ld _y = 0) : x(_x), y(_y) {}
  bool operator<(const pt &p2) const {
    if (fabs(x - p2.x) > EPS) return x < p2.x;
    return y < p2.y - EPS;
  }
  pt operator-(const pt &p2) const { return pt(x - p2.x, y - p2.y); }
  int operator*(const pt &p2) const { return sgn(x * p2.y - y * p2.x); }
  ld dist2() const { return x * x + y * y; }
};

pt root;
bool rcmp(pt a, pt b) {
  a = a - root; b = b - root;
  int s = a * b;
  if (s) return s > 0;
  return a.dist2() < b.dist2() - EPS;
}
bool good(const pt &a, const pt &b, const pt &c) { return (b - a) * (c - a) > 0; }
vector<pt> convex(vector<pt> pts) {
  swap(*pts.begin(), *min_element(pts.begin(), pts.end()));
  root = pts[0];
  sort(pts.begin() + 1, pts.end(), rcmp);
  vector<pt> st(sz(pts)); int ste = 0;
  for (int i = 0; i < sz(pts); i++) {
    while (ste >= 2 && !good(st[ste - 2], st[ste - 1], pts[i]))
      ste--;
    st[ste++] = pts[i];
  }
  st.resize(ste);
  return st;
}

struct line {
  ld a, b, c;
  ld cang;

  line() {}
  line(const pt &p1, const pt &p2) {
    a = p1.y - p2.y;
    b = p2.x - p1.x;
    c = -a * p1.x - b * p1.y;

    ld d = sqrt(a * a + b * b);
    a /= d;
    b /= d;
    c /= d;
    cang = atan2(b, a);
  }
  int side(const pt &p) const {
    return sgn(a * p.x + b * p.y + c);
  }

  bool operator||(const line &l2) const { return sgn(a * l2.b - b * l2.a) == 0; }

  pt operator&(const line &l2) const {
    ld d = a * l2.b - b * l2.a;
    assert(fabs(d) > EPS);
    pt res(
      (b * l2.c - c * l2.b) / d,
      (a * l2.c - c * l2.a) / -d
    );
    return res;
  }
  bool angEq(const line &l2) const {
    ld d = fabs(cang - l2.cang);
    if (2 * M_PI - d < d) d = 2 * M_PI - d;
    return d < EPS;
  }
  bool operator<(const line &l2) const {
    ld d = fabs(cang - l2.cang);
    if (d < EPS) return c < l2.c - EPS;
    return cang < l2.cang;
  }
};

// BEGIN ALGO
vector<pt> cross(vector<line> all) {
  { // \emph{PANIC. Better to rewrite!!!}
    const ld BBOXC = 5.5e10;
    vector<pt> bbox;
    bbox.pb(pt( BBOXC, -BBOXC));
    bbox.pb(pt( BBOXC,  BBOXC));
    bbox.pb(pt(-BBOXC,  BBOXC));
    bbox.pb(pt(-BBOXC, -BBOXC));
    bbox.pb(bbox[0]);
    for (int i = 0; i < 4; i++)
      all.pb(line(bbox[i], bbox[i + 1]));
  }

  sort(all.begin(), all.end());
  int off = 0;
  for (int i = 1; i < sz(all); i++)
    if (all[i - 1].angEq(all[i])) {
      off++;
    } else
      all[i - off] = all[i];
  all.resize(sz(all) - off);

  vector<pt> pts;
  vector<line> ss;
  ss.pb(all[0]);
  int deleted = 0;

  for (int i = 1; i < sz(all); i++) {
    int pcnt = sz(pts);
    while (sz(pts) > deleted && /*BOXNEXT*/
           all[i].side(pts[sz(pts) - 1]) <= 0) {
      pts.erase(pts.end() - 1);
      ss.erase(ss.end() - 1);
    }
    if (sz(pts) == deleted && pcnt) {  /*BOXNEXT*/
      if (pt(ss[sz(ss) - 1].a, ss[sz(ss) - 1].b) *
          pt(all[i].a, all[i].b) <= 0)
        return vector<pt>();
    } else {
      while (sz(pts) > deleted &&
             all[i].side(pts[deleted]) <= 0)
      	  ++deleted;
    } /*BOXNEXT*/
    if (ss[sz(ss) - 1] || all[i]) // parallel !!!
      return vector<pt>();  /*BOXNEXT*/
    pt cpt = ss[sz(ss) - 1] & all[i]; // intersect
    if (ss[deleted].side(cpt) >= 0) {
      pts.pb(cpt);
      ss.pb(all[i]);
    }
  }
  ss.erase(ss.begin(),ss.begin()+deleted);
  pts.erase(pts.begin(),pts.begin()+deleted);
  if (sz(ss) == 1) return pts;
  pts.pb(ss[0] & ss[sz(ss) - 1]);
  return pts;
}
// END ALGO

int main() {
  freopen(TASKNAME ".in", "r", stdin);
  freopen(TASKNAME ".out", "w", stdout);

  int n;
  while (scanf("%d", &n) >= 1) {
    vector<pt> spts(n);
    for (int i = 0; i < n; i++) {
      double x, y;
      scanf("%lf%lf", &x, &y);
      spts[i] = pt(x, y);
    }

    vector<pt> cpts = convex(spts);

    vector<line> all;
    for (int i = 0; i + 1 < n; i++)
      all.pb(line(spts[i], spts[i + 1]));

    vector<pt> pts = cross(all);

    if (!sz(pts)) printf("Impossible\n");
    else {
      for (int i = 0; i < sz(cpts); i++) {
        int ne = (i + 1 == sz(cpts)) ? 0 : i + 1;
        all.pb(line(cpts[i], cpts[ne]));
      }
      vector<pt> pts2 = cross(all);
      sort(pts.begin(), pts.end());
      sort(pts2.begin(), pts2.end());

      vector<pt> res(sz(pts));
      res.resize(set_difference(pts.begin(), pts.end(), pts2.begin(), pts2.end(), res.begin()) - res.begin());
      if (sz(res)) {
        printf("Possible\n%.18lf %.18lf\n", double(res[0].x), double(res[0].y));
      } else
        printf("Impossible\n");
    }
  }
  return 0;
}
	