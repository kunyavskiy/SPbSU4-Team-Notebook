#include <cstdio>
#include <cassert>
#include <cmath>
#include <vector>

using namespace std;

#define pb push_back
#define sz(x) ((int)(x).size())
#define eprintf(...) fprintf(stderr, __VA_ARGS__)

const double eps = 1e-9;

int sgn(double x) { return x < -eps ? -1 : x > eps; }
struct pt {
  double x, y;
  pt() : x(0), y(0) {}
  pt(double x, double y) : x(x), y(y) {}
  pt operator+(const pt &p2) const { return pt(x + p2.x, y + p2.y); }
  pt operator-(const pt &p2) const { return pt(x - p2.x, y - p2.y); }
  pt operator*(const double b) const { return pt(x * b, y * b); }
  int operator*(const pt &p2) const { return sgn(x * p2.y - y * p2.x); }
  int operator^(const pt &p2) const { return sgn(x * p2.x + y * p2.y); }
  double dist2() const { return x * x + y * y; }
  bool operator<(const pt &p2) const { return fabs(x - p2.x) > eps ? x < p2.x : y < p2.y - eps; }
  bool operator==(const pt &p2) const { return fabs(x - p2.x) < eps && fabs(y - p2.y) < eps; }
};
struct line {
  double a, b, c;
  line() : a(0), b(0), c(1) {}
  line(const pt &p1, const pt &p2) {
    a = p1.y - p2.y;
    b = p2.x - p1.x;
    c = -a * p1.x - b * p1.y;
  }
  line(double a, double b, double c) : a(a), b(b), c(c) {}
  double distz(const pt &p) const { return a * p.x + b * p.y + c; }
  double norm2() const { return a * a + b * b; }
  double dists(const pt &p) const { return distz(p) / sqrt(norm2()); }
  double dist(const pt &p) const { return fabs(distz(p)) / sqrt(norm2()); }
  int side(const pt &p2) const { return sgn(a * p2.x + b * p2.y + c); }
};

#define sqr(x) ((x) * (x))

// BEGIN ALGO
// Assumption: line.norm2() > 0
// Assumption: r >= 0

// returns 0-2 _unique_ points, which lie both on circle and line
// res[0] lies on the left side of (l.a, l.b)
// for example (center=0,0, r=2, 0x+1y-1=0): (-sqrt(3), 1), (sqrt(3),1)

pt operator&(const line &l1, const line &l2) {
  double d = l1.a * l2.b - l1.b * l2.a;
  assert(fabs(d) > eps);
  pt res(
    (l1.b * l2.c - l1.c * l2.b) / d,
    (l1.a * l2.c - l1.c * l2.a) / -d
  );
  assert(l1.side(res) == 0);
  assert(l2.side(res) == 0);
  return res;
}

bool is_in(double l, double r, double x) {
  if (l > r) swap(l, r);
  return l <= x + eps && x - eps <= r;
}
/*BOXNEXT*/
bool contains(const pt &a, const pt &b, const pt &p) {
  if ((p - a) * (p - b) != 0) return false;
  return is_in(a.x, b.x, p.x) &&
         is_in(a.y, b.y, p.y);
}

vector<pt> cross(const pt &center, double r,
                 const line &l) {
  double di = l.distz(center);
  double d2 = l.norm2();
  assert(fabs(d2) > eps);
  pt mid = center + pt(l.a, l.b) * (-di / d2);
  #ifdef DEBUG
  assert(l.side(mid) == 0);
  #endif

  double s = r * r - di * di / d2;
  if (s < -eps) return vector<pt>();
  if (fabs(di * di - r * r * d2) < eps)
    return vector<pt>(1, mid);

  pt off = pt(-l.b, l.a) * sqrt(s / d2);
  assert(fabs(off.dist2() - s) < eps);
  vector<pt> res;
  res.pb(mid + off);
  res.pb(mid - off);
  return res;
}




// returns 0-2 _unique_ points.
// res[0] lies on the left side of line (a-->b)
vector<pt> cross(const pt &a, double r1, pt b,
                 double r2) {
  b = b - a;
  if (fabs(b.dist2()) < eps) {
    if (fabs(r1) < eps && fabs(r2) < eps)
      return vector<pt>(1, a);

    assert(fabs(r1 - r2) > eps);
    return vector<pt>();
  }
  vector<pt> res = cross(b, r2,
      line(2 * b.x, 2 * b.y,
           sqr(r2) - sqr(r1) -
           sqr(b.x) - sqr(b.y)));
  for (int i = 0; i < sz(res); i++)
    res[i] = res[i] + a;
  return res;
}

// returns tangent line to (cent,r) in point (p)
// p should lie on the circle. !! r > 0 !!
line get_tangent(const pt &cent, /*BOXNEXT*/
                 const double &r, const pt &p){ /*BOXNEXT*/
  assert(fabs((p - cent).dist2() - r * r) < eps);
  assert(fabs(r) > eps);
  pt vec = p - cent;
  return line(
    vec.x, vec.y, -vec.x * p.x - vec.y * p.y
  );
}

// returns 0-2 _unique_ points such that (p-->x) is
// tangent to the circle. When r=0, returns no more than 1 point
// res[0] lies on the left side of line (p-->cent)
vector<pt> tangent(const pt &p,
                   const pt &cent,
                   const double &r) {
  vector<pt> res = cross(
    pt(), r,
    line(cent.x - p.x, cent.y - p.y, r * r)
  );
  for (int i = 0; i < sz(res); i++)
    res[i] = res[i] + cent;
  return res;
}







// returns 0, 2 or 4 _not unique_ lines
// first two are outer (to the left of a-->b, to the right),
// last two are inner (to left
vector<line> tangent(pt a, double r1,
                     pt b, double r2) {
  if (fabs((b - a).dist2()) < eps) {
    assert(fabs(r1 - r2) > eps);
    return vector<line>();
  }
  vector<line> res;
  for (int t = 1; t >= -1; t -= 2)
  for (int s1 = -1; s1 <= 1; s1 += 2) {
    int s2 = t / s1;
    vector<pt> base = cross(
      pt(), 1, /*BOXNEXT*/
      line(b.x - a.x, b.y - a.y, s1 * r1 - s2 * r2)
    );
    if (!base.empty()) {
      line l;
      l.a = base[0].x;
      l.b = base[0].y;
      l.c = -l.a * a.x - l.b * a.y + s1 * r1;
      res.pb(l);
    }
  }
  return res;
}

//cross [a1,b1] with [a2,b2]
vector<pt> cross(pt a1, pt b1, pt a2, pt b2) {
  if ((b1 - a1) * (b2 - a2) == 0) {
    if (b1 < a1) swap(a1, b1);
    if (b2 < a2) swap(a2, b2);
    if (a2 < a1) swap(a1, a2), swap(b1, b2);
    if (b1 < a2) return vector<pt>();
    if (b1 == a2) return vector<pt>(1, b1);
    vector<pt> res;
    res.pb(a2);
    res.pb(b1);
    return res;
  }
  pt res = line(a1, b1) & line(a2, b2);
  if (!contains(a1, b1, res) ||
      !contains(a2, b2, res))
        return vector<pt>();
  return vector<pt>(1, res);
}
/*BOXNEXT*/
double dist(const pt &a, const pt &b, const pt &p) {
  if (a == b) return sqrt((a - p).dist2());
  line l(a, b);
  pt mid = p + pt(l.a, l.b) *
               (-l.distz(p) / l.norm2());
  if (contains(a, b, mid))
    return sqrt((mid - p).dist2());
  return sqrt(min((a - p).dist2(),
                  (b - p).dist2()));
}

//cross [a1,b1] && [a2,b2]
double dist(const pt &a1, const pt &b1, const pt &a2, const pt &b2) {
  double ans = dist(a1, b1, a2);
  ans = min(ans, dist(a1, b1, b2));
  ans = min(ans, dist(a2, b2, a1));
  ans = min(ans, dist(a2, b2, b1));
  return ans;
}
// END ALGO

void check_cross(const pt &cent, const double &r, const line &l, int need_cnt) {
  vector<pt> res = cross(cent, r, l);
  printf("check circle&line\n");
  for (int i = 0; i < sz(res); i++) {
    printf("  %.2lf %.2lf\n", res[i].x, res[i].y);
    assert(l.side(res[i]) == 0);
    assert(fabs((cent - res[i]).dist2() - r * r) < eps);
  }
  assert(sz(res) == need_cnt);
}
void check_cross(const pt &c1, const double &r1, const pt &c2, const double &r2, int need_cnt) {
  vector<pt> res = cross(c1, r1, c2, r2);
  printf("check two circles\n");
  for (int i = 0; i < sz(res); i++) {
    printf("  %.2lf %.2lf\n", res[i].x, res[i].y);
    assert(fabs((c1 - res[i]).dist2() - r1 * r1) < eps);
    assert(fabs((c2 - res[i]).dist2() - r2 * r2) < eps);
  }
  assert(sz(res) == need_cnt);
}
void check_tangent(const pt &p, const pt &c, const double &r, int need_cnt) {
  vector<pt> res = tangent(p, c, r);
  printf("check tangent of (pt,circle)\n");
  for (int i = 0; i < sz(res); i++) {
    printf("  %.2lf %.2lf\n", res[i].x, res[i].y);
    assert(fabs((c - res[i]).dist2() - r * r) < eps);
    assert(((res[i] - p) ^ (res[i] - c)) == 0);
  }
  assert(sz(res) == need_cnt);
}
void check_tangent(const pt &a, const double &r1, const pt &b, const double &r2, int need_cnt) {
  vector<line> res = tangent(a, r1, b, r2);
  printf("check tangent of (circle,cirlce)\n");
  for (int i = 0; i < sz(res); i++) {
    line l = res[i];
    printf("  %.2lf %.2lf %.2lf\n", l.a, l.b, l.c);
    assert(fabs(l.dist(a) - r1) < eps);
    assert(fabs(l.dist(b) - r2) < eps);
  }
  assert(sz(res) == need_cnt);
}
void check_cross(const pt &a1, const pt &b1, const pt &a2, const pt &b2, int need_cnt) {
  eprintf("check cross of segments\n");
  vector<pt> res = cross(a1, b1, a2, b2);
  for (int i = 0; i < sz(res); i++) {
    eprintf("  %.2lf %.2lf\n", res[i].x, res[i].y);
    assert(contains(a1, b1, res[i]));
    assert(contains(a2, b2, res[i]));
  }
  assert(sz(res) == need_cnt);
}

int main() {
  check_cross(pt(4, 10), 53, line(10, 6, -13), 2);
  check_cross(pt(4, 10), 10, line(1, 0, -14),  1);
  check_cross(pt(4, 10), 9.9, line(0, 1, 0),   0);

  check_cross(pt(4, 10), 5, pt(0, 0), 6,  2);
  check_cross(pt(0, 0), 1, pt(3, 0), 2,   1);
  check_cross(pt(0, 0), 1, pt(4, 0), 2,   0);

  check_tangent(pt(1, 2), pt(4, 0), 2,   2);
  check_tangent(pt(0, 4), pt(3, 0), 5,   1);
  check_tangent(pt(1, 2), pt(4, 0), 100, 0);

  check_tangent(pt(3, 0), 1.5, pt(0, 4), 2.5, 4);
  check_tangent(pt(3, 0), 2.5, pt(0, 4), 2.5, 4);
  check_tangent(pt(3, 0), 1, pt(0, 4), 5, 2);
  check_tangent(pt(3, 0), 1, pt(0, 4), 6, 2);
  check_tangent(pt(3, 0), 0.5, pt(0, 4), 6, 0);

  check_cross(pt(2, 0), pt(1, 0), pt(1.5, -1), pt(1.5, 0), 1);
  check_cross(pt(2, 0), pt(1, 0), pt(1.5, 0), pt(2.5, 0), 2);

  return 0;
}
