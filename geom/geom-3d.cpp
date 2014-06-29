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

using namespace std;

typedef double ld;
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

struct point{
  ld x,y,z;
  point(ld x,ld y, ld z):x(x),y(y),z(z){}
  point(){ x = y = z = 0;}
  bool operator==(const point& a) const{
    return fabs(x - a.x) < eps && fabs(y - a.y) < eps && fabs(z - a.z) < eps;
  }
  bool load(){
    double _x,_y,_z;
    if (scanf("%lf%lf%lf",&_x,&_y,&_z) != 3) return 0;
    x = _x, y = _y, z = _z;
    return 1;
  }
  void print(){
    printf("%lf %lf %lf\n",(double)x, (double)y, (double)z);
  }
  void eprint(){
    eprintf("%lf %lf %lf\n",(double)x, (double)y, (double)z);
  }
  ld dist2() const {
    return x*x+y*y+z*z;
  }
  ld dist() const{
    return sqrt(dist2());
  }
};

inline point operator+(const point& a,const point& b){
  return point(a.x+b.x,a.y+b.y,a.z+b.z);
}

inline point operator-(const point& a,const point& b){
  return point(a.x-b.x,a.y-b.y,a.z-b.z);
}

inline point operator*(const point& a,ld t){
  return point(a.x*t, a.y*t, a.z*t);
}

inline ld det(ld a,ld b,ld c,ld d){
  return a*d - b*c;
}

inline ld det(ld a11,ld a12,ld a13,ld a21,ld a22, ld a23, ld a31, ld a32,ld a33){
  return a11*det(a22,a23,a32,a33) - a12*det(a21,a23,a31,a33) + a13*det(a21,a22,a31,a32);
}

int sgn(ld x){
  return (x > eps) - (x < -eps);
}

inline ld det(const point& a,const point& b,const point& c){
  return det(a.x,a.y,a.z,b.x,b.y,b.z,c.x,c.y,c.z);
}

point vp(const point& a,const point& b){
  return point(det(a.y,a.z,b.y,b.z),-det(a.x,a.z,b.x,b.z),det(a.x,a.y,b.x,b.y));
}

ld sp(const point& a,const point& b){
  return a.x * b.x + a.y * b.y + a.z * b.z;
}

// BEGIN ALGO

// vp is vector product (point)
// sp is scalar product (ld)

struct line{
  point p,v;
  line(){}; /*BOXNEXT*/
  line(const point& p,const point& v):p(p),v(v){
    assert(!(v == point()));
  }
  bool on(const point& pt) const{
    return vp(pt - p, v) == point();
  }
};
struct plane {
  point n;
  ld d;
  plane() : d(0) {}
  plane(const point &p1, const point &p2,
        const point &p3) {
    n = vp(p2 - p1, p3 - p1);
    d = -sp(n, p1);
    assert(side(p1) == 0);
    assert(side(p2) == 0);
    assert(side(p3) == 0);
  }
  int side(const point &p) const {
    return sgn(sp(n, p) + d);
  }
};

int intersec(const line& l1, const line& l2,
             point& res){
  assert(!(l1.v == point()));
  assert(!(l2.v == point()));
  if (vp(l1.v,l2.v) == point()){
    if (vp(l1.v, l1.p - l2.p) == point())
      return 2; // same
    return 0; // parallel
  }
  point n = vp(l1.v,l2.v);
  point p = l2.p - l1.p;
  if (sgn(sp(n,p)))
    return 0; // skew
  ld t;
  if (sgn(n.x))
    t = (p.y * l2.v.z - p.z * l2.v.y) / n.x;
  else if (sgn(n.y))
    t = (p.z * l2.v.x - p.x * l2.v.z) / n.y;
  else if (sgn(n.z))
    t = (p.x * l2.v.y - p.y * l2.v.x) / n.z;
  else
    assert(false);
  res = l1.p + l1.v * t;
  assert(l1.on(res)); assert(l2.on(res));
  return 1; // intersects
}

ld dist(const line& l1,const line& l2){
  point ret = l1.p - l2.p; /*BOXNEXT*/
  ret = ret - l1.v * (sp(l1.v,ret) / l1.v.dist2());
  point tmp = l2.v - l1.v *
              (sp(l1.v,l2.v) / l1.v.dist2());
  if (sgn(tmp.dist2()))  /*BOXNEXT*/
    ret = ret - tmp * (sp(tmp,ret) / tmp.dist2());
  assert(fabs(sp(ret,l1.v)) < eps);
  assert(fabs(sp(ret,tmp)) < eps);
  assert(fabs(sp(ret,l2.v)) < eps);
  return ret.dist();
}

void closest(const line& l1,const line& l2,
             point& p1,point& p2){
  if (vp(l1.v,l2.v) == point()){
    p1 = l1.p;
    p2 = l2.p - l1.v * /*BOXNEXT*/
         (sp(l1.v,l2.p - l1.p) / l1.v.dist2());
    return;
  }
  point p = l2.p   - l1.p;
  ld t1 = (
           sp(l1.v,p) * l2.v.dist2() -
           sp(l1.v,l2.v) * sp(l2.v,p)
          ) / vp(l1.v,l2.v).dist2();
  ld t2 = (
             sp(l2.v,l1.v) * sp(l1.v,p) -
             sp(l2.v,p) * l1.v.dist2()
          ) / vp(l2.v,l1.v).dist2();
  p1 = l1.p + l1.v * t1;
  p2 = l2.p + l2.v * t2;
  assert(l1.on(p1));
  assert(l2.on(p2));
}

int cross(const line &l, const plane &pl,
          point &res) {
  ld d = sp(pl.n, l.v);
  if (sgn(d) == 0) {
    return (pl.side(l.p) == 0) ? 2 : 0;
  }
  ld t = (-sp(pl.n, l.p) - pl.d) / d;
  res = l.p + l.v * t;
  #ifdef DEBUG
  assert(pl.side(res) == 0);
  #endif
  return 1;
}

bool cross(const plane& p1,const plane& p2,
           const plane& p3, point& res){
  ld d = det(p1.n,p2.n,p3.n);
  if (sgn(d) == 0) {
     return false;
  }
  point px(p1.n.x, p2.n.x, p3.n.x);
  point py(p1.n.y, p2.n.y, p3.n.y);
  point pz(p1.n.z, p2.n.z, p3.n.z);
  point p(-p1.d,-p2.d,-p3.d);
  res = point(
              det(p,py,pz)/d,
              det(px,p,pz)/d,
              det(px,py,p)/d
             );
  #ifdef DEBUG
  assert(p1.side(res) == 0);
  assert(p2.side(res) == 0);
  assert(p3.side(res) == 0);
  #endif
  return true;
}

int cross(const plane &p1, const plane &p2,
          line &res) {
  res.v = vp(p1.n, p2.n);
  if (res.v == point()) {
    if ( (p1.n * (p1.d / p1.n.dist2())) ==
         (p2.n * (p2.d / p2.n.dist2())) )
      return 2;
    else
      return 0;
  }
  plane p3;
  p3.n = res.v;
  p3.d = 0;
  bool ret = cross(p1, p2, p3, res.p);
  assert(ret);
  assert(p1.side(res.p) == 0);
  assert(p2.side(res.p) == 0);
  return 1;
}
// END ALGO



int main(){
  freopen(TASKNAME".in","r",stdin);
  #ifdef LOCAL
  freopen(TASKNAME".out","w",stdout);
  #endif

  {
    line l;
    l.p = point(1, 1, 1);
    l.v = point(1, 0, -1);
    plane p(point(10, 11, 12), point(9, 8, 7), point(1, 3, 2));
    point res;
    assert(cross(l, p, res) == 1);
  }
  {
    plane p1(point(1, 2, 3), point(4, 5, 6), point(-1, 5, -4));
    plane p2(point(3, 2, 1), point(6, 5, 4), point(239, 17, -42));
    line l;
    assert(cross(p1, p2, l) == 1);
  }
  {
    plane p1(point(1, 2, 3), point(4, 5, 6), point(-1, 5, -4));
    plane p2(point(1, 2, 3), point(7, 8, 9), point(3, -1, 10));
    line l;
    assert(cross(p1, p2, l) == 2);
  }
  {
    plane p1(point(1, 2, 3), point(4, 5, 6), point(-1, 5, -4));
    plane p2(point(1, 2, 4), point(4, 5, 7), point(-1, 5, -3));
    line l;
    assert(cross(p1, p2, l) == 0);
  }

  line l1,l2;
  while (l1.p.load()){
    l1.v.load(); l1.v = l1.v - l1.p;
    l2.p.load();
    l2.v.load(); l2.v = l2.v - l2.p;
    if (l1.v == point() || l2.v == point()) continue;
    point res;
    int cnt = intersec(l1,l2,res);
    ld d = dist(l1,l2);
    if (fabs(d) < eps)
      assert(cnt >= 1);
    else
      assert(cnt == 0);
    point p1,p2;
    closest(l1,l2,p1,p2);
    assert(fabs((p1-p2).dist() - d) < eps);
  }
  plane a(point(1,0,0),point(0,1,0),point(0,0,1));
  plane b(point(-1,0,0),point(0,-1,0),point(0,0,-1));
  line l;
  assert((cross(a,b,l))==0);
  TIMESTAMP(end);
  return 0;
}

