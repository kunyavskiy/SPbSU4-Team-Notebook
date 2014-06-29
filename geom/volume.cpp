#include <stdio.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <assert.h>

#define mp make_pair

using namespace std;

typedef long long ll;
typedef pair <int, int> pii;

struct point {
  ll x, y, z;
  void read ()
  {
    int _x, _y, _z;
    scanf("%d%d%d", &_x, &_y, &_z);
    x=_x, y=_y, z=_z;
  }
};

bool operator < (const point& a, const point& b)
{
  if (a.x!=b.x)
    return a.x<b.x;
  if (a.y!=b.y)
    return a.y<b.y;
  return a.z<b.z;
}

bool operator == (const point& a, const point& b)
{
  return a.x==b.x && a.y==b.y && a.z==b.z;
}

//BEGIN ALGO

ll det (point a, point b, point c)
{
  return a.x*(b.y*c.z-b.z*c.y)-
         a.y*(b.x*c.z-b.z*c.x)-
         a.z*(b.y*c.x-b.x*c.y);
}

//MAXG - maximal number of faces AFTER triangulation
//complexity - O(n*logn)
const int MAXF = 11000;

point A[MAXF][3], B[3*MAXF]; /*BOXNEXT*/
int n, m, G[3*MAXF][1000], P[MAXF][3], e, ef, ee;
pair <pii, int> E[3*MAXF];
bool u[MAXF];

void dfs (int g)
{
  int j, v1, v2, i, ng;
  u[g]=1;
  for (i=0; i<3; i++)
  {
    v1=P[g][i];
    v2=P[g][(i+1)%3]; /*BOXNEXT*/
    ng=lower_bound(E,E+ee,mp(mp(min(v1,v2),max(v1,v2)),g))-E; /*BOXNEXT*/
    assert(E[ng]==mp(mp(min(v1,v2),max(v1,v2)),g));
    ng=E[ng^1].second;
    if (u[ng])
      continue;
    for (j=0; j<3; j++)
      if (P[ng][j]==v1 && P[ng][(j+1)%3]==v2)
      {
        swap(P[ng][j],P[ng][(j+1)%3]);
        break;
      }
    dfs(ng);
  }
}

int main()
{
  int i, l, r, j, cnt, x, y;
  ll ans=0;
  freopen ("volume.in", "r", stdin);
  freopen ("volume.out", "w", stdout);
  scanf("%d", &n);
  for (i=0; i<n; i++)
  {
    /*    //NOT ALGO
    scanf("%d", &cnt);
    */    //NOT ALGO
    cnt=3;  //NOT ALGO
    for (j=0; j<cnt; j++)
      B[e+j].read();
    for (j=2; j<cnt; j++) /*BOXNEXT*/
      A[ef][0]=B[e], A[ef][1]=B[e+j-1], A[ef][2]=B[e+j], ef++;
    e+=cnt;
  }
  sort(B,B+e);
  e=unique(B,B+e)-B;
  for (i=0; i<ef; i++)
  {
    for (j=0; j<3; j++)
    {
      l=0, r=e;
      while (r-l>1)
        (A[i][j]<B[(l+r)/2])?(r=(l+r)/2):(l=(l+r)/2);
      P[i][j]=l;
    }
    for (j=0; j<3; j++)
    {
      x=P[i][j], y=P[i][(j+1)%3];
      E[ee++]=mp(mp(min(x,y),max(x,y)),i);
    }
  }
  sort(E,E+ee), dfs(0);
  for (i=0; i<ef; i++)
    ans+=det(B[P[i][0]],B[P[i][1]],B[P[i][2]]);
  printf("%.8lf", ((double)(abs(ans)))/6.0);
  return 0;
}
//END ALGO
