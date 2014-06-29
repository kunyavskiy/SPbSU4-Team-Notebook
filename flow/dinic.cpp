#include <stdio.h>
#include <ctype.h>
#include <iostream>
#include <math.h>
#include <string.h>
#include <algorithm>
#include <stdlib.h>
#include <time.h>
#include <queue>
#include <set>
#include <map>
#include <string>
#include <vector>
#include <utility>
#include <assert.h>

//#define M_PI 3.141592653589793238462643
#define eps 1e-8
#define inf ((int)1e9)
#define pb push_back
#define mp make_pair
#define taskname "flow2"

#define maxn 500
#define maxe 10000

using namespace std;

//BEGIN ALGO

//finds maxflow from "st" to "end"
//V - throughput of edges
//reverse edge has num "edge_num^1"

int Head[maxn], Next[2*maxe];
int Dest[2*maxe], e, mfl;
int V[2*maxe];
int u[maxn], uit, Lv[maxn], end, st, B[maxn];
long long flow;
int Q[maxn], bq, eq;

void add (int v1, int v2, int v)
{
  Dest[e]=v2, Next[e]=Head[v1];
  Head[v1]=e, V[e]=v, e++;
  assert(e<=2*maxe);
}

int dfs (int v, int fl)
{
  int nfl, nv, ne;
  u[v]=uit;
  if (v==end)
  {
    flow+=fl;
    return fl;
  }
  for ( ; B[v]!=-1; B[v]=Next[B[v]])
  {
    ne=B[v], nv=Dest[ne]; /*BOXNEXT*/
    if (u[nv]!=uit && V[ne]>=mfl && Lv[nv]==Lv[v]+1)
    {
      nfl=dfs(nv,min(V[ne],fl));
      if (nfl)
      {
        V[ne]-=nfl, V[ne^1]+=nfl;
        return nfl;
      }
    }
  }
  return 0;
}

long long dinic (int n, int _st, int _end)
{
  int v, r, i;
  flow=0;
  st=_st, end=_end, mfl=(1<<29);
  while (mfl)
  {
    while (true)
    {
      for (i=0; i<n; i++)
        Lv[i]=inf, B[i]=Head[i];
      Lv[st]=0, Q[0]=st, bq=0, eq=1;


      while (bq<eq)
      {
        v=Q[bq++];
        for (r=Head[v]; r!=-1; r=Next[r])
          if (V[r]>=mfl && Lv[Dest[r]]==inf) {
            Lv[Dest[r]]=Lv[v]+1;
            Q[eq++]=Dest[r];
          }
      }
      if (Lv[end]==inf)
        break;
      uit++;
      while (dfs(st,inf))
        uit++;
    }
    mfl>>=1;
  }
  return flow;
}

//END ALGO

int main()
{
  int i, v1, v2, v, n, m;
  freopen(taskname".in", "r", stdin);
  freopen(taskname".out", "w", stdout);
  memset(Head,-1,sizeof(Head));
  scanf("%d%d", &n, &m);
  assert(1<=n && n<=500);
  assert(1<=m && m<=10000);
  for (i=0; i<m; i++)
    scanf("%d%d%d", &v1, &v2, &v), v1--, v2--, add(v1,v2,v), add(v2,v1,0);
  printf("%I64d\n", dinic(n,0,n-1));
  return 0;
}

