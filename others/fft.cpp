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
#include <sstream>
#include <iostream>

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
#define strstr strstr_wregthrtu

using namespace std;

typedef long double ld;
typedef long long ll;
typedef vector<ll> vll;
typedef vector<int> vi;
typedef vector<vi> vvi;
typedef vector<bool> vb;
typedef vector<vb> vvb;
typedef pair<int, int> pii;
typedef pair<ll,ll> pll;
typedef pair<ld,ld> pld;

const int inf = 1e9;
const double eps = 1e-9;
const double INF = inf;
const double EPS = eps;

//BEGIN ALGO

struct cd{ld x,y;};

const int n = 19;
const int m = (1<<n);
const int maxm = m+20;

cd a[2][maxm], c[maxm]; 
int h, q, rev[maxm];
int x[2][maxm], d[maxm], z[maxm];
ld si[maxm*3/4+1];

void mktab () { //precalc
  ld ta, tb;
  int i;
  rev[0] = 0;
  for (q=1; q<m; q<<=1)
    for (i=0; i<q; i++)
      rev[i] <<= 1, rev[i+q]=rev[i]+1;
  ta=2.0*M_PI/m;
  q=m>>2;
  h=m>>1;
  for (i=0; i<=(m>>3); i++) /*BOXNEXT*/
    tb=ta*i, si[q-i]=cos(tb), si[i]=sin(tb);
  for (i=1; i<=q; i++)
    si[q+i] = si[q-i];
  for (i=0; i<=q; i++)
    si[h+i]=-si[i];
}

void step (int k, cd *from, cd *to, ld kk) {
  int st=m>>k, cc=st, dd=0, i;
  cd t, u, v, r;
  r.x=1.0;  r.y=0.0;
  for (i=0; i<h; i++) {
    t=from[2*i];
    u=from[2*i+1];
    v.x=u.x*r.x-u.y*r.y;
    v.y=u.x*r.y+u.y*r.x;
    to[i].x=t.x+v.x;
    to[i].y=t.y+v.y;
    to[h+i].x=t.x-v.x;
    to[h+i].y=t.y-v.y;
    if (!--cc) /*BOXNEXT*/
      dd+=st, cc=st, r.x=si[dd+q], r.y=si[dd]*kk;
  }
}

void fourier (cd *from, cd *to, ld kk) {
  int i;
  for (i=0; i<n; i++)
    if (i&1) step(i+1,to,from,kk);
    else step(i+1,from,to,kk);
  if (!(n&1))
    memcpy(to,from,m*sizeof(to[0]));
}


/*BOXNEXT*/
void tuda2in1 (int *from1, int *from2, cd *to, cd *tmp) {
  for (int i=0; i<m; i++)  {
    tmp[rev[i]].x=from1[i];
    tmp[rev[i]].y=from2[i];
  }
  fourier(tmp,to,1.0);
}

void tuda (int *from, cd *to, cd *tmp) {
  memset(tmp,0,m*sizeof(tmp[0]));
  for (int i=0; i<m; i++)
    tmp[rev[i]].x=from[i];
  fourier(tmp,to,1.0);
}

void obratno (cd *src, cd *tmp, int *ans) {
  int i;
  for (i=0; i<m; i++)
    tmp[rev[i]] = src[i];
  fourier(tmp,src,-1.0);
  memset(ans,0,m*sizeof(ans[0]));
  for (i=0; i<m; i++)
    ans[i]=(int)(src[i].x/m+0.5);
}

void multiply(){
  int i;
  mktab();
  
  //simple fft
  
  /* // NOT ALGO
  for (it=0; it<2; it++)
    tuda(x[it],a[it],c);
  for (i=0; i<m; i++)
  {
    ld tmp=a[0][i].x*a[1][i].x-a[0][i].y*a[1][i].y;
    a[0][i].y=a[0][i].x*a[1][i].y+a[0][i].y*a[1][i].x;
    a[0][i].x=tmp;
  }
  obratno(a[0],c,z)
  */ // NOT ALGO
  
  //2in1 fft
  
  tuda2in1(x[0],x[1],a[0],c);
  cd fs, sc;
  for (i=0; i<m; i++)
  {
    int p=(m-i)&(m-1);
    fs.x=(a[0][i].x+a[0][p].x)/2.0;
    fs.y=(a[0][i].y-a[0][p].y)/2.0;
    sc.x=(a[0][i].y+a[0][p].y)/2.0;
    sc.y=(-a[0][i].x+a[0][p].x)/2.0;
    a[1][i].x=fs.x*sc.x-fs.y*sc.y;
    a[1][i].y=fs.x*sc.y+fs.y*sc.x;
  }
  obratno(a[1],c,z);
}
//END ALGO


int zn=1;
char s[maxm];

int main() {
  int it, len, i;
  for (it=0; it<2; it++)
  {
    scanf("%s", s);
    int slen=strlen(s);
    if (s[0]=='-')
    {
      slen--, zn*=-1;
      for (int j=0; j<slen; j++)
        s[j]=s[j+1];
      s[slen]=0;
    }
    for (i=0; i<slen; i++)
      x[it][slen-1-i]=(int)(s[i]-'0');
  }  

  TIMESTAMP(beginfourie);
  multiply();
  TIMESTAMP(endfourie);

  bool fl=0;        
  len=0;
  for (i=0; i<m-2; i++)
  {
    z[i+1]+=z[i]/10, z[i]%=10;
    assert(0<=z[i] && z[i]<10);
    if (z[i])
      fl=1;
  }
  if (fl && zn==-1)
    printf("-");
  memset(s,0,sizeof(s));
  for (i=m-1; i>=0; i--)
    if (i==0 || len || z[i]>0)
      s[len++]=z[i]+'0';
  printf("%s\n", s);
  TIMESTAMP(end);
  return 0;
}
