// BEGIN ALGO
//   a[i][j] -- initial data
//   a[mt[i]][i] is result set
//   sum is minimized
//   TODO works on rectangular?

int a[maxn][maxn], fl[maxn],fr[maxn];
int p[maxn],d[maxn], mt[maxn];
bool used[maxn];
int n;

void iterate(int nv){
  memset(d,100,sizeof(d));
  memset(used,0,sizeof(used));
  
  d[n] = 0;
  int l = n;
  mt[n] = nv;
  
  while (mt[l] != -1){
    int k = mt[l];
    int nl = -1;
    used[l] = true;
    for (int i = 0; i < n; i++)
      if (!used[i]){ /*BOXNEXT*/
        if (d[i] > d[l] + a[k][i] + fl[k] + fr[i]) { /*BOXNEXT*/
          d[i] = d[l] + a[k][i] + fl[k] + fr[i];
          p[i] = l;
        }
        if (nl == -1 || d[nl] > d[i])
          nl = i;
      }
    
    int add = d[nl];
    
    for (int i = 0; i <= n; i++)
      if (used[i])
        fr[i] += add, fl[mt[i]] -= add;
      else
        d[i] -= add;
    l = nl;
  }
  
  while (l != n){
    mt[l] = mt[p[l]];
    l = p[l];
  }
}

void hungrian(){
  memset(mt,-1,sizeof(mt));
  for (int i = 0; i < n; i++)
    iterate(i);
  int ans = 0;
  for (int i = 0; i < n; i++)
    ans += a[mt[i]][i];
}
// END ALGO