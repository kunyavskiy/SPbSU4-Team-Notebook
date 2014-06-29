#include <cstdio>
#include <cstring>
#include <cassert>
#include <cctype>
#include <ctime>
#include <algorithm>
#include <set>

// http://acm.timus.ru/problem.aspx?space=1&num=1517

using namespace std;

#define eprintf(...) fprintf(stderr, __VA_ARGS__)

// BEGIN ALGO
// calculates cyclic shifts
const int MAXL = 200002; // without \0
const int ALPHA = 128;  // |alphabet|, should be<=MAXL
char s[MAXL + 1];
/*BOXNEXT*/
int *sarr; // cyclic shifts, in sorted order
int sapos[MAXL]; // sarr ^ -1

#define nlen(x) ((x) >= len ? (x) - len : (x))
void build() {
  int len = strlen(s);
  
  static int _sarr1[MAXL], _sarr2[MAXL];
  static int _cols1[MAXL], _cols2[MAXL];
  static int cnt[MAXL];
  
  sarr = _sarr1;
  int *nsarr = _sarr2;
  int *cols = _cols1, *ncols = _cols2;
  
  for (int i = len - 1; i >= 0; i--) {
    cols[i] = s[i];
    sarr[i] = i;
  }
  
  int colcnt = ALPHA;
  for (int clen = 0; clen < len;
           clen = clen ? (clen << 1) : 1) {
    memset(cnt, 0, sizeof(cnt[0]) * colcnt);
    for (int i = 0; i < len; i++)
      cnt[cols[i]]++;
    for (int i = 1; i < colcnt; i++)
      cnt[i] += cnt[i - 1];
    
    for (int i = len - 1; i >= 0; i--) {
      int a = nlen(sarr[i] + len - clen);
      nsarr[--cnt[cols[a]]] = a;
    }
    
    colcnt = 0;
    for (int i = 0; i < len; i++) {
      int a = nsarr[i], b = nlen(a + clen); /*BOXNEXT*/
      if (i == 0 || cols[a] != cols[nsarr[i - 1]] || cols[b] != cols[nlen(nsarr[i - 1] + clen)])
        colcnt++;
      ncols[a] = colcnt - 1;
    }
    swap(cols, ncols);
    swap(sarr, nsarr);
    
    if (colcnt == len) break;
  }
  if (colcnt < len) { // required, otherwise equal shifts will be sorted by (i+wtf)
    memset(cnt, 0, sizeof(cnt[0]) * colcnt);
    for (int i = 0; i < len; i++)
      cnt[cols[i]]++;
    for (int i = 1; i < colcnt; i++)
      cnt[i] += cnt[i - 1];
    for (int i = len - 1; i >= 0; i--)
      sarr[--cnt[cols[i]]] = i;
  }
  for (int i = 0; i < len; i++)
    sapos[sarr[i]] = i;
}

int lcps[MAXL]; // lcp[i] = lcp(sarr[i], sarr[i+1])

void build_lcp() {
  int len = strlen(s);
  if (len <= 1) return;
  
  int i = nlen(sarr[len - 1] + 1), clcp = 0;
  for (int step = 0; step < len - 1; step++) {
    int j = sarr[sapos[i] + 1];
    
    while (clcp < len) {
      int pos1 = nlen(i + clcp),
          pos2 = nlen(j + clcp);
      if (s[pos1] != s[pos2]) break;
      clcp++;
    }
    lcps[sapos[i]] = clcp;
    
    if (j == len - 1) clcp = 0;
    if (++i >= len) i = 0;
    clcp = max(clcp - 1, 0);
  }
}
#undef nlen
// END ALGO

char a[MAXL + 1], b[MAXL + 1];
char out[MAXL + 1];

int main() {
  #ifdef DEBUG
  freopen(".in", "r", stdin);
  freopen(".out", "w", stdout);
  #endif

  int n;
  while (scanf("%d%s%s", &n, a, b) == 3) {
    for (int i = 0; i < n; i++)
      s[i] = a[i];
    s[n] = '#';
    for (int i = 0; i < n; i++)
      s[n + i + 1] = b[i];
    s[2 * n + 1] = '$';
    s[2 * n + 2] = 0;

    build();
    build_lcp();

    int ans = 0;
    const char *bs = "";
    set<int> las[2];

    for (int i = 0; i < 2 * n + 2; i++) {
      int cur = sarr[i];
      int cty = cur > n;
      las[cty].clear();
      if (!las[!cty].empty()) {
        int cans = *las[!cty].begin();
        if (cans > ans) {
          ans = cans;
          if (cty == 0) bs = a + cur;
          else bs = b + (cur - n - 1);
        }
      }

      las[0].insert(lcps[i]);
      las[1].insert(lcps[i]);
    }

/*    #ifdef DEBUG
    int len = strlen(s);
    for (int i = 0; s[i]; i++) {
      printf("i=%2d: %2d %2d, ", i, sarr[i], lcps[i]);
      for (int i2 = 0; i2 < len; i2++)
        printf("%c", s[(sarr[i] + i2) % len]);
      printf("\n");
    }
    #endif*/
    for (int i = 0; i < ans; i++)
      out[i] = bs[i];
    out[ans] = 0;
    printf("%s\n", out);
  }
  eprintf("%d\n", clock());
  return 0;
}
