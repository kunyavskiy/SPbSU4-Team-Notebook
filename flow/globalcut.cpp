#include <cassert>
#include <cstdlib>
#include <iostream>
using namespace std;

typedef long long ll;
const int maxn = 100;
const int INF = 1e9;
// BEGIN ALGO
// g[i][j] = g[j][i] is sum of edges between i and j
// ans is value of mincut

ll g[maxn][maxn], w[maxn];
ll ans;
bool ex[maxn], inA[maxn];
int n;

void iterate(int curN){
  int prev = -1;
  memset(inA, 0, sizeof(inA));
  memset(w, 0, sizeof(w));
  for (int it = 0; it < curN; it++) {
    int best = -1;
    for (int i = 0; i < n; i++)
      if (ex[i] && !inA[i]
          && (best == -1 || w[i] > w[best]))
         best = i;
    assert(best != -1);
    if (it == curN - 1) {
      ans = min(ans, w[best]);
      for (int i = 0; i < n; i++){
        g[i][prev] += g[best][i];
        g[prev][i]  = g[i][prev];
      }
      ex[best] = false;
    } else {
      inA[best] = true;
      for (int i = 0; i < n; i++)
        w[i] += g[best][i];
        prev = best;
    }
  }
}

void solve(){
  ans = INF;
  for (int i = n; i > 1; i--)
    iterate(i);
  cout << ans << endl;
}
//END ALGO

int main(){
	solve();
}