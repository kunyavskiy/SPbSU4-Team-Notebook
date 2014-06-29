//#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <cstring>
#include <string>
#include <cmath>
#include <cassert>
#include <ctime>
#include <algorithm>
#include <sstream>
#include <list>
#include <queue>
#include <deque>
#include <stack>
#include <cstdlib>
#include <cstdio>
#include <iterator>
#include <functional>
#include <bitset>
#define mp make_pair
#define pb push_back

#ifdef LOCAL
#define eprintf(...) fprintf(stderr,__VA_ARGS__)
#else
#define eprintf(...)
#endif

#define TIMESTAMP(x) eprintf("["#x"] Time : %.3lf s.\n", clock()*1.0/CLOCKS_PER_SEC)
#define TIMESTAMPf(x,...) eprintf("[" x "] Time : %.3lf s.\n", __VA_ARGS__, clock()*1.0/CLOCKS_PER_SEC)

#if ( _WIN32 || __WIN32__ )
  #define LLD "%I64d"
#else
  #define LLD "%lld"
#endif

using namespace std;

#define TASKNAME "shift"

typedef long long ll;
typedef long double ld;

const int MAXN = 110000;

// BEGIN ALGO
char s[MAXN*2 + 1];  //input string, will be twiced by min_shift

// returns vector<int> of begins of prime strings
// string should have zero simbol at end.
vector<int> duval(){
  vector<int> res;
  int len = strlen(s) + 1; // zero used here
  int start = 0, mid = 1, cur = 0;
  res.pb(0);
  for (int i = 0; i < len; i++){
    if (s[i] == s[cur]){
      cur++;
      if (cur == mid) cur = start;
    } else if (s[i] > s[cur]){
      mid = i+1;
      cur = start;
    } else if (s[i] < s[cur]){
      int temp = mid - start;
      while (start + temp <= i){
        start += temp;
        res.pb(start);
      }
      i = cur = start;
      mid = start + 1;
    }
    else
      assert(false);
  }
  return res;
}

int min_shift(){
  int len = strlen(s);
  memcpy(s+len,s,sizeof(char)*len);
  s[2*len] = 0;
  vector<int> v = duval();
  s[len] = 0;
  for (int i = 0; i < (int)v.size(); i++)
    if (i == (int)v.size()-1 || v[i+1] >= len)
      return v[i];
  assert(false);
}
//END ALGO


int main(){
  freopen(TASKNAME".in","r",stdin);
  freopen(TASKNAME".out","w",stdout);

  gets(s);
  int n;
  sscanf(s,"%d",&n);

  for (int i = 0; i < n; i++){
    gets(s);
    printf("%d\n",min_shift());
  }



  TIMESTAMP(end);
  return 0;
}