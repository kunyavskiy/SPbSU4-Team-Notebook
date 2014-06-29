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
#define TASKNAME "nenokku"

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

typedef long double ld;
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

const int MAXN = 1100000;

// BEGIN ALGO
struct node{
  node* next[26], *lnk, *p;
  int c,L,R;
  node(){
    c = L = R = -1;
    memset(next,0,sizeof(next));
    lnk = p = 0;
  }
  node(node* p,int c,int L,int R):
          p(p), c(c), L(L), R(R){
    memset(next,0,sizeof(next));
    lnk = 0;
  }
};

struct position{
  node* to;
  int pos;
  void next(int val);
  bool can(int val);
  position(){}
  position(node* to,int pos):to(to),pos(pos){}
};

node* root;
char s[MAXN];
int len;
position last;

void position::next(int val){
  if (!to) return;
  if (to->R == pos){
    to = to->next[val];
    if (!to) return;
    pos = to->L;
  }
  if (s[pos] - 'a' != val) to = 0;
  pos++;
}

bool position::can(int val){
  if (!to) return false;
  if (to->R == pos){
    if (!to->next[val]) return false;
    return true;
  }
  if (s[pos] -'a' != val) return false;
  return true;
}

node* split(node* v,int id){
  assert(v && v->L < id && id <= v->R);
  if (id == v->R) return v;
  node* nv = new node(v->p,v->c,v->L,id);
  v->p->next[nv->c] = nv;
  v->L = id;
  v->p = nv;
  nv->next[s[id] - 'a'] = v;
  v->c = s[id] - 'a';
  return nv;
}

node* addleaf(node *v,int id){
  assert(!v->next[s[id]-'a']);
  node* nv = new node(v,s[id]-'a',id,inf);
  v->next[s[id]-'a'] = nv;
  return nv;
}

position go(node* v,int L,int R){
  while (true){
  	int c = s[L] - 'a';
    assert(v->next[c]);
    int npos = v->next[c]->L + R - L;
    if (v->next[c]->R >= npos)
      return position(v->next[c], npos);
    v = v->next[c];
    L += v->R - v->L;
  }
}

void add(){
  node *prevlast = 0;
  while (last.pos == len && last.to != root)
    last = go(last.to->p->lnk,last.to->L,last.pos);
  while (!last.can(s[len]-'a')) {
    node* lastnode = split(last.to, last.pos);
    addleaf(lastnode,len);
    if (prevlast) prevlast->lnk = lastnode;
    prevlast = lastnode;
    if (lastnode == root) {
      last = position(root,0);
      break;
    }
    assert(lastnode->p->lnk);
    last = go(lastnode->p->lnk,lastnode->L,last.pos);
  }
  if (prevlast && !prevlast->lnk){
    assert(last.to && last.pos == last.to->R);
    prevlast->lnk = last.to;
  }
  assert(last.can(s[len]-'a'));
  last.next(s[len]-'a');
  len++;
}

void init(){
  root = new node;
  last = position(root,0);
  node* rroot = new node;
  for (int i = 0; i < 26; i++)
    rroot->next[i] = root;
  rroot->lnk = root->lnk = root->p = rroot;
  root->L = -1, root->R = 0;
}
// END ALGO

char buf[MAXN];

int main(){
  freopen(TASKNAME".in","r",stdin);
  freopen(TASKNAME".out","w",stdout);
  init();
  char c;
  while (scanf(" %c",&c) == 1){
    if (c == '?'){
      position pos;
      pos.to = root;
      pos.pos = 0;
      scanf(" %s",buf);
      for (int i = 0; buf[i]; i++){
        buf[i] = tolower(buf[i]);
        pos.next(buf[i]-'a');
      }
      if (!pos.to) printf("NO\n");
      else printf("YES\n");
    }
    else {
      scanf(" %s",s+len);
      while (s[len]){
        s[len] = tolower(s[len]);
        add();
     }
    }
  }

  TIMESTAMP(end);
  return 0;
}