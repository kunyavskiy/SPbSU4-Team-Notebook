#include <cstdio>
#include <cstring>

// http://acm.timus.ru/problem.aspx?space=1&num=1590

int ans = 0;

// BEGIN ALGO
const int MAXL = 100000; // maximal length, without trailing \0
char s[MAXL + 1]; // used in build() only.

struct node{
  node* next[26], *link;
  int len; bool term;
  node(){
    memset(next,0,sizeof(next));
    len = 0; link = 0; term = 0;
  }
};
node mem[2 * MAXL];
int memptr; node *root, *last;

void add(int c){
  node* nlast = &mem[memptr++];
  *nlast = node();
  nlast->len = last->len + 1;
  node* p = last;
  for (;p && p->next[c] == 0; p = p->link)
    p->next[c] = nlast;
  if (!p){
    nlast->link = root;
  } else {
    node* q = p->next[c];
    if (q->len == p->len + 1){
      nlast->link = q;
    } else {
      node* clone = &mem[memptr++];
      memcpy(clone, q, sizeof(node));
      clone->link = q->link;
      clone->len = p->len + 1;
      q->link = nlast->link = clone;
      for (;p && p->next[c] == q; p = p->link)
        p->next[c] = clone;
    }
  }
  last = nlast;
}
void build(){
  memptr = 0;
  last = root = &mem[memptr++];
  *root = node();
  for (int i = 0; s[i]; i++)
  { // NOT ALGO
    add(s[i]-'a');
    ans += last->len - last->link->len; // NOT ALGO
  } // NOT ALGO
  for (node* it = last; it; it=it->link)
    it->term = true;
}
/* // NOT ALGO
// build suffix tree of \emph{reversed} string
void BuildTree(){
  dfsany(root); // find length to any terminal
  for (int i = 1; i < memptr; i++)
    g[mem[i].lnk - mem].pb(Edge(
        i,mem[i].any + mem[i].link->len,
        mem[i].any + mem[i].len)); // to, l, r
}
*/ // NOT ALGO
// END ALGO

int main() {
  #ifdef DEBUG
  freopen(".in", "r", stdin);
  freopen(".out", "w", stdout);
  #endif
  while (gets(s)) {
    ans = 0;
    build();
    printf("%d\n", ans);
  }
  return 0;
}
