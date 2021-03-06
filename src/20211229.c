#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <strings.h>
#include <stdbool.h>
#include "struct.h"
#include "chash-p.c"

#define O 6859 // 1331 //2197,4913,6859
#define K 5
#define P 37

//using namespa
// sagemath上での原始多項式
unsigned long long pp[4][4] = {{0, 0, 9, 2}, {0, 0, 11, 2}, {0, 0, 16, 3}, {0, 0, 15, 2}};
// {0,0,9,2}, {1,0,11,2}, {1,0,16,3}, {1,0,15,2};
// GF(11^3,13^3,17^3,19^3)
// unsigned long long ff[2][7]={{1,0,0,0,0,2,0,2},{0,0,1,0,0,0,1,2}}; //GF(3^7,5^5)

unsigned long long gf[O] = {0}, fg[O] = {0};
// int N =0,M=0;
unsigned long long c[K + 1] = {0};

// OP型からベクトル型への変換
vec o2v(OP f)
{
  vec a = {0};
  int i;

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0)
      a.x[f.t[i].n] = f.t[i].a;
  }

  return a;
}

//ベクトル型からOP型への変換
OP v2o(vec a)
{
  int i, j = 0;
  OP f = {0};

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
    {
      f.t[i].n = i;
      f.t[i].a = a.x[i];
    }
  }

  return f;
}

void op_print_raw(const OP f)
{
  puts("op_print_raw:");
  for (int i = 0; i < DEG; i++)
  {
    if (f.t[i].a > 0)
      printf("[%d] %ux^%u\n", i, f.t[i].a, f.t[i].n);
  }
}

//多項式の次数(default)
int deg(vec a)
{
  int i, n = 0;

  //#pragma omp parallel for
  for (i = 0; i < DEG; i++)
  {
    if (a.x[i] > 0)
      n = i;
  }

  return n;
}

//配列からベクトル表現の多項式へ変換する
vec Setvec(int n)
{
  int i;
  vec v = {0};

  for (i = 0; i < n; i++)
  {
    v.x[n - 1 - i] = c[i];
  }

  return v;
}

//配列の値を係数として多項式に設定する
OP setpol(unsigned long long f[], int n)
{
  OP g;
  vec a;
  int i;

  memset(c, 0, sizeof(c));
  // memcpy (c, f, n);
  for (i = 0; i < n; i++)
  {
    c[i] = f[i];
    // printf("%d,",f[i]);
  }
  // exit(1);
  a = Setvec(n);

  g = v2o(a);

  return g;
}

//多項式を表示する(default)
void printpol(vec a)
{
  int i, n;

  n = deg(a);

  // printf ("baka\n");
  assert(("baka\n", n >= 0));

  for (i = n; i > -1; i--)
  {
    if (a.x[i] > 0)
    {
      printf("%d", a.x[i]);
      // if (i > 0)
      printf("x^%d", i);
      // if (i > 0)
      printf("+");
    }
  }
  //  printf("\n");

  return;
}

//多項式を表示する(default)
void printpoln(vec a)
{
  int i, n;

  n = deg(a);

  // printf ("baka\n");
  assert(("baka\n", n >= 0));

  for (i = n; i > -1; i--)
  {
    if (a.x[i] > 0)
    {
      printf("%d", a.x[i]);
      // if (i > 0)
      printf("x^%d", i);
      // if (i > 0)
      printf("+");
    }
  }
  printf("\n");

  return;
}

//多項式の代入値
unsigned long long
xtrace(OP f, unsigned long long x)
{
  int i, d;
  unsigned long long u = 0, v = 1;

  d = deg(o2v(f));
  // printpol(o2v(f));
  // printf(" =ff\n");

  for (i = 0; i < d + 1; i++)
  {
    v = 0;
    if (f.t[i].a > 0)
    {
      v = 1;
      for (int j = 0; j < f.t[i].n; j++)
      {
        v = (v * x);
      }
      v = (v * f.t[i].a) % O;

      // printf("\nv=%d",v);
    }
    u = (u + v);
  }
  // printf("u=%d\n",u%O);

  return u % O;
}

void makefg(int n)
{
  unsigned long long i, j, count = 0;

  for (i = 0; i < O; i++)
  {

    for (j = 0; j < O; j++)
    {
      if (gf[i] == j)
      {
        fg[j] = i;
        count++;
      }
    }
  }
  printf("unsigned long long fg[%d]={", O);
  for (i = 0; i < O; i++)
    printf("%d,", fg[i]);
  printf("};\n");
  printf("count=%d\n", count);
  exit(1);

  return;
}

//多項式の次数(degのOP型)
int odeg(OP f)
{
  int i, j = 0, k;

  // k=terms(f);
  for (i = 0; i < 512; i++)
  {
    if (j < f.t[i].n && f.t[i].a > 0)
      j = f.t[i].n;
  }

  return j;
}

OP minus(OP f)
{
  unsigned int i, j;

  j = deg(o2v(f));
  for (i = 0; i < j + 1; i++)
  {
    if (f.t[i].a > 0)
      f.t[i].a = P - f.t[i].a;
  }

  return f;
}

//リーディグタームを抽出(default)
oterm LT(OP f)
{
  int i, k;
  oterm t = {0};

  // k = deg (o2v (f));
  for (i = 0; i < DEG; i++)
  {
    // printf("a=%d %d\n",f.t[i].a,f.t[i].n);
    if (f.t[i].a > 0)
    {
      t.n = f.t[i].n;
      t.a = f.t[i].a % P;
    }
  }

  return t;
}

// OP型を正規化する
OP conv(OP f)
{
  vec v = {0};
  OP g = {0};

  v = o2v(f);
  g = v2o(v);

  return g;
}

// 20200816:正規化したいところだがうまく行かない
//多項式の足し算
OP oadd(OP f, OP g)
{
  vec a = {0}, b = {0}, c = {0};
  int i, j, k, l = 0;
  OP h = {0}, f2 = {0}, g2 = {0};

  // for(i=0;i<257;i++)
  //  printf("%d %d %d %d %d\n",i,f.t[i].a,f.t[i].n,g.t[i].a,g.t[i].n);

  //  exit(1);
  f = conv(f);
  g = conv(g);

  a = o2v(f);
  // exit(1);
  b = o2v(g);

  j = deg(o2v(f));
  l = deg(o2v(g));
  printpol(o2v(f));
  printf(" =f_& %d\n", j);
  printpol(o2v(g));
  printf(" =g_& %d\n", l);
  if(LT(f).a==0 && LT(g).a==0)
   return h;

  if (j >= l)
  {
    k = j + 1;
  }
  else
  {

    k = l + 1;
  }
  // for(i=0;i<k;i++)
  // printf("%d %d\n",i,b.x[i]);
  //   exit(1);

  for (i = 0; i < k; i++)
  {
    // if(a.x[i]>b.x[i])
    c.x[i] = (a.x[i] + b.x[i]) % P;
    
  }
  //
  h = v2o(c);
  printpol(o2v(h));
  printf(" =====in oadd\n");

  return h;
}

//多項式を項ずつ掛ける
OP oterml(OP f, oterm t)
{

  // assert (op_verify (f));
  int i, k, j;
  OP h = {0};
  vec test;
  unsigned long long n;

  // f=conv(f);
  k = odeg(f);
  j = 0;
  for (i = 0; i < k + 1; i++)
  {
    h.t[i].n = f.t[i].n + t.n;
    h.t[i].a = (f.t[i].a * t.a) % P;
  }

  // h=conv(h);
  // assert (op_verify (h));
  return h;
}

//多項式の掛け算
OP omul(OP f, OP g)
{
  f = conv(f);
  g = conv(g);
  // assert (op_verify (f));
  // assert (op_verify (g));
  int i, count = 0, k, l;
  oterm t = {0};
  OP h = {0}, e = {0}, r = {0};
  vec c = {0};

  k = odeg(f);
  l = odeg(g);
  if(LT(f).a==0 || LT(g).a==0)
    return h;

  if (l > k)
  {
    k = l;
  }

  for (i = 0; i < k + 1; i++)
  {
    t = g.t[i];
    e = oterml(f, t);
    h = oadd(h, e);
  }
  if(LT(h).a==0){
    printf("wh==0\n");
    exit(1);
  }
  // assert (op_verify (h));
  return h;
}

OP osub(OP f, OP g)
{
  vec a = {0}, b = {0}, d = {0};
  int i, k, l, m;
  OP ans = {0};

  a = o2v(f);
  b = o2v(g);
  l = deg(a);
  m = deg(b);
  if (l >= m)
  {
    k = l;
  }
  else
  {
    k = m;
  }
  for (i = 0; i < k + 1; i++)
  {
    if (a.x[i] >= b.x[i])
    {
      d.x[i] = a.x[i] - b.x[i];
    }
    else
    {
      d.x[i] = (P + (a.x[i] - b.x[i]))%P ;
    }
    
    printf("%d - %d = %d\n",a.x[i],b.x[i],d.x[i]);
    //if(d.x[i]<0){
    //  printf("%d\n",d.x[i]);
    //  d.x[i]+=P;
    //}
    
  }
//exit(1);
  ans = v2o(d);

  return ans;
}

OP confer(OP f, int a)
{
  vec r;
  int n, i;
  OP g;

  r = o2v(f);
  n = deg(r);
  for (i = 0; i < n + 1; i++)
    r.x[i] = (r.x[i] * a) % P;
  g = v2o(r);

  return g;
}

int oequ(OP f, OP g)
{
  vec v, x;
  int i, flg = 0;

  v = o2v(f);
  x = o2v(g);
  for (i = 0; i < 512; i++)
  {
    if (v.x[i] != x.x[i])
      return -1;
  }

  return 0;
}

void mkmf()
{
  int i, j, k, count = 0;
  OP f = {0}, g = {0}, h = {0}, w = {0}, s = {0}, u = {0};
  vec b = {0}, a = {0}, d = {0}, t = {0}, v = {0};
  oterm o;
  unsigned long long ccp[4] = {0};

  if (O == 1331)
  {
    for (i = 0; i < K + 1; i++)
      ccp[i] = pp[0][i];
  }
  if (O == 2197)
  {
    for (i = 0; i < K + 1; i++)
      ccp[i] = pp[1][i];
  }
  if (O == 4913)
  {
    for (i = 0; i < K + 1; i++)
      ccp[i] = pp[2][i];
  }
  if (O == 6859)
  {
    for (i = 0; i < K + 1; i++)
      ccp[i] = pp[3][i];
  }

  g = setpol(ccp, 4);
  // b.x[0]=2;
  // b.x[1]=9;
  a.x[1] = 1;
  v.x[3] = 1;
  u = v2o(v);
  d.x[3] = 1;
  printpol(o2v(g));
  printf(" =g\n");
  // exit(1);

  // g=v2o(b);
  s = v2o(a);
  // s.t[1].a=1;
  // s.t[1].n=1;
  // gf[12]=P;
  // gf[13]=P*P;

  w = g;
  printpol(o2v(w));
  printf(" =w\n");
  printpol(o2v(g));
  printf(" =g\n");
  printpol(o2v(s));
  printf(" =s\n");

  printf("\n");
  // for(i=0;i<P;i++)
  gf[0] = 0;
  gf[1] = 1;
  gf[2] = P;
  gf[3] = P * P;
  gf[4] = xtrace(g, P);
  printf("\naa=%d\n", gf[P]);
  // exit(1);
  // w=omul(w,s);
  // gf[12]=1111;
  count = 4 + 1;
  while (1)
  {
    g = omul(g, s);
    printpol(o2v(g));
    printf(" =g\n\n");
    printf(" papaya\n");

    // exit(1);

    o = LT(g);
    memset(d.x, 0, sizeof(d));
    if (o.n == K)
    {
      d.x[o.n] = o.a;
      h = v2o(d);
      g = osub(g, h);
      f = confer(w, o.a);
      g = oadd(g, f);
      printpol(o2v(f));
      printf("\n");
    }

    gf[count] = xtrace(g, P);
    printf("count=%d %d ", count, gf[count]);
    printpol(o2v(g));
    printf(" =gg\n\n");
    if (gf[count] == 1)
    {
      printf("count=%d\n", count);
      break;
    }

    count++;
  }

  printf("unsigned long long gf[%d]={", O);
  for (i = 0; i < O; i++)
    printf("%d,", gf[i]);
  printf("};");
  printf("\n");

  // exit(1);
}

// nを法とする逆数
unsigned int inv(unsigned int a, unsigned int n)
{
  unsigned int d, x, s, q, r, t, gcd;
  d = n;
  x = 0;
  s = 1;
  while (a != 0)
  {
    q = d / a;
    r = d % a;
    d = a;
    a = r;
    t = (x - q * s);
    x = s;
    s = t;
  }
  gcd = d;

  return ((x + n) % (n / d));
}


// aに何をかけたらbになるか
unsigned long long
equ(unsigned long long a, unsigned long long b)
{
  int i=inv(a,P);

return i*b;

}


//多項式を単行式で割る
oterm LTdiv(OP f, oterm t)
{
  oterm tt = {0}, s = {
                      0};

  tt = LT(f);
  if (tt.n < t.n)
  {
    s.n = 0;
    s.a = 0;
  }
  else if (tt.n == t.n)
  {
    s.n = 0;
    s.a = equ(t.a, tt.a);
  }
  else if (tt.n > t.n)
  {
    s.n = tt.n - t.n;
    s.a = equ(t.a, tt.a);
    // printf("%u\n",s.a);
  }
  else if (t.n == 0 && t.a > 0)
  {
    s.a = (inv(t.a,P)*tt.a)%P; //(tt.a * inv(t.a, P)) % P;
    s.n = tt.n;
  }
  else
  {
    printf("debug in LTdiv\n");
    exit(1);
  }

  return s;
}

//多項式の剰余を取る
OP omod(OP f, OP g)
{
  int i = 0, j, n, k;
  OP h = {0}, o = {0}, e = {0};
  oterm a, b = {0}, c = {0};

  n = LT(g).n;
  if (LT(f).a == 0)
  {
    return h;
  }
  if(LT(g).a==0){
    printf("g de dib 0!\n");
    exit(1);

  }
  //  assert (("baka^\n", LT (f).n != 0));

  //  assert (("baka(A)\n", LT (g).n != 0));

  if (LT(f).n < LT(g).n)
  {
    //    exit(1);
    return f;
  }

  // printf ("in omod\n");
  // exit(1);

  k = LT(g).n;
  b = LT(g);
  OP ll;

  //assert(("double baka\n", b.a > 0 && b.n > 0));
  while(LT(g).a != 0)
  {
    c = LTdiv(f, b);
    h = oterml(g, c);
    printpol(o2v(f));
    printf("======f_before_omod\n");
    printpol(o2v(h));
    printf("======h_before_omod\n");
    f = osub(f, (h));
    printpol(o2v((h)));
    printf(" =====h_minus_omod\n");
    printpol(o2v(f));
    printf(" =====f_after_omod\n");
    // exit(1);
    if (LT(f).a == 0 )
    {
      //      printf("blake1\n");
      break;
    }

    if (c.a == 0 || b.a==0)
      break;
  }
  printpol(o2v(f));
  printf("\n");
  // exit(1);

  return f;
}

//項の数
int terms(OP f)
{
  int i, count = 0;

  for (i = 0; i < DEG; i++)
    if (f.t[i].a > 0)
      count++;

  return count;
}

//モニック多項式にする
OP coeff(OP f)
{
  int i, j, k;
  vec a, b;
  oterm t;

  t = LT(f);
  // f = conv(f);
  k = odeg((f)) + 1;
  for (i = 0; i < k; i++)
    f.t[i].a = (f.t[i].a * inv(t.a, P)) % P;

  return f;
}


OP cdiv(int a,OP f){
  vec v;
  int i,l,k;

  v=o2v(f);
  l=odeg(f);
  k=inv(a,P);
  printpol(o2v(f));
  printf(" ==kokko %d\n",a);
  //a=equ(a,LT(f).a);
  for(i=0;i<l+1;i++)
  v.x[i]=(k*v.x[i])%P;
  printpol(v);
// exit(1);
  f=v2o(v);

  return f;
}

//多項式の商を取る
OP odiv(OP f, OP g)
{

  f = conv(f);
  g = conv(g);
  // assert (op_verify (f));
  // assert (op_verify (g));
  int i = 0, j, n, k;
  OP h = {0}, e = {0}, tt = {0}, o = {0};
  oterm a, b = {0}, c = {0};

  printpol(o2v(f));
  printf("\n");
  printpol(o2v(g));
  printf("\n");
  // exit(1);

  if (LT(f).a == 0)
  {
    return h;
  }
  if(LT(g).a==0){
    printf("baka^\n");
    exit(1);
  }
   

  k = 0;//odeg(f) - odeg(g);
  b = LT(g);
  
  if (b.a > 0 && b.n == 0){
    e=cdiv(b.a,f);
    printpol(o2v(e));
    //exit(1);
    return e;
  }
  
  if (b.a == 0 && b.n == 0)
  {
    printf("baka in odiv\n");
    exit(1);
  }
  OP null = {0};
  if (odeg((f)) < odeg((g)))
  {
    return null;
  }

  i = 0;
  k=0;
  while (LT(f).a != 0 || LT(g).a != 0)
  {
    c = LTdiv(f, b);
    c.a = c.a % P;
    assert(c.n < DEG);
    tt.t[k] = c;
    k++;

    printf("%d", c.a);
    printf(" ccccccccccccccccc\n");
    printpol(o2v(g));
    printf(" ===before g in_odiv\n");
    printpol(o2v(f));
    printf(" ===before f in_odiv\n");
    h = oterml(g, c);
    f = osub(f, (h));
    printpol(o2v(h));
    printf(" ===h in_odiv\n");
    printpol(o2v(g));
    printf(" ===g in_odiv\n");
    printpol(o2v(f));
    printf(" ===f in_odiv\n");
    if (LT(f).a == 0 || LT(g).a == 0)
    {
      printf("blake2\n");
      break;
    }
    if (oequ(f, g) == 0)
    {
      printpol(o2v(tt));
      printf("\n");
      break;
      // exit(1);
    }
    if (c.a == 0)
      break;
  }

  // tt は逆順に入ってるので入れ替える
  OP ret = {0};
  
  int tt_terms = terms (tt);
  for (i = 0; i < tt_terms; i++)
    {
      ret.t[i] = tt.t[tt_terms - i - 1];
    }
  
  ret = conv(ret);
  printpol(o2v(ret));
  printf("  return\n");
  // exit(1);

  // assert (op_verify (ret));
  return ret;
}


OP scr(unsigned long long d, OP f)
{
  int i, n;
  vec v = {0};

  n = deg(o2v(f));
  v = o2v(f);
  for (i = 0; i < n+1 ; i++){
    printf("v[%d]=%d\n",i,v.x[i]);
    v.x[i] = (v.x[i] * d) % P;
  }
  f = v2o(v);

  return f;
}

OP monique(OP f){
  int e1;
  e1 = inv(LT(f).a, P);
  printf("e=%d\n",e1);
  f=scr(e1,f);

return f;
}

EX monic(EX X){
  int e1;
  e1=inv(LT(X.d).a,P);
  X.d=scr(e1,X.d);
  X.h=scr(e1,X.h);
  X.u=scr(e1,X.u);
  X.v=scr(e1,X.v);

  return X;
}

//拡張ユークリッドアルゴリズム
EX xgcd(OP f, OP g)
{
  OP h[10] = {0}, ww[10] = {0}, *v, *u,T={0};
  oterm a, b;
  int i = 0, j,  flg = 0,k;
  EX e = {0}, ee = {0};
/*
if(odeg(f)<odeg(g)){
T=f;
f=g;
g=T;
}
*/
  v = (OP *)malloc(sizeof(OP) * (DEG));
  u = (OP *)malloc(sizeof(OP) * (DEG));
  memset(v, 0, sizeof(OP) * DEG);
  memset(u, 0, sizeof(OP) * DEG);

  u[0].t[0].a = 1;
  u[0].t[0].n = 0;
  u[1].t[0].a = 0;
  u[1].t[0].n = 0;
  u[2].t[0].a = 1;
  u[2].t[0].n = 0;

  v[0].t[0].a = 0;
  v[0].t[0].n = 0;
  v[1].t[0].a = 1;
  v[1].t[0].n = 0;

  printpol(o2v(f));
  printf(" f===============\n");
  printpol(o2v(g));
  printf(" s===============\n");
  // exit(1);
  if(LT(g).a==0){
    printf("g ==0\n");
    exit(1);
  }

k = odeg(g);
if(k==0 && LT(g).a!=0){
  printf("use cdiv\n");
  //f=cdiv(LT(g).a,f);
  //printpoln(o2v(g));
  //exit(1);
}
 i = 1;
  while (1)
//  for (i = 1; i < 1+ k+1; i++)
  {
    printpol(o2v(f));
    printf(" fffffffffffffffff\n");

    printpol(o2v(g));
    printf(" ggggggggggggggggg\n");

    h[i] = omod(f, g);
    printpol(o2v(h[i]));
    printf(" %d hhhhhhhhhhhhhh\n", i);

    ww[i] = odiv(f, g);
    printpol(o2v(ww[i]));
    printf(" %d wwwwwwwwwwwwww\n", i);
  
    v[i + 1] = osub(v[i - 1], (omul(ww[i], v[i])));
    printpol(o2v(v[i]));
    printf(" vvvvvvvvvvvvvvvv[%d]\n", i);
    printpol(o2v(v[i + 1]));
    printf(" ==vv[%d]\n", i + 1);
    // exit(1);
    u[i + 1] = osub(u[i - 1], (omul(ww[i], u[i])));
    printpol(o2v(u[i]));
    printf(" %d uuuuuuuuuuuuuu\n", i);
    f = g;
    g = h[i];

if(LT(g).a==0){
printpol(o2v(g));
printf("ggggggggggggg\n");
}
if(LT(g).a==0)
  break;
  i++;
  }
//f=g;
//g=h[i];
  if(LT(g).a!=0)
    ww[i] = odiv(f, g);
   printpol(o2v(ww[i]));
  printf(" %d wwwwwwwwwwwwww\n", i);
  printpol(o2v(f));
  printf(" =========fvi@\n");
  printpol(o2v(v[i]));
  printf(" =========vi@\n");
  printpol(o2v(u[i]));
  printf(" =========uvi@\n");
  printpol(o2v(h[i]));
  printf(" =========mvi@\n");
  printpol(o2v(ww[i]));
  printf(" =========wvi@\n");
  printpol(o2v(g));
  printf(" =========gvi@\n");
//exit(1);
  

  e.d = f;
  printpoln(o2v(e.d));
  e.u = u[i];
  e.v = v[i];
  e.h = ww[i];

  free(v);
  free(u);

  printf("end of fnc\n");
  //exit(1);
  return e;
  //  wait ();
}


OP qinv(OP uu1, OP uu2)
{
  EX tt,V;
  OP v;
  int e;

//return cdiv(LT(g).a,g);
  tt = xgcd(uu1, uu2);
  printpol(o2v(tt.v));
  printf("\n");
  printpol(o2v(tt.u));
  printf("\n");
  printpol(o2v(tt.d));
  printf("\n");
  printpol(o2v(tt.h));
  printf("\n");
  //exit(1);
  v=scr(LT(tt.d).a,tt.u);
  printpol(o2v(v));
  printf(" in qinv's v\n");

  //exit(1);
if(LT(omod(omul(v,uu1),uu2)).a==0){
  return v;
}
printf("can't\n");
//exit(1);
}


int chkdiv(Div d, OP f)
{
  OP t;

  t=omod(osub(omul(d.v, d.v), (f)), d.u);
  printpol(o2v(t));
  printf(" 00000000000 t\n");
  if (LT(t).a==0){
    return 1;
  }else{
      printpoln(o2v(d.u));
      printpoln(o2v(d.v));
      //exit(1);
      }
  return -1;
}

vec diviser(OP o, OP m)
{
  int t1[2][3], cc[2];
  vec c1 = {0};
  int i, j, k;

  t1[0][0] = o.t[1].a;
  t1[1][0] = o.t[0].a;
  t1[0][1] = m.t[1].a;
  t1[1][1] = m.t[0].a;
  t1[0][2] = 0;
  t1[1][2] = 1;

  if (t1[0][0] == 0)
  {
    for (i = 0; i < 2; i++)
      t1[0][i] = (t1[0][i] + t1[1][i]) % P;
  }
  if (t1[1][1] == 0)
  {
    for (i = 0; i < 2; i++)
      t1[1][i] = (t1[1][i] + t1[0][i]) % P;
  }

  cc[0] = inv(t1[0][0], P);
  printf("%d\n", cc[0]);
  // exit(1);

  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < 3; j++)
      printf("%d,", t1[i][j]);
    printf("\n");
  }
  printf("\n");

  for (j = 0; j < 3; j++)
  {
    t1[0][j] = (t1[0][j] * cc[0]) % P;
    printf("%d,", t1[0][j]);
  }
  printf("\n");
  // exit(1);
  int z;
  z = t1[1][0];
  for (j = 0; j < 3; j++)
  {
    t1[1][j] = t1[1][j] - (t1[0][j] * z) % P;

    if (t1[1][j] < 0)
      t1[1][j] = P + t1[1][j];
    printf("%d,", t1[1][j]);
  }
  printf("\n\n");
  // exit(1);

  cc[1] = inv(t1[1][1], P);
  for (j = 0; j < 3; j++)
  {
    t1[1][j] = (t1[1][j] * cc[1]) % P;
    printf("%d,", t1[1][j]);
  }
  printf("\n\n");
  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < 3; j++)
      printf("%d,", t1[i][j]);
    printf("\n");
  }
  printf("\n\n");
  // exit(1);
  for (i = 0; i < 3; i++)
    printf("b%d,", t1[0][i]);
  printf("\n");
  // exit(1);
  printf("A%d", t1[0][1]);

  int y = t1[0][1];

  for (j = 0; j < 3; j++)
  {
    t1[0][j] = (t1[0][j] - t1[1][j] * y) % P;
    if (t1[0][j] < 0)
      t1[0][j] += P;
    printf("a%d,", t1[1][j] * t1[1][1]);
  }
  printf("\n\n");
  // exit(1);

  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < 3; j++)
      printf("%d,", t1[i][j]);
    printf("\n");
  }
  printf("\n");
  //  exit(1);
  c1.x[0] = t1[0][2];
  c1.x[1] = t1[1][2];
  printf("%d %d\n", c1.x[0], c1.x[1]);
  // exit(1);

  return c1;
}

Div reduce(Div D3,OP f){
  Div D;
  OP ur,vr;

  ur=osub(f,odiv(omul(D3.v,D3.v),D3.u));
  vr=omod(minus(D3.v),ur);
  D.u=ur;
  D.v=vr;

return D;
}

Div g2dbl(Div D,OP ff){
unsigned long long u1,u0,v1,v0,vu,vv,mm1,mm2,mm3,mm4,z1,z2,t1,t2,t3,t4,l2,l1,l3,d,aa,bb,cc,dd,ee,ue1,ue0,uue1,ve1,ve0,uue0,uu1,uu0;
Div D2={0};
vec vx;
unsigned f2,f3;
OP u,v;

u=D.u;
v=D.v;
vx=o2v(u);
u1=vx.x[1];
u0=vx.x[0];
vx=o2v(v);
v1=vx.x[1];
v0=vx.x[0];

vx=o2v(ff);
f2=vx.x[2];
f3=vx.x[3];

uu1=(u1*u1)%P; uu0=(u0*u1)%P;

vv=(v1*v1)%P; vu=(((v1+u1)*(v1+u1))-vv-uu1)%P; mm1=(2*v0-2*vu)%P; mm2=(2*v1*(u0+2*uu1))%P;
mm3=(-2*v1)%P; mm4=(vu+2*v0)%P; z1=(f2+2*uu1*u1+2*uu0-vv)%P; z2=(f3-2*u0+3*uu1)%P;
t1=((mm2-z1)*(z2-mm1))%P; t2=((-z1-mm2)*(z2+mm1))%P; 
t3=((mm4-z1)*(z2-mm3))%P; t4=((-z1-mm4)*(z2+mm3))%P;
l2=(t1-t2)%P; l3=(t3-t4)%P; d=(t3+t4-t1-t2-2*(mm2-mm4)*(mm1+mm3))%P;
aa=inv(d*l3,P); bb=(d*aa)%P; cc=(d*bb)%P; cc=(d*bb)%P; dd=(l2*bb)%P; ee=(l3*l3*aa)%P;
ue1=(2*dd-cc*cc-2*u1)%P; ve0=((dd-u1)*(dd-u1)+2*cc*(v1+cc*u1))%P; uue1=(ue1*ue1)%P; uue0=(ue1*ue0)%P;
ve1=(ee*ve1+v1)%P; ve0=(ee*ve0+v0)%P;
memset(vx.x,0,sizeof(vx.x));
vx.x[2]=1;
vx.x[1]=ue1;
vx.x[0]=ue0;
D2.u=v2o(vx);

memset(vx.x,0,sizeof(vx.x));
vx.x[1]=ve1;
vx.x[0]=ve0;
D2.v=v2o(vx);


return D2;
}

Div cdbl(Div D,OP f){
OP k,s,l,u,v,u1,u2,v2,ut,vt,vc;
Div D2;

u=D.u;
v=D.v;
if(chkdiv(D,f)==-1){
  printf("bakayo\n");
  exit(1);
}

OP t={0};

k=odiv(osub(f,omul(v,v)),u);
s=omod(omul(k,qinv(scr(2,v),u)),u);
l=omul(s,u);
u1=osub(omul(s,s),omod(osub(omul(scr(2,v),s),k),u));
u2=monique(u1);
v2=omod(minus(oadd(l,v)),u2);

printpoln(o2v(u2));
printpoln(o2v(omul(u,u)));
D2.u=u2;
D2.v=v2;

//if(LT(omod(omul(u,u),u2)).a==0 && 
//if(oequ(omod(v2,u),omod(v,u))==0 || LT(omod(osub(omul(v2,v2),f),v2)).a==0){
t=omod(osub(omul(v2, v2), (f)), u2);
if(LT(t).a==0){
printf("ii!\n");
return D2;
}
printf("dame\n");
exit(1);
}


Div cadd(OP ff,OP uu1,OP uu2,OP vv1,OP vv2){
  EX V;
  Div D3,D, null={0};
  vec vx;
  OP e1,e2,d1,s1,s2,s3,c1,c2,d,u;

V=xgcd(uu1,uu2);
V=monic(V);
e1=V.u;
printpol(o2v(V.u));
printf("  e1\n");
e2=V.v;
printpol(o2v(V.v));
printf("  e2\n");
d1=V.d;
printpol(o2v(V.d));
printf("  d1\n");

printpol(o2v(V.h));
printf("  Uh\n");
//exit(1);
V=xgcd(oadd(vv1,vv2),V.d);
V=monic(V);
c2=V.u;
printpol(o2v(V.u));
printf("  c2\n");
c1=V.v;
printpol(o2v(V.v));
printf("  c1\n");
d=V.d;
printpol(o2v(V.d));
printf("  d\n");

printpol(o2v(V.h));
printf("  Uh\n");
//exit(1);

s1=omul(c1,e1);
printpol(o2v(s1));
printf(" ==s1\n");
s2=omul(c1,e2);
printpol(o2v(s2));
printf(" ==s2\n");
s3=c2;
printpol(o2v(s3));
printf(" ==s3\n");
//exit(1);
int count=0;
OP v;

u=odiv(omul(uu1,uu2),omul(d,d));
printpol(o2v(u));
printf(" ==u\n");

reduct:
count++;
v=omod(oadd(oadd(omul(omul(s1,uu1),vv2),omul(omul(s2,uu2),vv1)),omul(s3,oadd(omul(vv1,vv2),ff))),u);
printpol(o2v(v));
printf(" ==v\n");
//exit(1);

OP ud,vd;
ud=odiv(osub(ff,omul(v,v)),u);
vd=omod(minus(v),ud);
if(odeg(ud)>2){
  if(count>100){
    printf("over 100\n");
    D3.u=ud;
    D3.v=vd;
    return D3;
  //exit(1);
  }
  u=ud;
  v=vd;
  printf("==================\n");
  goto reduct;
}
ud=monique(ud);
printpol(o2v(ud));
printf(" @@ud\n");
printpol(o2v(vd));
printf(" @@udv\n");
//printpoln(o2v(oadd(vv1,vv2)));
D3.u=ud;
D3.v=vd;

return D3;
}

Div dobule(Div D, OP ff){
vec vx=o2v(ff);
unsigned long long d0,uv010,uu1,uu0,uuv11,d1,d2,d3,d4,d5,f2,e0,e1,f1,f3,mm0,mm1,mm3,mm4,s1,s2,l1,l0,l2,l3,u31,u30,v31,v30,uue0,uue1,iv,uuv01,uuv10,uuv00,u1,u0,v1,v0;
OP u,v;
Div D2;
f3=vx.x[3];
f2=vx.x[2];
f1=vx.x[1];

u=D.u;
v=D.v;
memset(vx.x,0,sizeof(vx.x));
vx=o2v(u);
u1=vx.x[1];
u0=vx.x[0];

memset(vx.x,0,sizeof(vx.x));
vx=o2v(v);
v1=vx.x[1];
v0=vx.x[0];

uu1=(u1*u1)%P; uu0=(u1*u0)%P; uuv01=(u0*v1)%P; uuv10=(u1*v0)%P; uuv11=(u1*v1)%P;
uuv00=(u0*v0)%P;

  d0=(6*v1*uu1-uuv10)%P; d1= (-4*uuv11+4*v0)%P; d2= (2*v1)%P;
  d3=(6*v1*uu0-6*uuv00)%P; d4= (-4*uuv01)%P; d5=(2*v0)%P;
  e0=(5*(-u1*uu1+2*uu0)-3*f3*u1+2*f2)%P;
  e1=(5*(-u0+uu0+u0*u0)-3*f3*u0+f1)%P;
  mm0=(d3-d5*(uu1-u0))%P; mm1=(d4-d5*(-u1))%P;
  mm3=(d0-d2*(uu1-u0))%P; mm4=(d1-d2*(-u1))%P;
  s1=(e1-d5*v1)%P; s2=(e0-d2*v1)%P;
  iv=inv((mm0*mm4-mm1*mm3),P);
  l3= (iv*(mm4*s1-mm1*s2))%P; l2=(iv*(mm0*s2-mm3*s1))%P;
  l1=(v1+u1*l2-(uu1-u0)*l3)%P; l0=(v0+u0*l2-(uu0)*l3)%P;
  u31=(2*l3*l2-2*u1-1)%P;
  u30=(2*l3*l1+l2*l2-2*u0-uu1-2*u31*u1)%P;
  uue1=(u31*u31)%P; uue0=(u31*u30)%P;
  v31=((uu1-u30)*l3-u31*l2+l1)%P;
  v30=(uue0*l3-u30*l2+l0)%P;

memset(vx.x,0,sizeof(vx.x));
vx.x[1]=u31;
vx.x[0]=u30;
vx.x[2]=1;
u=v2o(vx);
memset(vx.x,0,sizeof(vx.x));
vx.x[1]=v31;
vx.x[0]=v30;
vx.x[2]=1;
v=v2o(vx);

D2.u=u;
D2.v=v;

return D2;
}

Div g2add(OP ff, OP uu1, OP uu2, OP vv1, OP vv2)
{
  OP ll, u;
  OP v, s, l, k, v3, u3;
  Div X;

  ll = (omul(vv2, vv2));
  printpol(o2v(ll));
  printf("\n");
  ll = (osub(ff, (ll)));
  printpol(o2v(ll));
  printf("\n");
  k = (odiv(ll, uu2));
  printpol(o2v(k));
  printf(" ('A`)\n");
  //exit(1);
  ll=qinv(uu2,uu1);
  s=omod(omul(osub(vv1,vv2),ll),uu1);
  printpol(o2v(ll));
  printf("========inv\n");
l=omul(s,uu2);
u=odiv(osub(k,omul(s,oadd(l,omul(s,vv2)))),uu1);
printpol(o2v(u));
printf(" ===u\n");
u3=monique(u);
printpol(o2v(u3));
printf(" ===u3\n");
v3=omod(minus(oadd(l,vv2)),u3);
printpol(o2v(v3));
printf(" ===v3\n");
exit(1);


  X.u = u3;
  X.v = v3;

  return X;
}


PO tr1e(int f4, int f3, int f2, int f1, int f0, int p)
{
  int x,y,f,g;
  PO aa = {0};

  while (1)
  {
    x = rand() % p;
    y = rand() % p;
    f=((x * x * x * x * x)%p + (f4 * x * x * x * x)%p + (f3 * x * x * x)%p + (f2 * x * x)%p + f1 * x + f0)%p;
    g=(y*y)%p;
    if(f == g){
      aa.x=x;
      aa.y=y;
      return aa;
    }
  }
  //  return -1;
}


Div gendiv(OP f)
{
  PO a = {0}, b = {0}, e = {0};
  OP d1 = {0}, d2 = {0}, c = {0}, d = {0};
  int x, y, i, j, k, count;
  Div D = {0};
  vec v1 = {0}, v2 = {0}, z1 = {0}, z2 = {0}, ff = {0};
  count = 0;

  ff = o2v(f);
    v1.x[1] = 1;
    v2.x[1] = 1;
if(deg(v1)==0 || deg(v2)==0){
  printf("ee!?\n");
  exit(1);
}
  do
  {  
    a = tr1e(ff.x[4], ff.x[3], ff.x[2], ff.x[1], ff.x[0], P); // cofficient of function
    e = tr1e(ff.x[4], ff.x[3], ff.x[2], ff.x[1], ff.x[0], P); // cofficient of function
    b = tr1e(ff.x[4], ff.x[3], ff.x[2], ff.x[1], ff.x[0], P); // cofficient of function
  
    v1.x[0] = a.x;
    printpol(v1);
    printf("ppppppppppppp\n");
    c = v2o(v1);

    v2.x[0] = b.x;
    d = v2o(v2);

    z1.x[1] = rand()%P;
    z1.x[0] = rand()%P;


    d2 = v2o(z1);
    d1 = omul(c, d);

  } while (LT((omod(osub(omul(d2, d2), (f)), d1))).a != 0);
printf("debug mode\n");
  printpol(o2v(d1));
  printf(" ==u\n");
  printpol(o2v(d2));
  printf(" ==v\n");
  // exit(1);

  D.u = d1;
  D.v = d2;
if(chkdiv(D,f)==-1){
  printf("so buggy!\n");
  exit(1);
}
  return D;
}


EX manford(OP a, OP b)
{
  EX V;

  V = xgcd(a, b);
  printpol(o2v(V.u));
  printf(" =====u3\n");
  printpol(o2v((V.v)));
  printf(" =====v3\n");
  printpol(o2v(V.d));
  printf(" =====d3\n");
  printpol(o2v((V.h)));
  printf(" =====h3\n");

  return V;
}




int main()
{
  unsigned int i, count = 0;
  unsigned long long  aaa[O] = {0};
  

  unsigned long long  f[K + 1] = {1, 7, 6, 2, 8, 2};
/*
  unsigned long long  u2[K + 1] = {0, 0, 0, 1, 21, 16};
  unsigned long long  u1[K + 1] = {0, 0, 0, 1, 19, 20};
  unsigned long long  v2[K + 1] = {0, 0, 0, 0, 21, 21};
  unsigned long long  v1[K + 1] = {0, 0, 0, 0, 12, 8};
*/
// unsigned long long  f[K+1]={1 ,0, 2,  30,  5,  1};
/*
  //unsigned long long  f[K+1]= {1, 1597 , 1041 ,5503 , 6101 , 1887 };
//f1 = x + 28555025517563816 and f2 = x + 74658844563359755 ;
unsigned long long  u2[K+1]={0,0,0,1,1571353025997967 , 12198441063534328};
unsigned long long  v2[K+1]={0,0,0,0,32227723250469108 , 68133247565452990};
unsigned long long  u1[K+1]={0,0,0,1, 70887725815800572 , 94321182398888258};
unsigned long long  v1[K+1]={0,0,0,0, 42016761890161508 , 3182371156137467 };
*/
unsigned long long u2[K+1]={0,0,0,0,1,22};
unsigned long long v2[K+1]={0,0,0,0,13,8};
unsigned long long u1[K+1]= {0,0,0,1,22,24};
unsigned long long v1[K+1]= {0,0,0,0,5,21};

//  unsigned long long  f[K + 1] = {1, 0, 3, 7, 1, 2};
/*
  unsigned long long  u2[K + 1] = {0, 0, 0, 1, 7, 10};
  unsigned long long  u1[K + 1] = {0, 0, 0, 1, 0, 10};
  unsigned long long  v2[K + 1] = {0, 0, 0, 0, 1, 9};
  unsigned long long  v1[K + 1] = {0, 0, 0, 0, 7, 9};
*/
  OP ff, k, uu1, uu2, vv1, vv2, s, l, u3, v3, u, ll, t, m, o, d, c;
  unsigned long long  tst1[K + 1] = {0, 0, 0, 0, 8, 7};
  unsigned long long  tst2[K + 1] = {0, 0, 0, 0, 0, 10};
  unsigned tmp[2][3]={0};
  EX V;
  oterm a;
  OP b={0};
//unsigned long long  q1 = to_unsigned long long ("1208925819614629174708801");
int a1 = 1331;
//J1 =to_unsigned long long ("1461501637326815988079848163961117521046955445901");
//e y2 = x5+a, a ∈ Fp

//unsigned long long  q2 = to_unsigned long long ("1208925819614629174709941");
int a2 = 2;
//J2 = to_unsigned long long ("1461501637331762771847359428275278989652932675771");

vec vx={0},xv={0};
Div G0,G1,X;
  ff = setpol(f, K + 1);

  uu1 = setpol(u1, K + 1);
  uu2 = setpol(u2, K + 1);
  vv1 = setpol(v1, K + 1);
  vv2 = setpol(v2, K + 1);
  o=setpol(tst1,K+1);
  m=setpol(tst2,K+1);
G1.u=uu1;
G1.v=vv1;
X.u=uu2;
X.v=vv2;
if(chkdiv(G1,ff)==-1 || chkdiv(X,ff)==-1){
    printf("erro!\n");
exit(1);
}
count=0;
int xount=0;

PO xx;
unsigned long long rr=0;

G0=cadd(ff,uu1,uu2,vv1,vv2);
if(chkdiv(G0,ff)==-1){
    printf("bug\n");
    exit(1);
}
/*
for(rr=0;rr<1000;rr++){
  printf("%llu %llu\n",root(rr,P),rr);
}
exit(1);
vx=o2v(ff);
xx = tr1e(vx.x[4], vx.x[3], vx.x[2], vx.x[1], vx.x[0], P); // cofficient of function
printf("%llu, %llu\n",xx.x, xx.y);
exit(1);
*/
/*
for(i=0;i<1000;i++){
  x=i;
  
}
*/
G1=gendiv(ff);
X=gendiv(ff);
uu1=G1.u;
vv1=G1.v;
srand(clock());
while(1){
G1=gendiv(ff);
X=gendiv(ff);
G0=cadd(ff,G1.u,X.u,G1.v,X.v);
//if(oequ(G1.u,uu1)==0)
//exit(1);
//G0=dobule(X,ff);
/*
while((LT(xgcd(X.v,X.u).d).a!=1)){
  X=gendiv(ff);
  //G0=cdbl(X,ff);
}
*/
if(chkdiv(G0,ff)==-1)
{
printf("baka\n");
//    exit(1);
  V=xgcd(X.u,G1.u);
  if(LT(V.d).n>0){
  printf("gcd!\n");
  }else if(chkdiv(G1,ff)==-1 || chkdiv(X,ff)==-1){
    printf("gen div err\n");
    break;
  }else{
    printpoln(o2v(G1.u));
    printpoln(o2v(G1.v));
    printpoln(o2v(X.u));
    printpoln(o2v(X.v));
  printf("why?\n");
  count++;
  exit(1);
  }
}else{
  printf("ウホッ！いい因子。\n");
  xount++;
  //exit(1);
}
if(xount>100)
break;
printf("%u xount=%u\n",count,xount);
}
printf("%u %u\n",count,xount);
exit(1);

/*
G0=cdbl(X,ff);
if(chkdiv(G0,ff)==-1)
{
  printf("naze?\n");
}else{
  printf("イイっ！この因子すげえいいっ！\n");
}
*/
//exit(1);
/*
G0=g2add(ff,uu1,uu2,vv1,vv2);
if(chkdiv(G0,ff)==-1){
  printf("nani?\n");
}else{
printf("だろうな\n");
}
exit(1);
*/
//b=odiv(o,m);
//exit(1);
/*
X=g2add(ff,uu1,uu2,vv1,vv2);
printpol(o2v(X.u));
printf("  Xu\n");
printpol(o2v(X.v));
printf("  Xv\n");
//printf("%llu\n",chkdiv(X,ff));
exit(1);
*/
srand(clock());

/*
EX F;
Div U={0},U0={0};
while(1)
{
U=gendiv(ff);
U0=gendiv(ff);
printf("1111111111111111111111111111\n");
//F=xgcd(U.u,U0.u);
if(chkdiv(U,ff)==-1 || chkdiv(U0,ff)==-1)
{
  printf("ee?!\n");
  exit(1);
}
X=cadd(ff,U.u,U0.u,U.v,U0.v);
if(chkdiv(X,ff)==-1)
{
  printpoln(o2v(U.u));
  printpoln(o2v(U0.u));
  printpoln(o2v(U.v));
  printpoln(o2v(U0.v));
  printf("baka^^\n");
  exit(1);
}else{
  printf("happy!\n");
  //exit(1);
}
}
printf("isdiv=%llu\n",chkdiv(X,ff));
printpol(o2v(F.u));
printf("  Fuwwwwwwww\n");
printpol(o2v(F.v));
printf("  Fvwwwwwwww\n");
printpol(o2v(F.d));
printf("  Fdwwwwwwww\n");
printpol(o2v(F.h));
printf("  Fhwwwwwwww\n");
printpol(o2v(X.u));
printf("  Xuwwwwwwww\n");
printpol(o2v(X.v));
printf("  Xvwwwwwwww\n");
printpol(o2v(U.u));
printf("  Uu wwwwwwww\n");
printpol(o2v(U0.u));
printf("  U0u vwwwwwwww\n");
printpol(o2v(U.v));
printf("  Uv wwwwwwww\n");
printpol(o2v(U0.v));
printf("  U0v vwwwwwwww\n");
*/
//exit(1);
/*
//printf("isdiv=%llu\n",chkdiv(U,ff));
//printf("isdiv=%llu\n",chkdiv(U0,ff));
X=g2add(ff,U.u,U0.u,U.v,U0.v);
printpol(o2v(U.u));
printf("  isUwwwwwwww\n");
printpol(o2v(U0.u));
printf("  isU0wwwwwwww\n");
  printf("isdiv=%llu\n",chkdiv(X,ff));
  exit(1);
  */
OP oi,rem;
Div D1,D2;
/*
while(1){
  D1=gendiv(ff);
  D2=gendiv(ff);
  printpol(o2v(D1.u));
  printf(" ===D1.u\n");
  printpol(o2v(D1.v));
  printf(" ===D1.v\n");
  printpol(o2v(D2.u));
  printf(" ===D2.u\n");
  printpol(o2v(D2.v));
  printf(" ===D2.v\n");
  oi=qinv(D2.u,D1.u);
  printpol(o2v(oi));
  printf(" ===Doi\n");
  rem=omod(omul(oi,D2.u),D1.u);
  printpol(o2v(rem));
  printf(" ==D rem\n");
  exit(1);
  printf("isdiv=%llu D\n",chkdiv(D1,ff));
  printf("isdiv=%llu D\n",chkdiv(D2,ff));
  U = g2add(ff, D1.u, D2.u, D1.v, D2.v);
  printf("isdiv=%llu D\n",chkdiv(U,ff));
  //exit(1);
  
  if(chkdiv(U,ff)==-1){
    printf("buggy\n");
  exit(1);
  }
}
  //U = g2add(ff, uu1, uu2, vv1, vv2);
  //printf("isdiv=%llu\n",chkdiv(U,ff));
//    printf("isdiv=%llu\n",chkdiv(D2,ff));
//  exit(1);
  //V=xgcd(uu1,uu2,2);
  printpol(o2v(V.u));
  printf(" =====u3\n");
  printpol(o2v(V.v));
  printf(" =====v3\n");
*/

D1.u=uu1;
D1.v=vv1;
D2.u=uu2;
D2.v=vv2;

printf("isdiv=%llu\n",chkdiv(D1,ff));
printf("isdiv=%llu\n",chkdiv(D2,ff));
//exit(1);

D1=gendiv(ff);
D2=gendiv(ff);

printf("isdiv=%llu\n",chkdiv(D1,ff));
printf("isdiv=%llu\n",chkdiv(D2,ff));
//exit(1);
//U=g2add(ff,D1.u,D2.u,D1.v,D2.v);

  o=oadd(vv1,vv2);
  printpol(o2v(o));
  printf("\n");

  Div W;
  // below undercondtruction
  k = odiv(oadd(ff, minus(omul(D1.v, D1.v))), D1.u);
  s = omod(odiv(k, scr(2, D1.v)), D1.u);
  l = omul(s, uu1);
  u3 = omod(osub(omul(s, s), (osub(scr(2, omul(D1.v, s)), (k)))), D1.u);
  v3 = omod(minus(oadd(l, D1.v)), u3);
  W.u=u3;
  W.v=v3;
  printpol(o2v(u3));
  printf("======du3\n");
  printpol(o2v(v3));
  printf("======dv3\n");
  printf("W's isdiv=%llu\n",chkdiv(W,ff));
//  exit(1);


  return 0;
}
