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
#define P 11 //37

// using namespa
//  sagemath上での原始多項式
unsigned short pp[4][4] = {{0, 0, 9, 2}, {0, 0, 11, 2}, {0, 0, 16, 3}, {0, 0, 15, 2}};
// {0,0,9,2}, {1,0,11,2}, {1,0,16,3}, {1,0,15,2};
// GF(11^3,13^3,17^3,19^3)
// unsigned short ff[2][7]={{1,0,0,0,0,2,0,2},{0,0,1,0,0,0,1,2}}; //GF(3^7,5^5)

unsigned short gf[O] = {0}, fg[O] = {0};
// int N =0,M=0;
unsigned short c[K + 1] = {0};

uint32_t xor64(void) {
  static uint64_t x = 88172645463325252ULL;
  x = x ^ (x << 13); x = x ^ (x >> 7);
  return x = x ^ (x << 17);
}


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
      f.t[j].n = i;
      f.t[j++].a = a.x[i];
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
OP setpol(unsigned short f[], int n)
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
//    printf("\n");

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
unsigned short
xtrace(OP f, unsigned short x)
{
  int i, d;
  unsigned short u = 0, v = 1;

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
  unsigned short i, j, count = 0;

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
  printf("unsigned short fg[%d]={", O);
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
  // exit(1);

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
  unsigned short n;

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
      d.x[i] = (P + (a.x[i] - b.x[i]));
    }

    printf("%d - %d = %d\n", a.x[i], b.x[i], d.x[i]);
    // if(d.x[i]<0){
    //   printf("%d\n",d.x[i]);
    //   d.x[i]+=P;
    // }
  }
  // exit(1);
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

// aに何をかけたらbになるか
unsigned short
equ(unsigned short a, unsigned short b)
{
  int i;

  for (i = 0; i < P; i++)
  {
    if ((a * i) % P == b)
      break;
  }
  return i;
}

void mkmf()
{
  int i, j, k, count = 0;
  OP f = {0}, g = {0}, h = {0}, w = {0}, s = {0}, u = {0};
  vec b = {0}, a = {0}, d = {0}, t = {0}, v = {0};
  oterm o;
  unsigned short ccp[4] = {0};

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
      // w=omod(w,u);
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

  printf("unsigned short gf[%d]={", O);
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
    s.a = (tt.a * inv(t.a, P)) % P;
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
  if (LT(g).a == 0 || LT(f).a == 0)
  {
    printf("g & f is 0!\n");
    return o;
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

  // assert(("double baka\n", b.a > 0 && b.n > 0));
  while (LT(f).n > -1 && LT(g).n > -1)
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
    if (odeg((f)) == 0 || odeg((g)) == 0)
    {
      //      printf("blake1\n");
      break;
    }

    if (c.n == 0 || b.n == 0)
      break;
    if (LT(f).a == 4 && deg(o2v(f)) == 0)
      exit(1);
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

  if (LT(f).a == 0 || LT(g).a == 0)
  {
    printf("baka^\n");
    // return f;
    exit(1);
  }
  // if (LT (g).n == 0 && LT (g).a > 1)
  //   return g; //coeff (f);

  k = 0; // odeg(f) - odeg(g);
  b = LT(g);
  if (b.a == 1 && b.n == 0)
    return f;
  if (b.a == 0 && b.n == 0)
  {
    printf("baka in odiv\n");
    exit(1);
  }
  if (odeg((f)) < odeg((g)))
  {
    return o;
  }
  OP null = {0};
  i = 0;
  k = 0;
  while (LT(f).a > -1 || LT(g).a > -1)
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

  int tt_terms = terms(tt);
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

// invert of polynomial
OP inv3(OP a, OP n)
{
  OP d = {0}, x = {0}, s = {0}, q = {0}, r = {0}, t = {0}, u = {0}, v = {0}, w = {0}, tt = {0}, gcd = {0}, tmp = {0};
  oterm b = {0};
  vec vv = {0}, xx = {0}, aa = {0}, bb = {0}, cc = {0};

  if (odeg((a)) > odeg((n)))
  {
    tmp = a;
    a = n;
    n = tmp;
    printf("baka_i\n");
    // exit (1);
  }
  if (LT(a).a == 0)
  {
    printf(" a ga 0\n");
    exit(1);
  }

  tt = n;

  d = n;
  x.t[0].a = 0;
  x.t[0].n = 0;
  s.t[0].a = 1;
  s.t[0].n = 0;
  while (odeg((a)) > 1)
  {
    if (odeg((a)) > 0)
      r = omod(d, a);
    if (LT(a).a == 0)
      break;
    if (LT(a).a > 0)
      q = odiv(d, a);

    d = a;
    a = r;
    t = osub(x, omul(q, s));
    ////printpol (o2v (a));
    // printf ("\nin roop a==================%d\n", odeg ((a)));
    // printf ("\n");

    x = s;
    s = t;
  }
  // exit(1);
  //  if(LT(a).a>0){
  d = a;
  a = r;
  ////printpol (o2v (a));
  // printf ("\nin roop a|==================%d\n", odeg ((a)));
  // printf ("\n");

  x = s;
  s = t;

  ////printpol (o2v (d));
  // printf ("\nout1================\n");
  gcd = d; // $\gcd(a, n)$
  printpol(o2v(gcd));
  printf(" =========gcd\n");
  // exit(1);
  // printf ("\n");
  ////printpol (o2v (n));
  // printf ("\n");
  // printf ("out2===============\n");

  printf("before odiv\n");
  // w=tt;

  b = LT(w);
  ////printpol (o2v (w));
  // printf ("\nw=======%d %d\n", b.a, b.n);
  // w=tt;
  aa = o2v(x);
  bb = o2v(n);
  v = oadd(x, n);

  ////printpol (o2v (v));
  // printf ("\n");
  /*
     if (LT (v).a == 0)
     {
     printf ("v=============0\n");
     }
     printf ("d==============\n");
   */
  //  } //end of a>0
  w = tt;
  ////printpol (o2v (v));
  // printf ("\n");
  // printf ("ss==============\n");
  //        exit(1);
  //  if(odeg((w))>0)
  if (LT(v).n > 0 && LT(w).n > 0)
  {
    u = omod(v, w);
  }
  else
  {
    printpol(o2v(v));
    printf(" v===========\n");
    printpol(o2v(x));
    printf(" x==0?\n");
    printpol(o2v(n));
    printf(" n==0?\n");

    exit(1);
  }
  // caution !!
  if (LT(u).a > 0 && LT(d).a > 0)
  {
    u = odiv(u, d);
  }

  if (LT(u).a == 0 || LT(d).a == 0)
  {
    printf("inv div u or d==0\n");
    // exit(1);
  }
  // u=coeff(u,d.t[0].a);
  ////printpol (o2v (u));
  // printf ("\nu==================\n");
  if (LT(u).a == 0)
  {
    printf("no return at u==0\n");
    exit(1);
  }
  // return ((x + n) % (n / d));
  return u;
}

OP scr(unsigned short d, OP f)
{
  int i, n;
  vec v = {0};

  n = deg(o2v(f));
  v = o2v(f);
  for (i = 0; i < n + 1; i++)
    v.x[i] = (v.x[i] * d) % P;
  f = v2o(v);

  return f;
}

OP monic(OP f)
{
  int e1;
  e1 = inv(LT(f).a, P);
  printf("e=%d\n", e1);
  f = scr(e1, f);

  return f;
}

OP cdiv(int a, OP f)
{
  vec v;
  int i, l;

  v = o2v(f);
  l = odeg(f);
  for (i = 0; i < l; i++)
    v.x[i] = (a * v.x[i]) % P;

  f = v2o(v);

  return f;
}

//拡張ユークリッドアルゴリズム
EX xgcd(OP f, OP g)
{
  OP h[10] = {0}, ww[10] = {0}, *v, *u;
  oterm a, b;
  int i = 0, j, flg = 0, k;
  EX e = {0}, ee = {0};

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
  if (LT(f).a == 0 || LT(g).a == 0)
  {
    printf("f or g ==0\n");
    exit(1);
  }

  k = odeg(g);
  if (k == 0 && LT(g).a > 0)
  {
    printf("use cdiv\n");
    exit(1);
  }
  i = 1;
  while (LT(g).a != 0)
  //  for (i = 1; i < 1+ k+1; i++)
  {
    printpol(o2v(f));
    printf(" fffffffffffffffff\n");

    // if(LT(g).a>-0){
    printpol(o2v(g));
    printf(" ggggggggggggggggg\n");
    // if (LT(g).a > 0)
    h[i] = omod(f, g);
    printpol(o2v(h[i]));
    printf(" %d hhhhhhhhhhhhhh\n", i);

    // if (LT(g).a > 0)
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

    // if (LT(g).a == 0)
    //   break;
    i++;
  }
  // f=g;
  // g=h[i];
  if (LT(g).a > 0)
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
  // exit(1);

  e.d = f;
  e.u = u[i];
  e.v = v[i];
  e.h = ww[i];

  free(v);
  free(u);

  printf("end of fnc\n");
  // exit(1);
  return e;
  //  wait ();
}

OP qinv(OP uu1, OP uu2)
{
  EX tt, V;
  OP v;
  int e;

  // return cdiv(LT(g).a,g);
  tt = xgcd(uu1, uu2);
  printpol(o2v(tt.v));
  printf("\n");
  printpol(o2v(tt.u));
  printf("\n");
  printpol(o2v(tt.d));
  printf("\n");
  printpol(o2v(tt.h));
  printf("\n");
  // exit(1);
  v = omul(tt.d, tt.v);
  printpol(o2v(v));
  printf(" in qinv's v\n");

  // exit(1);

  return v;
}

int chkdiv(Div d, OP f)
{
  OP t;

  t=(omod(oadd(omul(d.v, d.v), minus(f)), d.u));
  printpol(o2v(t));
  printf(" 00000000000\n");
  if (LT(t).a==0)
    return 1;

  return -1;
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
  ll = odiv(uu1, uu2);
  printpol(o2v(ll));
  printf("========div\n");
  ll = omod(uu2, uu1);
  printpol(o2v(ll));
  printf("=======mod\n");
  // exit(1);
  EX tt = {0};

  // tt=muri(uu2,uu1);
  v = qinv(uu1, uu2);
  printpol(o2v(v));
  // omod(omul(t,uu2),uu1);
  printf(" ===inv\n");
  //   exit(1);

  // v=oadd(vv1,vv2);
  // printpol(o2v(v));
  // printf(" v\n");
  // exit(1);
  // tt=xgcd(uu1,uu2);
  if(LT(vv1).a!=0 && LT(vv2).a!=0){
  tt = xgcd(oadd(vv1, vv2), v);
  printpol(o2v(tt.d));
  printf(" d@\n");
  printpol(o2v(tt.u));
  printf(" u@\n");
  printpol(o2v(tt.v));
  printf(" v@\n");
  printpol(o2v(tt.h));
  printf(" h@\n");
  // exit(1);
  }else{
    printf("v==0\n");
    exit(1);
  }
  ll = osub(vv1, (vv2));
  printpol(o2v(ll));
  printf("\n");
  s = omod(omul(ll, v), uu1);
  printpol(o2v(s));
  printf("\n");
  // exit(1);

  l = omul(s, uu2);
  printpol(o2v(l));
  printf("\n");
  u = odiv(oadd(k, minus(omul(s, oadd(l, omul(s, vv2))))), uu1);
  printpol(o2v(u));
  printf("\n");
  u3 = coeff(u);
  printpol(o2v(u3));
  printf(" =======u3\n");
  v3 = omod(minus(oadd(l, vv2)), u3);
  printpol(o2v(v3));
  printf(" =========v3\n");
   exit(1);

  X.u = u3;
  X.v = v3;

if(chkdiv(X,ff)==-1){
  printf("failed!\n");
exit(1);
}else{
  printf("g00d\n");
  exit(1);
}

  return X;
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

int bit(unsigned b, int i)
{
  int k = 1;

  if ((b & (1 << i)) > 0)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

// jj=aa^bb mod oo
int Exp(int aa, int bb, int oo)
{
  int ii, jj, kk[8192];
  int j, c[8192], count = 0, i;
  ii = oo;
  j = 0;
  jj = 0;
  //  kk[4096]; //prime is 4096 bit table
  //  c[8192]  //mod is 8192 bit table
  count = 0;

  for (i = 0; i < 8192; i++)
  {
    kk[i] = 0;
  }
  while (ii > 0)
  {
    ii = (ii >> 1);
    j = j + 1;
  }

  kk[0] = aa;

  //  cout << j << "\n";

  // ex.1000=2**3+2**5+2**6+2**7+2**8+2**9 makes a array c=[3,5,6,7,8,9]
  for (i = 0; i < j + 1; i++)
  {
    if (bit(bb, i) != 0)
    { // testbit(bb,i)
      c[count] = i;
      count = count + 1;
    }
  }
  //    cout << bb << endl;
  //    cout << count << "\n";
  // exit(1);
  for (i = 1; i < c[count - 1] + 1; i++)
  {
    kk[i] = kk[i - 1] * kk[i - 1] % oo;
  }

  jj = 1;
  for (i = 0; i < count; i++)
  {
    jj = kk[c[i]] * jj % oo;
    if (jj == 0)
    {
      //	print i,"\n"
    }
  }

  return jj;
}

long pow_mod(long x, long n, long p)
{
  if (n == 0)
    return 1;
  if (n & 1)
    return (pow_mod(x, n - 1, p) * x) % p;
  x = pow_mod(x, n / 2, p);
  return (x * x) % p;
}

/* Takes as input an odd prime p and n < p and returns r
 * such that r * r = n [mod p]. */
long tonelli_shanks(long n, long p)
{
  long s = 0;
  long q = p - 1;
  while ((q & 1) == 0)
  {
    q /= 2;
    ++s;
  }
  if (s == 1)
  {
    long r = pow_mod(n, (p + 1) / 4, p);
    if ((r * r) % p == n)
      return r;
    return 0;
  }
  // Find the first quadratic non-residue z by brute-force search
  long z = 1;
  while (pow_mod(++z, (p - 1) / 2, p) != p - 1)
    ;
  long c = pow_mod(z, q, p);
  long r = pow_mod(n, (q + 1) / 2, p);
  long t = pow_mod(n, q, p);
  long m = s;
  while (t != 1)
  {
    long tt = t;
    long i = 0;
    while (tt != 1)
    {
      tt = (tt * tt) % p;
      ++i;
      if (i == m)
        return 0;
    }
    long b = pow_mod(c, pow_mod(2, m - i - 1, p - 1), p);
    long b2 = (b * b) % p;
    r = (r * b) % p;
    t = (t * b2) % p;
    c = b2;
    m = i;
  }
  if ((r * r) % p == n)
    return r;
  return 0;
}

unsigned root(unsigned a, unsigned p)
{
  int c, b;

  printf("p mod = %d , %d\n", p % 4, p % 8);
  if (p % 4 == 3 || p % 8 == 5)
  {
    if (p % 4 == 3)
    {
      b = (p + 1) / 4;
      c = Exp(a, b, p);
      if (c * c % p != a)
      {
        printf("baka1\n");
        return -1;
        //#exit()
      }
      if (c * c % p == a)
      {
        printf("good\n");
        printf("c=%d\n", c);
        return c;
      }
    }
    if (p % 8 == 5)
    {
      c = Exp(a, (p + 3) / 8, p);
      if (c * c % p != a)
        printf("baka\n");
        return -1;
    }
    if (c * c % p == a)
    {
      printf("good\n");
      printf("%d\n", c);
      return c;
    }
  }
  if (p % 8 == 5)
  {
    c = 2 * a * Exp(4 * a, (p - 5) / 8, p);
    if (c * c % p != a)
    {
      printf("dangerous\n");
      return -1;
    }
    if (c * c % p == a)
    {
      printf("good\n");
      printf("%d\n", c);
      return c;
    }
    if (p % 8 == 1){
      c = tonelli_shanks(a, p);
      if(c*c%p != a){
        printf("fail!\n");
        return -1;
      }else{
      return c;
      }
    }
    return 0;
  }

  return -1;
}

PO tr1e(int f4, int f3, int f2, int f1, int f0, int p)
{
  int x,y;
  PO aa = {0};

  while (1)
  {
    x = rand() % p;
    y = rand() % p;
    if((y*y)%p == (x * x * x * x * x + f4 * x * x * x * x + f3 * x * x * x + f2 * x * x + f1 * x + f0) % p){
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

  do
  {
    a = tr1e(ff.x[4], ff.x[3], ff.x[2], ff.x[1], ff.x[0], P); // cofficient of function
    e = tr1e(ff.x[4], ff.x[3], ff.x[2], ff.x[1], ff.x[0], P); // cofficient of function
    b = tr1e(ff.x[4], ff.x[3], ff.x[2], ff.x[1], ff.x[0], P); // cofficient of function

    if(a.x!=-1){ 
    v1.x[0] = a.x;
    v1.x[1] = 1;
    c = v2o(v1);
    }
    if(b.x!= -1){
    v2.x[0] = b.x;
    v2.x[1] = 1;
    d = v2o(v2);
    }
    if(e.x!= -1){
    z1.x[1] = e.x;
    z1.x[0] = e.y;
    d2 = v2o(z1);
    }
    d1 = omul(c, d);
    
  } while (LT((omod(oadd(omul(d2, d2), minus(f)), d1))).a != 0);

  printpol(o2v(d1));
  printf(" ==u\n");
  printpol(o2v(d2));
  printf(" ==v\n");
  //  exit(1);

  D.u = d1;
  D.v = d2;

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
  unsigned int i, count = 0, c2 = 0;
  unsigned short aaa[O] = {0};
  /*
    unsigned short f[K + 1] = {1, 7, 6, 2, 8, 2};
    unsigned short u1[K + 1] = {0, 0, 0, 1, 21, 16};
    unsigned short u2[K + 1] = {0, 0, 0, 1, 19, 20};
    unsigned short v1[K + 1] = {0, 0, 0, 0, 21, 21};
    unsigned short v2[K + 1] = {0, 0, 0, 0, 12, 8};
  */

  unsigned short f[K + 1] = {1, 0, 3, 7, 1, 2};
  unsigned short u1[K + 1] = {0, 0, 0, 1, 7, 10};
  unsigned short u2[K + 1] = {0, 0, 0, 1, 0, 10};
  unsigned short v1[K + 1] = {0, 0, 0, 0, 1, 9};
  unsigned short v2[K + 1] = {0, 0, 0, 0, 7, 9};

  OP ff, k, uu1, uu2, vv1, vv2, s, l, u3, v3, u, ll, t, m, o, d, c;
  unsigned short tst1[K + 1] = {0, 0, 0, 0, 1, 29};
  unsigned short tst2[K + 1] = {0, 0, 0, 0, 2, 29};
  unsigned tmp[2][3] = {0};
  Div V;
  oterm a;
  OP b = {0};
  // int q1 = to_int("1208925819614629174708801");
  int a1 = 1331;
  // J1 =to_int("1461501637326815988079848163961117521046955445901");
  // e y2 = x5+a, a ∈ Fp

  // int q2 = to_int("1208925819614629174709941");
  int a2 = 2;
  // J2 = to_int("1461501637331762771847359428275278989652932675771");
  int j, t1[2][3] = {0}, c1[2] = {0}, cc[2] = {0};
  vec vx = {0};
  Div Z;
  
  
  srand(clock());

  ff = setpol(f, K + 1);
  uu1 = setpol(u1, K + 1);
  uu2 = setpol(u2, K + 1);
  vv1 = setpol(v1, K + 1);
  vv2 = setpol(v2, K + 1);
  o = setpol(tst1, K + 1);
  m = setpol(tst2, K + 1);

Z.u=uu1;
Z.v=vv1;

/*
  printpol(o2v(vv1));
  printf("\n");
  printpol(o2v(vv2));
  printf("\n");
  // exit(1);

  vx = diviser(vv1, vv2);
  printpol(vx);
  printf("\n");
  V.d = v2o(vx);
  //exit(1);

  //  exit(1);
  V = manford(uu1, uu2);
  d = oadd(vv1, vv2);
  printpol(o2v(d));
  printf(" v1+v2\n");
  // exit(1);

  vx = diviser(V.d, d);
  printpol(vx);
  printf("\n");
  V.d = v2o(vx);
  //exit(1);

  V = manford(V.d, d);
  // exit(1);

  unsigned short f1[K + 1] = {1, 0, 2, 30, 5, 1};
  unsigned short d1[K + 1] = {0, 0, 0, 1, 3, 0};
  unsigned short l1[K + 1] = {0, 0, 0, 1, 9, 20};
  unsigned short d2[K + 1] = {0, 0, 0, 0, 11, 14};
  unsigned short l2[K + 1] = {0, 0, 0, 0, 15, 13};

  OP ff1, dd1, dd2, v0, vv;


  ff1 = setpol(f1, K + 1);
  dd1 = setpol(d1, K + 1);
  dd2 = setpol(d1, K + 1);
  v0 = setpol(l1, K + 1);
  vv = setpol(l2, K + 1);

  D1 = gendiv(ff1);
  D2 = gendiv(ff1);
  //exit(1);
*/
  Div D1, D2;
/*
  D1 = gendiv(ff);
  D2 = gendiv(ff);

  if(chkdiv(D1,ff)==-1){
    printf("vaka\n");
  exit(1);
  }else{
    printf("good! === ");
    printpol(o2v(D1.u));
    printpoln(o2v(D1.v));
    //exit(1);
  }
  if(chkdiv(D2,ff)==-1){
    printf("vaka\n");
  exit(1);
  }else{
    printf("good! === ");
    printpol(o2v(D2.u));
    printpoln(o2v(D2.v));
    //exit(1);
  }
  */
  V = g2add(ff, uu1 , uu2, vv1, vv2);
  // V=xgcd(uu1,uu2,2);
  printpol(o2v(V.u));
  printf(" =====u3\n");
  printpol(o2v(V.v));
  printf(" =====v3\n");
  // printf("%d\n",equ(5,8));
  exit(1);
  if(chkdiv(V,ff)==- 1){
  printf("baka1^\n");
  exit(1);
  }
exit(1);
  Div C={0};
  // below undercondtruction
  k = odiv(osub(ff, (omul(D1.v, D1.v))), D1.u);
  s = omod(odiv(k, scr(2, D1.v)), D1.u);
  l = omul(s, D1.u);
  u3 = omod(osub(omul(s, s), (osub(scr(2, omul(D1.v, s)), (k)))), D1.u);
  v3 = omod(minus(oadd(l, D1.v)), u3);
  C.u=u3;
  C.v=v3;
  if(chkdiv(C,ff)==-1){
  printf("baka\n");
  exit(1);
  }else{
    printf("That's True!\n");
  }
  printpol(o2v(u3));
  printf("======du3\n");
  printpol(o2v(v3));
  printf("======dv3\n");
  exit(1);

  mkmf();

  makefg(O);

  count = 0;
  c2 = 0;
  for (i = 0; i < O; i++)
  {
    if (gf[i] > 0)
      count++;
    if (fg[i] > 0)
    {
      c2++;
    }
    else
    {
      // printf("i=%d\n",i);
    }
  }

  printf("%d %d\n", count, c2);
  // exit(1);

  return 0;
}
