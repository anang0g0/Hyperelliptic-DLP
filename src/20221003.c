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
#define P 31
#define J 1412  // https://eprint.iacr.org/2011/306.pdf  example.4


// 20211231 GPL HyperElliptic Curve DLP (･∀･) ﾔｺﾋﾞﾔｰﾝ!!
// 院卒失業中(小飼弾と同い年)
// 本格的なおおきな素体上の曲線については、
// このプログラムを元にNTLの古いバージョンを使って作る予定。

//unsigned long long PP = 100000000000000003LLU; // Harley's example

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
  int i;
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
      printf("[%d] %llux^%llu\n", i, f.t[i].a, f.t[i].n);
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

OP init_pol(OP f){
  int i=0;
  vec v=o2v(f);

   memset(v.x,0,sizeof(v.x));

  f=v2o(v);

return f;
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
    // printf("%llu,",f[i]);
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
  assert((n >= 0));

  for (i = n; i > -1; i--)
  {
    if (a.x[i] > 0)
    {
      printf("%llu", a.x[i]);
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
  assert((n >= 0));

  for (i = n; i > -1; i--)
  {
    if (a.x[i] > 0)
    {
      printf("%llu", a.x[i]);
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

      // printf("\nv=%llu",v);
    }
    u = (u + v);
  }
  // printf("u=%llu\n",u%O);

  return u % O;
}

//多項式の次数(degのOP型)
int odeg(OP f)
{
  int i, j = 0;

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
  int i;
  oterm t = {0};

//  k = deg(o2v(f));
  for (i = 0; i < DEG; i++)
  {
    // printf("a=%llu %llu\n",f.t[i].a,f.t[i].n);
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
  OP h = {0};

  // for(i=0;i<257;i++)
  //  printf("%llu %llu %llu %llu %llu\n",i,f.t[i].a,f.t[i].n,g.t[i].a,g.t[i].n);

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
  if (LT(f).a == 0 && LT(g).a == 0)
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
  // printf("%llu %llu\n",i,b.x[i]);
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
  int i, k;
  OP h = {0};
  

  // f=conv(f);
  k = odeg(f);
  
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
  printf("qqqqqqqqqqqqqqqqqq\n");
  // assert (op_verify (f));
  // assert (op_verify (g));
  int i, k, l;
  oterm t = {0};
  OP h = {0}, e = {0};
  

  k = odeg(f);
  l = odeg(g);
  if (LT(f).a == 0 || LT(g).a == 0)
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
  if (LT(h).a == 0)
  {
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
      d.x[i] = (P + (a.x[i] - b.x[i])) % P;
    }

    printf("%llu - %llu = %llu\n", a.x[i], b.x[i], d.x[i]);
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
  int i;

  v = o2v(f);
  x = o2v(g);
  for (i = 0; i < 512; i++)
  {
    if (v.x[i] != x.x[i])
      return -1;
  }

  return 0;
}

// nを法とする逆数
unsigned long long inv(unsigned long long a, unsigned long long n)
{
  unsigned long long d, x, s, q, r, t;
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
  //gcd = d;

  return ((x + n) % (n / d));
}

// aに何をかけたらbになるか
unsigned long long
equ(unsigned long long a, unsigned long long b)
{
  unsigned long long i = inv(a, P);

  return (i * b) % P;
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
    s.a = (inv(t.a, P) * tt.a) % P; //(tt.a * inv(t.a, P)) % P;
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
  OP h = {0};
  oterm  b = {0}, c = {0};

  //n = LT(g).n;
  if (LT(f).a == 0)
  {
    return h;
  }
  if (LT(g).a == 0)
  {
    printf("g de dib 0!\n");
    exit(1);
  }

  if (LT(f).n < LT(g).n)
  {
    //    exit(1);
    return f;
  }

  //k = LT(g).n;
  b = LT(g);


  // assert(("double baka\n", b.a > 0 && b.n > 0));
  while (LT(g).a != 0)
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
    if (LT(f).a == 0)
    {
      //      printf("blake1\n");
      break;
    }

    if (c.a == 0 || b.a == 0)
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
  int i, k;
  
  oterm t;

  t = LT(f);
  // f = conv(f);
  k = odeg((f)) + 1;
  for (i = 0; i < k; i++)
    f.t[i].a = (f.t[i].a * inv(t.a, P)) % P;

  return f;
}

//　スカラーで割る
OP cdiv(unsigned long long a, OP f)
{
  vec v;
  int i, l, k;

  v = o2v(f);
  l = odeg(f);
  k = inv(a, P);
  printpol(o2v(f));
  printf(" ==kokko %llu\n", a);
  // a=equ(a,LT(f).a);
  for (i = 0; i < l + 1; i++)
    v.x[i] = (k * v.x[i]) % P;
  printpol(v);
  // exit(1);
  f = v2o(v);

  return f;
}

//多項式の商を取る
OP odiv(OP f, OP g)
{

  f = conv(f);
  g = conv(g);
  // assert (op_verify (f));
  // assert (op_verify (g));
  int i = 0, k;
  OP h = {0}, e = {0}, tt = {0};
  oterm  b = {0}, c = {0};
  vec vx={0};

  printpol(o2v(f));
  printf("@@@@@@f\n");
  printpol(o2v(g));
  printf("@@@@@@@@@@g\n");
  // exit(1);

  if (LT(f).a == 0)
  {
    return h;
  }
  if (LT(g).a == 0)
  {
    printf("baka^\n");
    exit(1);
  }

  k = 0; // odeg(f) - odeg(g);
  b = LT(g);

  if (b.a > 0 && b.n == 0)
  {
    //vx.x[0] = (inv(LT(g).a,P)*LT(f).a)%P;
    //e=o2v(vx);
    e=cdiv(b.a, f);
    printpol(o2v(e));
    printf(" cdiv\n");
    // exit(1);
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
  k = 0;
  while (LT(g).a != 0)
  {
    c = LTdiv(f, b);
    c.a = c.a % P;
    assert(c.n < DEG);
    if (c.a > 0)
    {
      tt.t[k] = c;
      k++;
    }
    printf("%llu", c.a);
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
    if (LT(f).a == 0)
    {
      printf("blake2\n");
      // c.a=1;
      break;
    }
    
    if ((oequ(f, g)) == 0)
    {
      printpol(o2v(tt));
      printf("\n");
      c.a = 1;
      // break;
      //  exit(1);
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

// スカラー倍する
OP scr(unsigned long long d, OP f)
{
  int i, n;
  vec v = {0};

  n = deg(o2v(f));
  v = o2v(f);
  for (i = 0; i < n + 1; i++)
  {
    printf("v[%d]=%llu\n", i, v.x[i]);
    v.x[i] = (v.x[i] * d) % P;
  }
  f = v2o(v);

  return f;
}

// モニック多項式にする
OP monique(OP f)
{
  unsigned long long e1;
  e1 = inv(LT(f).a, P);
  printf("e=%llu\n", e1);
  f = scr(e1, f);

  return f;
}

//　構造体ごとモニックにする
EX monic(EX X)
{
  int e1;
  e1 = inv(LT(X.d).a, P);
  X.d = scr(e1, X.d);
  X.h = scr(e1, X.h);
  X.u = scr(e1, X.u);
  X.v = scr(e1, X.v);

  return X;
}

// ある多項式の倍数になっているか
int isideal(OP f, OP g)
{
  int a, b;

  a = inv(LT(f).a, P);
  f = scr(a, f);
  b = LT(g).a;
  f = scr(b, f);
  if (oequ(f, g) == 0)
    return b * a % P;

  return -1;
}

//拡張ユークリッドアルゴリズム
EX xgcd(OP f, OP g)
{
  OP h[10] = {0}, ww[10] = {0}, *v, *u;
  int i = 0, k;
  EX e = {0};
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
   //exit(1);
  if (LT(g).a == 0)
  {
    printf("g ==0\n");
    exit(1);
  }

  k = odeg(g);
  if (k == 0 && LT(g).a != 0)
  {
    printf("use cdiv\n");
    // f=cdiv(LT(g).a,f);
    // printpoln(o2v(g));
    // exit(1);
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

    if (LT(g).a == 0)
    {
      printpol(o2v(g));
      printf("ggggggggggggg\n");
    }
    if (LT(g).a == 0)
      break;
    i++;
  //exit(1);
  }
  // f=g;
  // g=h[i];
  if (LT(g).a != 0)
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
  printpoln(o2v(e.d));
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

// uu2を法とするuu1の逆元
OP qinv(OP uu1, OP uu2)
{
  EX tt;
  OP v;
  

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
  v = scr(LT(tt.d).a, tt.u);
  printpol(o2v(v));
  printf(" in qinv's v\n");

  // exit(1);
  if (LT(omod(omul(v, uu1), uu2)).a == 0)
  {
    return v;
  }
  printf("can't\n");
  // exit(1);
}

// 因子の条件をチェック
int chkdiv(Div d, OP f)
{
  OP t;

  t = omod(osub(omul(d.v, d.v), (f)), d.u);
  printpol(o2v(t));
  printf(" 00000000000 t\n");
  if (LT(t).a == 0)
  {
    return 1;
  }
  else
  {
    printpoln(o2v(d.u));
    printpoln(o2v(d.v));
    // exit(1);
  }
  return -1;
}
// 使うかもしれない
Div reduce(Div D3, OP f)
{
  Div D;
  OP ur, vr;

  ur = osub(f, odiv(omul(D3.v, D3.v), D3.u));
  vr = omod(minus(D3.v), ur);
  D.u = ur;
  D.v = vr;

  return D;
}

// 因子の加法
Div cadd(OP ff, OP uu1, OP uu2, OP vv1, OP vv2)
{
  EX V;
  Div D3;
  OP e1, e2, s1, s2, s3, c1, c2, d, u,d1;

  V = xgcd(uu1, uu2);
  V = monic(V);
  e1 = V.u;
  printpol(o2v(V.u));
  printf("  e1\n");
  e2 = V.v;
  printpol(o2v(V.v));
  printf("  e2\n");
  d1 = V.d;
  printpol(o2v(V.d));
  printf("  d1\n");

  printpol(o2v(V.h));
  printf("  Uh\n");
  // exit(1);
  V = xgcd(oadd(vv1, vv2), d1);
  V = monic(V);
  c2 = V.u;
  printpol(o2v(V.u));
  printf("  c2\n");
  c1 = V.v;
  printpol(o2v(V.v));
  printf("  c1\n");
  d = V.d;
  printpol(o2v(V.d));
  printf("  d\n");

  printpol(o2v(V.h));
  printf("  Uh\n");
  // exit(1);

  s1 = omul(c1, e1);
  printpol(o2v(s1));
  printf(" ==s1\n");
  s2 = omul(c1, e2);
  printpol(o2v(s2));
  printf(" ==s2\n");
  s3 = c2;
  printpol(o2v(s3));
  printf(" ==s3\n");
  // exit(1);
  int count = 0;
  OP v;
  Div D1;

  // u=odiv(omul(uu1,uu2),omul(d,d));
  // u=omul(d,d);
  u = odiv(omul(uu1, uu2), omul(d, d));
  printpol(o2v(u));
  printf(" ==u3@\n");
  // exit(1);

  count++;
  v = omod(odiv(oadd(oadd(omul(omul(s1, uu1), vv2), omul(omul(s2, uu2), vv1)), omul(s3, oadd(omul(vv1, vv2), ff))), d), u);
  printpol(o2v(v));
  printf(" ==vu3@\n");
  // exit(1);
  D1.u = u;
  D1.v = v;
  printf("%d\n", chkdiv(D1, ff));
  // exit(1);
  OP ud, vd;
  if (odeg(u) > 2)
  {

  reduct:
    printpol(o2v(u));
    printf(" =======UUUUUUUU\n");
    // printpoln(o2v(odiv(osub(ff,omul(v,v)),u)));
    // printpoln(o2v(u));
    // exit(1);
    ud = odiv(osub(ff, omul(v, v)), u);
    // exit(1);
    printpoln(o2v(u));
    printpoln(o2v(v));
    printpoln(o2v(ud));
    printpoln(o2v(ff));
    // exit(1);
    vd = omod(minus(v), ud);
    D1.u = monique(ud);
    D1.v = vd;

    printf("%d\n", chkdiv(D1, ff));
    // exit(1);
    if (odeg(ud) > 2)
    {
      if (count > 100)
      {
        printf("over 100\n");
        exit(1);
      }
      u = ud;
      v = vd;
      printf("==================\n");
      goto reduct;
    }
  }
  else
  {
    ud = u;
    vd = v;
  }
  ud = monique(ud);
  printpol(o2v(ud));
  printf(" @@ud\n");
  printpol(o2v(vd));
  printf(" @@udv\n");
  // printpoln(o2v(oadd(vv1,vv2)));
  // exit(1);
  D3.u = ud;
  D3.v = vd;
  /*
  printf("debug point\n");
  if(chkdiv(D3,ff)==1){
  return D3;
  }else if(chkdiv(D3,ff)==-1){
  printf("ptr\n");
  V=xgcd(ud,vd);
  if(odeg(V.d)==0){
    printf("kasu\n");
    exit(1);
  }
  }
  */

  return D3;
}

// bit count
int bit(unsigned long long b, int i)
{

  if (((b & (1 << i)) >> i) % 2 == 1)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

// 誰かが書いたやつ (https://gist.github.com/LaurentMazare/6745649)
unsigned long long pow_mod(__int128_t x, __int128_t n, __int128_t p)
{
  if (n == 0)
    return 1;
  if (n & 1)
    return (pow_mod(x, n - 1, p) * x) % p;
  x = pow_mod(x, n / 2, p);
  return (unsigned long long)((x * x) % p);
}

/* Takes as input an odd prime p and n < p and returns r
 * such that r * r = n [mod p]. */
unsigned long long tonelli_shanks(unsigned long long n, unsigned long long p)
{
  unsigned long long s = 0;
  unsigned long long q = p - 1;
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
  unsigned long long z = 1;
  while (pow_mod(++z, (p - 1) / 2, p) != p - 1)
    ;
  unsigned long long c = pow_mod(z, q, p);
  unsigned long long r = pow_mod(n, (q + 1) / 2, p);
  unsigned long long t = pow_mod(n, q, p);
  unsigned long long m = s;
  while (t != 1)
  {
    unsigned long long tt = t;
    unsigned long long i = 0;
    while (tt != 1)
    {
      tt = (tt * tt) % p;
      ++i;
      if (i == m)
        return 0;
    }
    unsigned long long b = pow_mod(c, pow_mod(2, m - i - 1, p - 1), p);
    unsigned long long b2 = (b * b) % p;
    r = (r * b) % p;
    t = (t * b2) % p;
    c = b2;
    m = i;
  }
  if ((r * r) % p == n)
    return r;
  return 0;
}

// cf.A cource in computational number theory p.31（素数を法とする平方根を求める）
unsigned long long root(unsigned long long a, unsigned long long p)
{
  __int128_t c, b;

  printf("p mod = %llu == %llu , %llu\n", a, p % 4, p);
//  if (p % 4 == 3 || p % 8 == 5)
  {
    if (p % 4 == 3)
    {
      b = (p + 1) / 4;
      c = pow_mod(a, b, p);
      {
        // printf("good\n");
        printf("good c= %llu\n", (unsigned long long)c);
        return (unsigned long long)c;
      }
    }
    if (p % 8 == 5)
    {
      c = pow_mod(a, (p + 3) / 8, p);
      if ((c * c) % p != a){
        if (p % 8 == 5)
        {
          c = 2 * a * pow_mod(4 * a, (p - 5) / 8, p);
        if (c * c % p != a)
        {
          printf("dangerous\n");
          return -1;
        }
    if (c * c % p == a)
    {
      printf("good\n");
      printf("%llu\n", (unsigned long long)c);
      return c;
    }
      }
        printf("baka2\n");
      return -1;
    }
    if (c * c % p == a)
    {
      printf("good\n");
      printf("%llu\n", (unsigned long long)c);
      return (unsigned long long)c;
    }
  }  
    if (p % 8 == 1 || p%4==1)
    {
      c = tonelli_shanks(a, p);
      if (c * c % p != a)
      {
        printf("fail!\n");
        return -1;
      }
      else
      {
        return c;
      }
    }
    return 0;
  }

  return -1;
}

// 曲線に代入した値を計算する
PO tr1e(unsigned long long f4, unsigned long long f3, unsigned long long f2, unsigned long long f1, unsigned long long f0, unsigned long long p)
{
  unsigned long long x, y, f, g;
  PO aa = {0};

  while (1)
  {
    x = rand() % p;
    // y = rand() % p;
    f = (pow_mod(x, 5, P) + (f4 * pow_mod(x, 4, P)) % P + (f3 * pow_mod(x, 3, P)) % P + (f2 * pow_mod(x, 2, P)) % p + f1 * x + f0) % p;
    y = root(f, P);
    g = (y * y) % p;
    printf("%d %d\n", f , g); 
    printf("%d %d %d\n", x , y, f0); 
    //exit(1);

    if (f == g)
    {
      aa.x = x;
      aa.y = y;
      return aa;
    }
  }
  //  return -1;
}


OP diviser(PO o, PO m)
{
  unsigned long long t1[2][3], cc[2];
  vec c1 ,c2;
  int i, j, k;
  unsigned long long a,b;
  OP f;

  t1[0][0] = o.x;//o.t[1].a;
  t1[1][0] = m.x; //t[0].a;
  t1[0][1] = 1; //t[1].a;
  t1[1][1] = 1; //t[0].a;
  t1[0][2] = o.y;
  t1[1][2] = m.y;


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
  printf("%llu\n", cc[0]);
  unsigned long long z;
  z = t1[1][0];
  printf("z=%d\n",z);
  
cc[0]=(cc[0]*z)%P;
printf("c0=%d\n",cc[0]);
//exit(1);

  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < 3; j++)
      printf("%llu,", t1[i][j]);
    printf("\n");
  }
  printf("\n");

  for (j = 0; j < 3; j++)
  {
    t1[1][j] =((t1[0][j] * cc[0])%P - t1[1][j] ) % P;
    printf("%llu,", t1[1][j]);
  }
  printf("\n");
  b=(t1[1][2]*inv(t1[1][1],P)%P);
  printf("b=%d\n",b);
  for(i=0;i<3;i++)
  a=(t1[0][2]-(t1[0][1]*b)%P);
  if(a<0)
  a+=P;
  a=(inv(t1[0][0],P)*a)%P;
  printf("a=%d\n",a);
 
  c1.x[0] = b;
  c1.x[1] = a;
  //cout << c1.x[0] << " , " <<  c1.x[1] << endl;
  // exit(1);
  f=v2o(c1);
  printpoln(o2v(f));
  //exit(1);

  return f;
}



// ランダムな因子の生成
Div gendiv(OP f)
{
  int count=0;
  PO a , b , e ;
  OP d1 = {0}, d2 = {0}, c = {0}, d = {0},vv1={0},vv2={0},uu1,v;
  //  ZZ  x, y, i, j, k,

  Div D = {0};
  vec v1 , v2, z1, z2, ff,vx;
  EX V;

vx=o2v(f);
///srand(clock());
memset(v1.x,0,sizeof(v1.x));
memset(v2.x,0,sizeof(v2.x));
v1.x[1]=1;
v2.x[1]=1;
do
{
  do{
    a = tr1e(vx.x[4], vx.x[3], vx.x[2], vx.x[1], vx.x[0], P); // cofficient of function
    b = tr1e(vx.x[4], vx.x[3], vx.x[2], vx.x[1], vx.x[0], P); // cofficient of function
    count++;
    if(count>100){
      printf("try over 100\n");
    exit(1);
    }
  }while(a.y > P || b.y > P);

    printf("P= %d\n", P );
    printf("P-x= %d\n", P-a.x);
    printf("ax= %d\n", a.x);
    printf("ay= %d\n", a.y);
    printf("bx= %d\n", a.x);
    printf("by= %d\n", b.y);

    v1.x[0] = P-a.x;
    printpol(v1);
    printf("ppppppppppppp\n");
    //init_pol(c);
    //init_pol(d);
    c = v2o(v1);
    printpol(o2v(c));
    printf("--------------c\n");

    v2.x[0] = P-b.x;
    d = v2o(v2);
    printpol(o2v(d));
    printf("--------------d\n");

    d2=diviser(a,b);
    printpol(o2v(d2));
    printf("============d2\n");
    //exit(1);

    //d2 = v2o(z1);
    d1 = omul(c, d);
    printpol(o2v(d1));
    printf("==============d1\n");
    //exit(1);
    D.u = d1;
    D.v = d2;
    if(chkdiv(D,f)!=-1)
    {
      printf("buggy\n");
      //exit(1);
      //return D;
    }

  //exit(1);
}while(chkdiv(D,f)== -1);

if (chkdiv(D, f) == -1)
  {
    printf("so buggy!\n");
    exit(1);
  }
  return D;
}


// test function
EX munford(EX V)
{

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

//　2倍点の計算
Div cdbl(Div D, OP f)
{
  Div D2;
  EX V;
  OP a, b, uu, vv;
  int count = 0;

  V = xgcd(D.u, scr(2, D.v));
  munford(V);

  a = odiv(omul(D.u, D.u), omul(V.d, V.d));
  b = odiv(oadd(omul(V.u, omul(D.u, D.v)), omul(V.v, oadd(omul(D.v, D.v), f))), V.d);

  while (deg(o2v(a)) > 2) // || odeg(b) > odeg(a))
  {
    count++;
    printf("count=%d\n", count);
    if (count > 10)
      break;
    uu = odiv(osub(omul(b, b),f), a);
    vv = omod(minus(b), uu);
    a = uu;
    b = vv;
  }
  printpol(o2v(a));
  printf("aaaaaaaaaaaaaaaa\n");
  printpol(o2v(b));
  printf("bbbbbbbbbbbbbbbbb\n");

  // if(deg(o2v(a))>0)
  // exit(1);
  D2.u = a;
  D2.v = b;

  // printf("chk==%d\n",chkdiv(D2,f));

  return D2;
}

Div tbl[640] = {0};

// 演算テーブルを作る（繰り返し2乗法）
void mktbl(Div D, OP f)
{
  int i;

  tbl[0] = D;
  if (chkdiv(D, f) == -1)
    exit(1);
  for (i = 0; i < 256; i++)
  {
    tbl[i + 1] = cdbl(tbl[i], f);
    if (chkdiv(tbl[i + 1], f) == -1)
    {
      printf("bakayo\n");
      exit(1);
    }
  }
}

// 因子のスカラー倍
Div jac(unsigned long long n, OP f)
{
  int i, j = 0, tmp[640]= {0}, k;
  Div L = {0}, D, G;

  G=gendiv(f);
//mktbl(G,f,tbl);
  k = n;
  i = 0;
  while (k > 0)
  {
    if (k % 2 == 1)
    {
      tmp[j++] = i;
      // printf("i=%d\n",i);
    }
    k = (k >> 1);
    i++;
  }
  for (i = 0; i < j; i++)
  {
    if (chkdiv(tbl[tmp[i]], f) == -1)
    {
      printf("tbl is bad %d\n", i);
      exit(1);
    }
  }
  L = tbl[tmp[0]];
  D = L;
  // printf("j=%d\n",j);
  for (i = 1; i < j; i++)
  {
    G = L;
    if (chkdiv(tbl[tmp[i]], f) == -1)
    {
      printf("before\n");
      printpoln(o2v(tbl[tmp[i]].u));
      printpoln(o2v(tbl[tmp[i]].v));
      exit(1);
    }

    L = cadd(f, tbl[tmp[i]].u, L.u, tbl[tmp[i]].v, L.v);
    if (chkdiv(L, f) == -1)
    {
      printf("dame %d\n", i);
      printpoln(o2v(tbl[tmp[i]].u));
      printpoln(o2v(tbl[tmp[i]].v));
      printpoln(o2v(G.u));
      printpoln(o2v(G.v));
      printpoln(o2v(L.u));
      printpoln(o2v(L.v));

      exit(1);
    }
    if (oequ(D.u, L.u) == 0 && oequ(D.v, L.v) == 0)
    {
      printf("infinity devide! %llu\n", n);
      exit(1);
    }
  }
  printpoln(o2v(D.u));

  return L;
}

// 例
int main()
{
  unsigned int i, count = 0;

  //unsigned long long f[K+1] = {1, 0,1184, 1846, 956, 560};


   unsigned long long f[K + 1] = {1, 7, 6, 2, 8, 2};
  /*
    unsigned long long  u2[K + 1] = {0, 0, 0, 1, 21, 16};
    unsigned long long  u1[K + 1] = {0, 0, 0, 1, 19, 20};
    unsigned long long  v2[K + 1] = {0, 0, 0, 0, 21, 21};
    unsigned long long  v1[K + 1] = {0, 0, 0, 0, 12, 8};
  */
  //unsigned long long f[K + 1] = {1, 0, 2, 30, 5, 1};
  /*
    //unsigned long long  f[K+1]= {1, 1597 , 1041 ,5503 , 6101 , 1887 };
  //f1 = x + 28555025517563816 and f2 = x + 74658844563359755 ;
  unsigned long long  u2[K+1]={0,0,0,1,1571353025997967 , 12198441063534328};
  unsigned long long  v2[K+1]={0,0,0,0,32227723250469108 , 68133247565452990};
  unsigned long long  u1[K+1]={0,0,0,1, 70887725815800572 , 94321182398888258};
  unsigned long long  v1[K+1]={0,0,0,0, 42016761890161508 , 3182371156137467 };

  unsigned long long u2[K + 1] = {0, 0, 0, 1, 26, 20};
  unsigned long long v2[K + 1] = {0, 0, 0, 0, 29, 26};
  unsigned long long u1[K + 1] = {0, 0, 0, 1, 9, 27};
  unsigned long long v1[K + 1] = {0, 0, 0, 0, 29, 16};

  unsigned long long u2[K + 1] = {0, 0, 0, 1, 30, 3};
  unsigned long long v2[K + 1] = {0, 0, 0, 0, 12, 8};
  unsigned long long u1[K + 1] = {0, 0, 0, 1, 0, 12};
  unsigned long long v1[K + 1] = {0, 0, 0, 0, 10, 4};
*/
  //  unsigned long long  f[K + 1] = {1, 0, 3, 7, 1, 2};
  /*
    unsigned long long  u2[K + 1] = {0, 0, 0, 1, 7, 10};
    unsigned long long  u1[K + 1] = {0, 0, 0, 1, 0, 10};
    unsigned long long  v2[K + 1] = {0, 0, 0, 0, 1, 9};
    unsigned long long  v1[K + 1] = {0, 0, 0, 0, 7, 9};
  */
  OP ff, k, uu1, uu2, vv1, vv2, s,o, d, c,m;
  unsigned long long tst1[K + 1] = {0, 0, 0, 0, 8, 7};
  unsigned long long tst2[K + 1] = {0, 0, 0, 0, 0, 10};
  unsigned tmp[2][3] = {0};
  EX V;
  oterm a;
  OP b = {0};
  // unsigned long long  q1 = to_unsigned long long ("1208925819614629174708801");
  int a1 = 1331;
  // J1 =to_unsigned long long ("1461501637326815988079848163961117521046955445901");
  // e y2 = x5+a, a ∈ Fp

  // unsigned long long  q2 = to_unsigned long long ("1208925819614629174709941");
  int a2 = 2;
  // J2 = to_unsigned long long ("1461501637331762771847359428275278989652932675771");

  vec vx = {0}, xv = {0};
  Div G0, G1, X;
  PO a11,b11;
  unsigned x,y;

  ff = setpol(f, K + 1);
  
  srand(clock());

  gendiv(ff);
//  exit(1);
  /*
  vx=o2v(ff);
  a11=tr1e(vx.x[4],vx.x[3],vx.x[2],vx.x[1],vx.x[0],P);
  b11=tr1e(vx.x[4],vx.x[3],vx.x[2],vx.x[1],vx.x[0],P);
  x=a11.x;
  y=a11.y;
  printf("f=%u\n",(x*x*x*x*x+vx.x[4]*x*x*x*x+vx.x[3]*x*x*x+vx.x[2]*x*x+vx.x[1]*x+vx.x[0])%P);
  printf("g=%u\n",(y*y)%P);
  //exit(1);

  uu1 = setpol(u1, K + 1);
  uu2 = setpol(u2, K + 1);
  vv1 = setpol(v1, K + 1);
  vv2 = setpol(v2, K + 1);
  o = setpol(tst1, K + 1);
  m = setpol(tst2, K + 1);

X=cadd(ff,uu1,uu2,vv1,vv2);
printf("%d\n",chkdiv(X,ff));
printpoln(o2v(X.u));
printpoln(o2v(X.v));
//exit(1);


printf("%d\n",chkdiv(G0,ff));
/exit(1);
X.u=uu1;
X.v=vv1;
V=xgcd(uu1,uu2);
V=monic(V);
munford(V);
exit(1);

PO ah,bh;
OP h,g;


ah.x=2;
bh.x=10;
ah.y=1;
bh.y=4;
//printf("%d\n",equ(8,21));
//exit(1);
//exit(1);
  //　ランダムな因子をヤコビ多様体の位数倍して無限遠点になれば正しい
i=1;
*/
X = gendiv(ff);
//exit(1);
if(chkdiv(X,ff)==-1){
  printf("end of baka\n");
  exit(1);
}
mktbl(X, ff);

for(i=1;i<P*P;i++)
{

  X = jac(i, ff);
  if (chkdiv(X, ff) == -1)
  {
    printf("baka\n");
     exit(1);
  }
}

  return 0;
}