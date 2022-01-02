#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdint.h>
#include <strings.h>
#include <stdbool.h>
#include "struct.h"
#include "chash-p.c"
#include <NTL/ZZ.h>

#define O 6859 // 1331 //2197,4913,6859
#define K 5
//#define P 37
//#define J 1412  // https://eprint.iacr.org/2011/306.pdf  example.4

NTL_CLIENT

//ZZ P=to_ZZ("37");

// 20211231 GPL HyperElliptic Curve DLP (･∀･) ﾔｺﾋﾞﾔｰﾝ!!
// 院卒失業中(小飼弾と同い年)
// 本格的なおおきな素体上の曲線については、
// このプログラムを元にNTLの古いバージョンを使って作る予定。

ZZ PP = to_ZZ("100000000000000003"); // Harley's example
ZZ P=to_ZZ("10000000000000000051");
//ZZ P =to_ZZ("5000000000000000008503491");

ZZ c[K + 1];

// OP型からベクトル型への変換
vec o2v(OP f)
{
  vec a;
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

//配列からベクトル表現の多項式へ変換する
vec Setvec(int n)
{
  int i;
  vec v;

  for (i = 0; i < n; i++)
  {
    v.x[n - 1 - i] = c[i];
  }

  return v;
}

//配列の値を係数として多項式に設定する
OP setpol(ZZ f[], int n)
{
  OP g;
  vec a;
  int i;

  memset(c, 0, sizeof(c));
  // memcpy (c, f, n);
  for (i = 0; i < n; i++)
  {
    c[i] = (f[i]);
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
  assert(("baka\n", n >= 0));

  for (i = n; i > -1; i--)
  {
    if (a.x[i] > 0)
    {
      cout << a.x[i];
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
      cout << a.x[i];
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
ZZ
xtrace(OP f, ZZ x)
{
  int i, d;
  ZZ u=to_ZZ("0") , v = to_ZZ("1");

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

  return u % P;
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

OP minu(OP f)
{
  unsigned int i, j;
  vec e;

  j = deg(e=o2v(f));
  for (i = 0; i < j + 1; i++)
  {
    if (e.x[i] > 0)
      e.x[i] = P - e.x[i];
  }
  f=v2o(e);

  return f;
}


//リーディグタームを抽出(default)
oterm LT(OP f)
{
  int i, k;
  oterm t = {0};

  k = deg(o2v(f));
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
  vec v;
  OP g = {0};

  v = o2v(f);
  g = v2o(v);

  return g;
}

// 20200816:正規化したいところだがうまく行かない
//多項式の足し算
OP oadd(OP f, OP g)
{
  vec a , b , c ;
  int i, j, k, l = 0;
  OP h = {0}, f2 = {0}, g2 = {0};

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
  int i, k, j;
  OP h = {0};
  vec test;
  ZZ n;

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
  printf("qqqqqqqqqqqqqqqqqq\n");
  // assert (op_verify (f));
  // assert (op_verify (g));
  int i, count = 0, k, l;
  oterm t = {0};
  OP h = {0}, e = {0}, r = {0};
  vec c ;

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
  vec a, b , d;
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

    cout << a.x[i] << " - " << b.x[i] << " = " << d.x[i] << endl;
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

// nを法とする逆数
ZZ inv(ZZ a, ZZ n)
{
  ZZ d, x, s, q, r, t, gcd;
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
ZZ
equ(ZZ a, ZZ b)
{
  ZZ i = inv(a, P);

  return i * b;
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
  int i = 0, j, n, k;
  OP h = {0}, o = {0}, e = {0};
  oterm a, b = {0}, c = {0};

  n = LT(g).n;
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

  k = LT(g).n;
  b = LT(g);
  OP ll;

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

//　スカラーで割る
OP cdiv(ZZ a, OP f)
{
  vec v;
  int i, l;
  ZZ k;

  v = o2v(f);
  l = odeg(f);
  k = inv(a, P);
  printpol(o2v(f));
  //printf(" ==kokko %llu\n", a);
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
  int i = 0, j, n, k;
  OP h = {0}, e = {0}, tt = {0}, o = {0};
  oterm a, b = {0}, c = {0};

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
    e = cdiv(b.a, f);
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
    int u;
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
OP scr(ZZ d, OP f)
{
  int i, n;
  vec v ;

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
  ZZ e1;
  e1 = inv(LT(f).a, P);
  printf("e=%llu\n", e1);
  f = scr(e1, f);

  return f;
}

//　構造体ごとモニックにする
EX monic(EX X)
{
  ZZ e1;
  e1 = inv(LT(X.d).a, P);
  X.d = scr(e1, X.d);
  X.h = scr(e1, X.h);
  X.u = scr(e1, X.u);
  X.v = scr(e1, X.v);

  return X;
}

// ある多項式の倍数になっているか
ZZ isideal(OP f, OP g)
{
  ZZ a, b, c;
  OP h;

  a = inv(LT(f).a, P);
  f = scr(a, f);
  b = LT(g).a;
  f = scr(b, f);
  if (oequ(f, g) == 0)
    return b * a % P;

  return to_ZZ("-1");
}

//拡張ユークリッドアルゴリズム
EX xgcd(OP f, OP g)
{
  OP h[10] = {0}, ww[10] = {0}, *v, *u, T = {0};
  oterm a, b;
  int i = 0, j, flg = 0, k;
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
    printf("dame\n");
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
  vr = omod(minu(D3.v), ur);
  D.u = ur;
  D.v = vr;

  return D;
}

// 因子の加法
Div cadd(OP ff, OP uu1, OP uu2, OP vv1, OP vv2)
{
  EX V;
  Div D3, D, null = {0};
  vec vx;
  OP e1, e2, d1, s1, s2, s3, c1, c2, d, u;

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
  //exit(1);
  V = xgcd(oadd(vv1, vv2), V.d);
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
   //exit(1);

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
    vd = omod(minu(v), ud);
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
int bit(ZZ b, int i)
{
  int k = 1;

  if (((b & (1 << i)) >> i) % 2 == 1)
  {
    return 1;
  }
  else
  {
    return 0;
  }
}





ZZ pow_mod(ZZ x, ZZ n, ZZ p)
{
  if (n == 0)
    return to_ZZ("1");
  if ((n & 1)==1)
    return (pow_mod(x, n - 1, p) * x) % p;
  x = pow_mod(x, n / 2, p);
  return (ZZ)((x * x) % p);
}

/* Takes as input an odd prime p and n < p and returns r
 * such that r * r = n [mod p]. */
ZZ tonelli_shanks(ZZ n, ZZ p)
{
 ZZ s = to_ZZ("0");
 ZZ q = p - 1;
  while ((q & 1) == 0)
  {
    q /= 2;
    ++s;
  }
  if (s == 1)
  {
    ZZ r = pow_mod(n, (p + 1) / 4, p);
    if ((r * r) % p == n)
      return r;
    return to_ZZ("0");
  }
  // Find the first quadratic non-residue z by brute-force search
 ZZ z = to_ZZ("1");
  while (pow_mod(++z, (p - 1) / 2, p) != p - 1)
    ;
 ZZ c = pow_mod(z, q, p);
 ZZ r = pow_mod(n, (q + 1) / 2, p);
 ZZ t = pow_mod(n, q, p);
 ZZ m = s;
  while (t != 1)
  {
   ZZ tt = t;
   ZZ i = to_ZZ("0");
    while (tt != 1)
    {
      tt = (tt * tt) % p;
      ++i;
      if (i == m)
        return to_ZZ("0");
    }
   ZZ b = pow_mod(c, pow_mod(to_ZZ("2"), m - i - 1, p - 1), p);
   ZZ b2 = (b * b) % p;
    r = (r * b) % p;
    t = (t * b2) % p;
    c = b2;
    m = i;
  }
  if ((r * r) % p == n)
    return r;
  return to_ZZ("0");
}

ZZ root(ZZ a, ZZ p)
{
  ZZ c, b;

  cout << "p mod = " << a << "@" << p%4 << p << endl;
  
  if (p % 4 == 3 || p % 8 == 5)
  {
    if (p % 4 == 3)
    {
      b = (p + 1) / 4;
      c = pow_mod(a, b, p);
      {
        //printf("good\n");
        cout << c << endl;
        return c;
      }
    }
    if (p % 8 == 5)
    {
      c = pow_mod(a, (p + 3) / 8, p);
      if ((c * c) % p != a)
        printf("baka2\n");
        return to_ZZ("-1");
    }
    if (c * c % p == a)
    {
      printf("good\n");
      cout << c << endl;
      return c;
    }
  }
  if (p % 8 == 5)
  {
    c = 2 * a * pow_mod(4 * a, (p - 5) / 8, p);
    if (c * c % p != a)
    {
      printf("dangerous\n");
      return to_ZZ("-1");
    }
    if (c * c % p == a)
    {
      printf("good\n");
      cout << c << endl;
      return c;
    }
    if (p % 8 == 1){
      c = tonelli_shanks(a, p);
      if(c*c%p != a){
        printf("fail!\n");
        return to_ZZ("-1");
      }else{
      return c;
      }
    }
    return to_ZZ("0");
  }

  return to_ZZ("-1");
}

// 曲線に代入した値を計算する
PO tr1e(ZZ f4, ZZ f3, ZZ f2, ZZ f1, ZZ f0, ZZ p)
{
  ZZ  x, y, f, g;
  PO aa;

  while (1)
  {
    x = xor128(); // % p;
    //y = xor128(); // % p;
    f = (pow_mod(x, to_ZZ("5"), p) + (f4 * pow_mod(x, to_ZZ("4"), p)) % p + (f3 * pow_mod(x, to_ZZ("3"), p)) % p + (f2 * pow_mod(x, to_ZZ("2"), p)) % p + f1 * x + f0) % p;
    y = root(f, p);
    g = (y * y) % p;
    if (f == g)
    {
      aa.x = x;
      aa.y = y;
      //exit(1);
      return aa;
    }
  }
  //  return -1;
}


OP diviser(PO o, PO m)
{
  ZZ  t1[2][3], cc[2];
  vec c1 ,c2;
  int i, j, k;
  ZZ a,b;
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
  ZZ z;
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
  printf("%d %d\n", c1.x[0], c1.x[1]);
  // exit(1);
  f=v2o(c1);
  printpoln(o2v(f));
  //exit(1);

  return f;
}


OP genv(PO a,PO b){
ZZ f,l,m,n;
vec v;
OP g;

l=(a.x-b.x)%P;
if(l<0)
  l+=P;
m=(a.y-b.y)%P;
if(m<0)
  m+=P;
n=inv(l,P);
l=(m*n)%P;
f=(a.x*l-a.y)%P;
if(f<0)
  f+=P;
v.x[0]=l;
v.x[1]=f;
g=v2o(v);

return g;
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
v1.x[1]=1;
v2.x[1]=1;
while(1){
    a = tr1e(vx.x[4], vx.x[3], vx.x[2], vx.x[1], vx.x[0], P); // cofficient of function
    b = tr1e(vx.x[4], vx.x[3], vx.x[2], vx.x[1], vx.x[0], P); // cofficient of function

    v1.x[0] = P-a.x;
    printpol(v1);
    printf("ppppppppppppp\n");
    c = v2o(v1);

    v2.x[0] = P-b.x;
    d = v2o(v2);
    d2=diviser(a,b);
    printpoln(o2v(d2));
    //exit(1);

    //d2 = v2o(z1);
    d1 = omul(c, d);
    printpoln(o2v(d1));
    //exit(1);
    D.u = d1;
    D.v = d2;
    if(chkdiv(D,f)!=-1)
    {
      printf("line\n");
      return D;
    }else if (chkdiv(D, f) == -1)
  {
    printf("so buggy!\n");
    exit(1);
  }

  //exit(1);
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

  V = xgcd(D.u, scr(to_ZZ("2"), D.v));
  munford(V);

  a = odiv(omul(D.u, D.u), omul(V.d, V.d));
  b = odiv(oadd(omul(V.u, omul(D.u, D.v)), omul(V.v, oadd(omul(D.v, D.v), f))), V.d);

  while (deg(o2v(a)) > 2 || odeg(b) > odeg(a))
  {
    count++;
    printf("count=%d\n", count);
    if (count > 100)
      break;
    uu = odiv(osub(f, omul(b, b)), a);
    vv = omod(minu(b), uu);
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

Div tbl[1024] = {0};
ZZ tmp[1024];
// 演算テーブルを作る（繰り返し2乗法）
void mktbl(Div D, OP f)
{
  int i;
printf("begin\n");
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
  printf("end\n");
}

// 因子のスカラー倍
Div jac(ZZ n, OP f)
{
  int i, j = 0, tmp[1024] = {0};
  ZZ k;
  Div L = {0}, D, G;

printf("in jac\n");
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
      //printpoln(o2v(tbl[tmp[i]].u));
      //printpoln(o2v(tbl[tmp[i]].v));
      exit(1);
    }

    L = cadd(f, tbl[tmp[i]].u, L.u, tbl[tmp[i]].v, L.v);
    if (chkdiv(L, f) == -1)
    {
      printf("dame %d\n", i);
      /*
      printpoln(o2v(tbl[tmp[i]].u));
      printpoln(o2v(tbl[tmp[i]].v));
      printpoln(o2v(G.u));
      printpoln(o2v(G.v));
      printpoln(o2v(L.u));
      printpoln(o2v(L.v));
    */
      exit(1);
    }
    if (oequ(D.u, L.u) == 0 && oequ(D.v, L.v) == 0)
    {
      cout << "infinity devide!" <<  n << endl;
      exit(1);
    }
  }
  printpoln(o2v(D.u));

  return L;
}

 /* definition of a curve */

void HEC(){
/*
//p2 = 2 ** 128 - 173
ZZ P2=to_ZZ("340282366920938463463374607431768211283");
ZZ r2 = to_ZZ("115792089237316195429342203801033554170931615651881657307308068079702089951781");
ZZ ef[K+1] = {1,0,318258242717201709453901384328569236653, 75380722035796344355219475510170298006, 129416082603460579272847694630998099237, 143864072772599444046778416709082679388};

//2**127-1
P3=to_ZZ("170141183460469231731687303715884105727");
ZZ p=to_ZZ("170141183460469231731687303715884105727");
ZZ r = to_ZZ("28948022309329048848169239995659025138451177973091551374101475732892580332259");
//which is 254 bits. A possible degree 5 model is C : y2 = x**5 + f3*x**3 + f2*x**2 + f1x + f0, where
ZZ f3[K+1] = {1,0,34744234758245218589390329770704207149,  132713617209345335075125059444256188021, 90907655901711006083734360528442376758,  6667986622173728337823560857179992816};


// p127m = (2**63 - 27443) · 2**64 + 1.
ZZ P4=to_ZZ("170141183460468725497689688904659107841");
//f:y^2 = x^5 + 17
fa[K+1]={1,0,0,0,0,17};
ZZ r2 = to_ZZ("28948022309328876595115567994214488524823328209723866335483563634241778912751");


ZZ P5=to_ZZ("28=340282366920938463463374607431768186521");
// y2 = x5+ 37
fb[K+1]={1,0,0,0,0,37};
ZZ r3 = to_ZZ("115792089237316195401210495125503591471546519982099914586091636775415022457661");
*/


/*
  //find by Harley
//  @q = 10 ** 19 + 51
q=10000000000000000051;
  ZZ a[5] = [3141592653589793238, 4626433832795028841, 9716939937510582097, 4944592307816406286, 2089986280348253421]
  //#@u1_=[13131182302866750318,6953593084278582387]
  //#@v1_=[@q,0] #infinity
  ZZ uh = [8940387226809150403, 3838225076702837943]
  ZZ vh = [8035450087728851271, 1893861348804881148]
  // ff=x^5+314159265358979338*x^4+4626433832795028841*x^3+9716939937510582097*x^2+4944592307816406286*x+2089986280348253421
  ZZ u1 = [10027301878627002813, 9681764174062850433]
  //#616419646419685014=a*3542790122851877922+b
  //#0=a*6484511755775124891+b
  //#616419646419685014=a*7058278367076753082
  ZZ v1 = [9406915506559133975, 920961725690419616]
  ZZ u2 = [15109848135481867673, 5563304430399854240]
  //#2935061693073737419=a*5239897978117534135+b
  //#3524464046627319761=a*9869950157364333538+b
  //#589402353553582342=a*4630052179246799403
  ZZ v2 = [7250939689363649434, 6461431514924022130]
  ZZ J = 99999999982871020671452277000281660080
  //#? 7054215880371151972602291562049
  
 
  #find by Lange
  //# f=xx^5+153834295433461683634059*xx^3+1503542947764347319629935*xx^2+1930714025804554453580068*xx+790992824799875905266969
  ZZ p3 = 1932005208863265003490787
  ZZ F1 = [0, 153834295433461683634059, 1503542947764347319629935, 1930714025804554453580068, 790992824799875905266969]
  ZZ uu0 = [1594018975878036024296315, 52552598504459997856285]
  ZZ uu1 = [1791061143796384566472590, 160038959612724914387201]
  //#504894935863953268767725=a*106028591185649525291891+b
  //#704210062398295465154981=a*1487990384692386499004424+b
  ZZ vv0 = [1288294670775269897135356, 942599769284250370891960]
  //#1=a*704665787761008893641614+b
  //#824513992484349685277891=a*1086395356035375672830976+b
  ZZ vv1 = [325694573428709528176000, 1410279347067078324080410]
  ZZ p4 = 3713820117856140824697372689
  //# ff=x^5+241216435998068557682742515*x^3+553011586465186980114036462*x^2+1456621446251091989731057514*x+3440013483680364963850133535
  ZZ F0 = [0, 241216435998068557682742515, 553011586465186980114036462, 1456621446251091989731057514, 3440013483680364963850133535]
  ZZ uua = [3090907731099435637713212933, 3430279740253146837327450789]
  //#1090190095529845640563737560=a*901573235033529767913345809+b
  //#1607209110223949051233778495=a*2189334496065905869799867124+b
  ZZ vva = [1516936926385660062377077660, 177161035616247877668423903]
  //#1044171804289905858226438503=a*10486172923382208811538764+b
  //#928996521279782754969642877=a*1874123356916514825274712274+b
  //#3598644834846017721440577063=a*1863637183993132616463173510
  ZZ uub = [1884609529839897034086251038, 265492627929763696542013047]
  ZZ vvb = [515973747200726989346030903, 2866067948417124660466300132]
*/
//  #find by Gaudry
//  #fg=x^5+2682810822839355644900736*x^3+226591355295993102902116*x^2+2547674715952929717899918*x+4797309959708489673059350
//  @P = 5 * 10 ** 24 + 8503491
 ZZ P6 =to_ZZ("5000000000000000008503491");
 ZZ FF[K+1] = {to_ZZ("1"),to_ZZ("0"), to_ZZ("2682810822839355644900736"), to_ZZ("226591355295993102902116"), to_ZZ("2547674715952929717899918"), to_ZZ("4797309959708489673059350")};
 ZZ Jga = to_ZZ("24999999999994130438600999402209463966197516075699");
 ZZ Jgb = to_ZZ("25000000000005869731468829402229428962794965968171");
  //#50724386855111482309402*2779199501981512279739817%P
  //#2055622596816515886446193*1553122609714208136553134%P
  //#1520505942073936921231867
  ZZ ug0[K+1] = {to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("1"),to_ZZ("1"), to_ZZ("-1713538969626908355896596")};
  ZZ vg0[K+1] = {to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("0"), to_ZZ("138905579055173741542118")};
  ZZ ug1[K+1] = {to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("1"),to_ZZ("1738366273424896804842766"), to_ZZ("3184841659043138633535652")};
  ZZ vg1[K+1] = {to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("2931056213155642836850986"), to_ZZ("402980510097415333052905")};

}


// 例
int main()
{
  unsigned int i, count = 0;


  ZZ f[K+1] = {to_ZZ("1"),to_ZZ("3141592653589793238"), to_ZZ("4626433832795028841"), to_ZZ("9716939937510582097"), to_ZZ("4944592307816406286"), to_ZZ("2089986280348253421")};
  //#@u1_=[13131182302866750318,6953593084278582387]
  //#@v1_=[@q,0] #infinity
  ZZ u1[K+1] = {to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("1"),to_ZZ("8940387226809150403"), to_ZZ("3838225076702837943")};
  ZZ v1[K+1] = {to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("8035450087728851271"), to_ZZ("1893861348804881148")};
  // ff=x^5+314159265358979338*x^4+4626433832795028841*x^3+9716939937510582097*x^2+4944592307816406286*x+2089986280348253421
 ZZ J = to_ZZ("99999999982871020671452277000281660080");
 
/*
  ZZ uh[K+1] = {0,0,0,1,10027301878627002813, 9681764174062850433};
  //#616419646419685014=a*3542790122851877922+b
  //#0=a*6484511755775124891+b
  //#616419646419685014=a*7058278367076753082
  ZZ vh[K+1] = {0,0,0,0,9406915506559133975, 920961725690419616};
  ZZ u2[K+1] = {0,0,0,1,15109848135481867673, 5563304430399854240};
  //#2935061693073737419=a*5239897978117534135+b
  //#3524464046627319761=a*9869950157364333538+b
  //#589402353553582342=a*4630052179246799403
  ZZ v2[K+1] = {0,0,0,0,7250939689363649434, 6461431514924022130};
*/
//  #find by Gaudry
//  #fg=x^5+2682810822839355644900736*x^3+226591355295993102902116*x^2+2547674715952929717899918*x+4797309959708489673059350
//  @P = 5 * 10 ** 24 + 8503491
/*
 ZZ f[K+1] = {to_ZZ("1"),to_ZZ("0"), to_ZZ("2682810822839355644900736"), to_ZZ("226591355295993102902116"), to_ZZ("2547674715952929717899918"), to_ZZ("4797309959708489673059350")};
 ZZ Jga = to_ZZ("24999999999994130438600999402209463966197516075699");
 ZZ Jgb = to_ZZ("25000000000005869731468829402229428962794965968171");
  //#50724386855111482309402*2779199501981512279739817%P
  //#2055622596816515886446193*1553122609714208136553134%P
  //#1520505942073936921231867
  ZZ u1[K+1] = {to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("1"),to_ZZ("1738366273424896804842766"), to_ZZ("3184841659043138633535652")};
  ZZ v1[K+1] = {to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("2931056213155642836850986"), to_ZZ("402980510097415333052905")};
*/

/*
   char f[K + 1] = {1, 7, 6, 2, 8, 2};
  
    char  u2[K + 1] = {0, 0, 0, 1, 21, 16};
    char  u1[K + 1] = {0, 0, 0, 1, 19, 20};
    char  v2[K + 1] = {0, 0, 0, 0, 21, 21};
    char  v1[K + 1] = {0, 0, 0, 0, 12, 8};
  */
  //char f[K + 1] = {1, 0, 2, 30, 5, 1};
/*  
  ZZ  f[K+1]= {to_ZZ("1"), to_ZZ("1597") , to_ZZ("1041") ,to_ZZ("5503") , to_ZZ("6101") , to_ZZ("1887") };
  //f1 = x + 28555025517563816 and f2 = x + 74658844563359755 ;
  ZZ  u2[K+1]={to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("1"),to_ZZ("1571353025997967")  , to_ZZ("12198441063534328")};
  ZZ  v2[K+1]={to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("32227723250469108") , to_ZZ("68133247565452990")};
  ZZ  u1[K+1]={to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("1"),to_ZZ("70887725815800572") , to_ZZ("94321182398888258")};
  ZZ  v1[K+1]={to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("42016761890161508") , to_ZZ("3182371156137467") };

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
  OP ff, k, uu1, uu2, vv1, vv2, s,o, d, c,m,d1,d2;
//  unsigned long long tst1[K + 1] = {0, 0, 0, 0, 8, 7};
//  unsigned long long tst2[K + 1] = {0, 0, 0, 0, 0, 10};
 // unsigned tmp[2][3] = {0};
  EX V;
  Div D;
  oterm a;
  OP b = {0};
  // unsigned long long  q1 = to_unsigned long long ("1208925819614629174708801");
  int a1 = 1331;
  // J1 =to_unsigned long long ("1461501637326815988079848163961117521046955445901");
  // e y2 = x5+a, a ∈ Fp

  // unsigned long long  q2 = to_unsigned long long ("1208925819614629174709941");
  int a2 = 2;
  // J2 = to_unsigned long long ("1461501637331762771847359428275278989652932675771");

  vec vx , xv ,v11,v22;
  Div G0, G1, X;
  PO a11,b11;
  ZZ x,y;

  ff = setpol(f, K + 1);
  printpoln(o2v(ff));
  //exit(1);
  ZZ A=root(f[5],P);
  cout << A << endl;
  //exit(1);
  vx=o2v(ff);
  a11=tr1e(vx.x[4],vx.x[3],vx.x[2],vx.x[1],vx.x[0],P);
  x=a11.x;
  y=a11.y;
  cout << "f= " << (x*x*x*x*x+vx.x[4]*x*x*x*x+vx.x[3]*x*x*x+vx.x[2]*x*x+vx.x[1]*x+vx.x[0])%P << endl;
  cout << " g= " << (y*y)%P << endl;
  cout << "x= " << x << "y= " << y << endl;
  //exit(1);
srand(clock());
/*
v11.x[1]=1;
v22.x[1]=1;
while(1){
    a11 = tr1e(vx.x[4], vx.x[3], vx.x[2], vx.x[1], vx.x[0], P); // cofficient of function
    b11 = tr1e(vx.x[4], vx.x[3], vx.x[2], vx.x[1], vx.x[0], P); // cofficient of function

    v11.x[0] = P-a11.x;
    printpol(v11);
    printf("ppppppppppppp\n");
    c = v2o(v11);

    v22.x[0] = P-b11.x;
    d = v2o(v22);
    d2=diviser(a11,b11);
    printpoln(o2v(d2));
    //exit(1);

    //d2 = v2o(z1);
    d1 = omul(c, d);
    printpoln(o2v(d1));
    //exit(1);
    D.u = d1;
    D.v = d2;
    if(chkdiv(D,ff)!=-1)
    {
      printf("line\n");
      exit(1);
    }
  exit(1);
}
*/
  uu1 = setpol(u1, K + 1);
  //uu2 = setpol(u2, K + 1);
  vv1 = setpol(v1, K + 1);
  //vv2 = setpol(v2, K + 1);
  //o = setpol(tst1, K + 1);
  //m = setpol(tst2, K + 1);

/*  
X=cadd(ff,uu1,uu2,vv1,vv2);
printf("%d\n",chkdiv(X,ff));
X.u=uu1;
X.v=vv1;
//exit(1);
V=xgcd(uu1,uu2);
V=monic(V);
munford(V);
*/
//exit(1);
  X.u=uu1;
  X.v=vv1;
if(LT(omod(osub(ff,omul(vv1,vv1)),uu1)).a==0)
printf("seikou\n");
//exit(1);

  //　ランダムな因子をヤコビ多様体の位数倍して無限遠点になれば正しい
  srand(clock());
   X = gendiv(ff);
   printf("%d\n",chkdiv(X,ff));
   //exit(1);

  mktbl(X, ff);

  X = jac(J+1, ff);
  if (chkdiv(X, ff) == -1)
  {
    printf("bakayo\n");
    // break;
  }
  exit(1);

  /*
    // V=xgcd(uu1,uu2);
    G0 = cadd(ff, uu1, uu2, vv1, vv2);
    printf("%d\n", chkdiv(G0, ff));
    // exit(1);
    srand(clock());
    G1.u = uu1;
    G1.v = vv1;
    X.u = uu2;
    X.v = vv2;
    G0 = gendiv(ff);

    while (1)
    {
        G0 = cdbl(G0, ff);
      if (chkdiv(G0, ff) == -1)
      {
        break;
      }
      else
      {
        printf("イイっ！この因子すげえいいっ！\n");
      }
    }
    exit(1);

    if (chkdiv(G1, ff) == -1 || chkdiv(X, ff) == -1)
    {
      printf("erro!\n");
      exit(1);
    }
    G0 = cadd(ff, G1.u, X.u, G1.v, X.v);
    if (chkdiv(G0, ff) == -1)
    {
      printpoln(o2v(G1.u));
      printpoln(o2v(G1.v));
      printpoln(o2v(X.u));
      printpoln(o2v(X.v));
      printpoln(o2v(G0.u));
      printpoln(o2v(G0.v));
      printf("bug\n");
      exit(1);
    }
  */
/*
  PO xx;
  unsigned long long rr = 0;
  count = 0;
  int xount = 0;

  G1 = gendiv(ff);
  X = gendiv(ff);
  uu1 = G1.u;
  vv1 = G1.v;
  srand(clock());
  while (1)
  {
    G1 = gendiv(ff);
    X = gendiv(ff);
    G0 = cadd(ff, G1.u, X.u, G1.v, X.v);

    if (chkdiv(G0, ff) == -1)
    {
      printf("baka\n");
      printpoln(o2v(G1.u));
      printpoln(o2v(G1.v));
      printpoln(o2v(X.u));
      printpoln(o2v(X.v));
      printpoln(o2v(G0.u));
      printpoln(o2v(G0.v));
      count++;
      exit(1);
      V = xgcd(X.u, G1.u);
      if (LT(V.d).n > 0)
      {
        printf("gcd!\n");

        // exit(1);
      }
    }
    else if (oequ(G1.u, uu1) == 0)
    {
      printf("order #J= %d\n", xount);
      exit(1);
    }
    else if (chkdiv(G1, ff) != -1)
    {
      printf("ウホッ！いい因子。\n");
      xount++;
      // exit(1);
    }
    else
    {
      printpoln(o2v(G0.u));
      printpoln(o2v(G0.v));
      printpoln(o2v(G1.u));
      printpoln(o2v(G1.v));
      printpoln(o2v(X.u));
      printpoln(o2v(X.v));
      printf("why?\n");
      // count++;
      exit(1);
    }
    if (count > 100)
      break;
    printf("%u xount=%u\n", count, xount);
  }
  printf("%u %u\n", count, xount);
  // exit(1);
*/

  return 0;
}