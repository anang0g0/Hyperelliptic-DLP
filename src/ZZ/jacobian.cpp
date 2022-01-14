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

NTL_CLIENT

// ZZ P=to_ZZ("2003");
//  20211231 GPL HyperElliptic Curve DLP (･∀･) ﾔｺﾋﾞﾔｰﾝ!!
//  このプログラムを元にNTLの古いバージョンを使っています。

// ZZ P = to_ZZ("340282366920938463463374607431768186521"); //aa=2187
// ZZ P = to_ZZ("170141183460468725497689688904659107841"); //aa=17
 ZZ P = to_ZZ("2923003274661805836407369665432566039311865180529"); // default
// ZZ P = to_ZZ("340282366920938463463374607431768211283"); //omake
//ZZ P=to_ZZ("340282366920938463463374607431756119641");


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
  int i;
  OP f = {};

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
      cout << "[" << i << "] " << f.t[i].a << "x^" << f.t[i].n << endl;
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
  assert(( n >= 0));

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
  assert(( n >= 0));

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
ZZ xtrace(OP f, ZZ x)
{
  int i, d;
  ZZ u = to_ZZ("0"), v = to_ZZ("1");

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
  int i, j = 0;

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

  j = deg(e = o2v(f));
  for (i = 0; i < j + 1; i++)
  {
    if (e.x[i] > 0)
      e.x[i] = P - e.x[i];
  }
  f = v2o(e);

  return f;
}

OP init_pol(OP f)
{
  vec v = o2v(f);

  memset(v.x, 0, sizeof(v.x));

  f = v2o(v);

  return f;
}

//リーディグタームを抽出(default)
oterm LT(OP f)
{
  int i;
  oterm t = {};


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
  OP g = {};

  v = o2v(f);
  g = v2o(v);

  return g;
}

// 20200816:正規化したいところだがうまく行かない
//多項式の足し算
OP oadd(OP f, OP g)
{
  vec a, b, c;
  int i, j, k, l = 0;
  OP h = {}, f2 = {}, g2 = {};

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
  OP h = {};
  vec test;
  ZZ n;

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
  int i,  k, l;
  oterm t = {};
  OP h = {}, e = {}, r = {};
  vec c;

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
  vec a, b, d;
  int i, k, l, m;
  OP ans = {};

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
ZZ inv(ZZ a, ZZ n)
{
  ZZ d, x, s, q, r, t, gcd;
  d = n;
  x = to_ZZ("0");
  s = to_ZZ("1");
  while (a != 0)
  {
    q = d / a;
    r = d % a;
    d = a;
    a = r;
    t = (x - q * s) % n;
    x = s;
    s = t;
  }
  gcd = d;

  return ((x + n) % (n / d)) % n;
}

// aに何をかけたらbになるか
ZZ equ(ZZ a, ZZ b)
{
  ZZ i = inv(a, P);

  return i * b;
}

//多項式を単行式で割る
oterm LTdiv(OP f, oterm t)
{
  oterm tt = {}, s = {
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
  OP h = {}, o = {}, e = {};
  oterm a, b = {}, c = {};


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


  b = LT(g);
  OP ll;

  // assert(("double baka\n", b.a > 0 && b.n > 0));
  while (1)
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

    if (c.a == 0 || LT(f).a==0)
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
  int i,  k;
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
  // printf(" ==kokko %llu\n", a);
  //  a=equ(a,LT(f).a);
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
  OP h = {}, e = {}, tt = {}, o = {};
  oterm a, b = {}, c = {};

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
  OP null = {};
  if (odeg((f)) < odeg((g)))
  {
    return null;
  }

  i = 0;
  k = 0;
  while (1)
  {
    c = LTdiv(f, b);
    c.a = c.a % P;
    assert(c.n < DEG);
    if (c.a > 0)
    {
      tt.t[k] = c;
      k++;
    }

    cout << c.a << endl;
    printf(" ccccccccccccccccc\n");
    printpol(o2v(g));
    printf(" ===before_g in_odiv\n");
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
    /*
    if (LT(f).a == 0)
    {
      printf("blake2\n");
      // c.a=1;
      //break;
    }
    */
    if ((oequ(f, g)) == 0)
    {
      printpol(o2v(tt));
      printf("\n");
      c.a = 1;
      //break;
      //exit(1);
    }

    if (c.a == 0)
      break;
  }

  // tt は逆順に入ってるので入れ替える
  OP ret = {};

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
  vec v;

  n = deg(o2v(f));
  v = o2v(f);
  for (i = 0; i < n + 1; i++)
  {
    cout << "v[" << i << "]=" << v.x[i] << endl;
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
  cout << "e=" << e1 << endl;
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
  OP h[10] = {}, ww[10] = {}, v[256]={}, u[256]={}, T = {};
  oterm a, b;
  int i = 0,  k;
  EX e = {}, ee = {};


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
      break;
    }
    i++;
    // exit(1);
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
   exit(1);
}

// 因子の条件をチェック
int chkdiv(Div d, OP f)
{
  OP t;

  t = omod(osub(omul(d.v, d.v), (f)), d.u);
  printpol(o2v(t));
  printf(" 00000000000 t\n");
  if (LT(t).a == to_ZZ("0"))
  {
    return 1;
  }
  else
  {
    printf("dame1\n");
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
  Div D3, D, null = {};
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
  // exit(1);
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
  if ((n & 1) == 1)
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

  cout << "p mod = " << a << "@" << p % 4 << p << endl;

  if (p % 4 == 3 || p % 8 == 5)
  {
    if (p % 4 == 3)
    {
      b = (p + 1) / 4;
      c = pow_mod(a, b, p);
      {
        // printf("good\n");
        cout << c << endl;
        return c;
      }
    }
    if (p % 8 == 5)
    {
      c = pow_mod(a, (p + 3) / 8, p);
      if ((c * c) % p != a)
      {
        printf("baka2\n");
        c = 2 * a * pow_mod(4 * a, (p - 5) / 8, p);
        if (c * c % P == a)
        {
          printf("good\n");
          return c;
        }
        // return to_ZZ("-1");
      }
      else if (c * c % p == a)
      {
        printf("good\n");
        cout << c << endl;
        return c;
      }
    }
  }
  if (p % 8 == 1 || p % 4 == 1)
  {
    c = tonelli_shanks(a, p);
    if (c * c % p != a)
    {
      printf("fail!\n");
      return to_ZZ("-1");
    }
    else
    {
      return c;
    }
  }
  //  return to_ZZ("0");
  printf("get back\n");
  exit(1);
}

// 曲線に代入した値を計算する
PO tr1e(ZZ f4, ZZ f3, ZZ f2, ZZ f1, ZZ f0, ZZ p)
{
  ZZ x, y, f, g;
  PO aa;

  while (1)
  {
    // while(x==0)
    {
      x = (ZZ)xor128() % p; // % p;
    }
    // y = xor128(); // % p;
    f = (pow_mod(x, to_ZZ("5"), p) + (f4 * pow_mod(x, to_ZZ("4"), p)) % p + (f3 * pow_mod(x, to_ZZ("3"), p)) % p + (f2 * pow_mod(x, to_ZZ("2"), p)) % p + f1 * x + f0) % p;
    y = root(f, p);
    g = (y * y) % p;
    cout << "f= " << f << " g= " << g << endl;
    cout << "x= " << x << " y= " << y << endl;
    if (f == g)
    {
      aa.x = x;
      aa.y = y;
      // exit(1);
      return aa;
    }
  }
  //  return -1;
}

OP diviser(PO o, PO m)
{
  ZZ t1[2][3] = {}, cc[2] = {};
  vec c1, c2;
  int i, j;
  ZZ a, b;
  OP f;

  t1[0][0] = o.x; // o.t[1].a;
  t1[1][0] = m.x; // t[0].a;
  t1[0][1] = 1;   // t[1].a;
  t1[1][1] = 1;   // t[0].a;
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

  cc[0] = inv(t1[0][0], P) % P;
  cout << cc[0] << endl;
  if (cc[0] > P)
    exit(1);
  ZZ z;
  z = t1[1][0] % P;
  cout << "z=" << z << endl;

  cc[0] = (cc[0] * z) % P;
  cout << "c0=" << cc[0] << endl;
  // exit(1);

  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < 3; j++)
      cout << t1[i][j] << endl;
    printf("\n");
  }
  printf("\n");

  for (j = 0; j < 3; j++)
  {
    t1[1][j] = (P + ((t1[0][j] * cc[0]) % P - t1[1][j])) % P;
    cout << "," << t1[1][j] << endl;
  }
  printf("\n");
  b = (t1[1][2] * inv(t1[1][1], P) % P);
  cout << "b=" << b << endl;
  for (i = 0; i < 3; i++)
    a = (t1[0][2] - (t1[0][1] * b) % P);
  if (a < 0)
    a += P;
  a = (inv(t1[0][0], P) * a) % P;
  cout << "a=" << a << endl;

  c1.x[0] = b;
  c1.x[1] = a;
  cout << c1.x[0] << " , " << c1.x[1] << endl;
  // exit(1);
  f = v2o(c1);
  printpoln(o2v(f));
  // exit(1);

  return f;
}

// ランダムな因子の生成
Div gendiv(OP f)
{

  PO a, b, e;
  OP d1 = {}, d2 = {}, c = {}, d = {}, vv1 = {}, vv2 = {}, uu1, v;
  //  ZZ  x, y, i, j, k,

  Div D = {};
  vec v1, v2, z1, z2, ff, vx;
  EX V;

  vx = o2v(f);
  /// srand(clock());
  v1.x[1] = to_ZZ("1");
  v2.x[1] = to_ZZ("1");
  //  while(1)
  {
    do
    {
      a = tr1e(vx.x[4], vx.x[3], vx.x[2], vx.x[1], vx.x[0], P); // cofficient of function
      b = tr1e(vx.x[4], vx.x[3], vx.x[2], vx.x[1], vx.x[0], P); // cofficient of function
    } while ((a.x == b.x) || (a.x == 0 || b.x == 0));
    cout << "P= " << P << endl;
    cout << "P-x= " << P - a.x << endl;
    cout << "x= " << a.x << endl;
    v1.x[0] = P - a.x;

    printpol(v1);
    printf("ppppppppppppp\n");
    c = v2o(v1);

    v2.x[0] = P - b.x;
    d = v2o(v2);
    d2 = diviser(a, b);
    printpoln(o2v(d2));
    // exit(1);

    // d2 = v2o(z1);
    d1 = omul(c, d);
    printpoln(o2v(d1));
    // exit(1);
    D.u = d1;
    D.v = d2;
    if (chkdiv(D, f) != -1)
    {
      printf("line\n");
      return D;
    }
    else
    {
      printf("just say buggy\n");
      exit(1);
    }
    //
  } // while (chkdiv(D, f) == -1);

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

Div tbl[1024] = {};
ZZ tmp[1024];
// 演算テーブルを作る（繰り返し2乗法）
void mktbl(Div D, OP f)
{
  int i;
  printf("begin\n");
  tbl[0] = D;
  if (chkdiv(D, f) == -1)
  {
    printf("what?\n");
    exit(1);
  }
  for (i = 0; i < 512; i++)
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

int flg = 0;
Div jac2(ZZ n, OP f, Div D)
{
  Div L;

  // D=gendiv(f);

  if (flg == 0 && n % 2 == 1)
  {
    L = D;
    n = (n >> 1);
    D = cdbl(D, f);
    flg = 1;
  }

  if (n % 2 == 1)
  {
    D = cadd(f, D.u, L.u, D.v, L.v);
    L = D;
    D = cdbl(D, f);
    n = (n >> 1);
  }
  else
  {
    D = cdbl(D, f);
    n = (n >> 1);
  }

return D;
}

// 因子のスカラー倍
Div jac(ZZ n, OP f)
{
  int i, j = 0, tmp[1024] = {};
  ZZ k;
  Div L = {}, D, G;

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
    // G = L;
    if (chkdiv(tbl[tmp[i]], f) == -1)
    {
      printf("before\n");
      // printpoln(o2v(tbl[tmp[i]].u));
      // printpoln(o2v(tbl[tmp[i]].v));
      exit(1);
    }

    L = cadd(f, tbl[tmp[i]].u, L.u, tbl[tmp[i]].v, L.v);
    if (chkdiv(L, f) == -1)
    {
      printf("dame2 %d\n", i);
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
      //printpoln(o2v(D.u));
      //printpoln(o2v(D.v));
      cout << "infinity devide! : " << n << endl;
      exit(1);
    }
  }
  printpoln(o2v(D.u));

  return L;
}

// 例
int main()
{

//  ZZ bb=to_ZZ("1419857");
//  ZZ f[K + 1] = {to_ZZ("1"), to_ZZ("0"), to_ZZ("0"), to_ZZ("0"), to_ZZ("0"), bb};
//  ZZ JC=to_ZZ("115792089237316195441280509386542044597668387405625706622058883245386161540656");

//  ZZ JC = to_ZZ("231584178474632390858684407602067108341863231303763314614616136159404179903562"); //omake
//  ZZ f[K + 1] = {to_ZZ("1"), to_ZZ("0"), to_ZZ("318258242717201709453901384328569236653"), to_ZZ("75380722035796344355219475510170298006"), to_ZZ("129416082603460579272847694630998099237"), to_ZZ("143864072772599444046778416709082679388")};

  // ZZ aa = to_ZZ("17");
  // ZZ JC=to_ZZ("28948022309328876598047865106864508015392396932332911323766668185591760466882");

  // ZZ aa=to_ZZ("2187");
  // ZZ JC=to_ZZ("115792089237316195413084848286254559191634403377276171765138133738613798321602");

  ZZ aa = to_ZZ("371293");                                                                                            // default
  ZZ JC = to_ZZ("8543948143683640329580084318401338115672828124663448275867130387651937373152534160174163969676194"); // 321
  // 2 · 4271974071841820164790042159200669057836414062331724137933565193825968686576267080087081984838097
  ZZ f[K + 1] = {to_ZZ("1"), to_ZZ("0"), to_ZZ("0"), to_ZZ("0"), aa, to_ZZ("0")};

  OP ff, k, uu1, uu2, vv1, vv2, s, o, d, c, m, d1, d2;
  EX V;
  Div D;
  oterm a;
  OP b = {};

  vec vx, xv, v11, v22;
  Div G0, G1, X;
  PO a11, b11;
  ZZ x, y, I;

  ff = setpol(f, K + 1);
  printpoln(o2v(ff));

  //　ランダムな因子をヤコビ多様体の位数倍して無限遠点になれば正しい
  srand(clock());
  // while(1)
  {
    X = gendiv(ff);
    if (chkdiv(X, ff) == -1)
    {
      printf("buggy may\n");
      exit(1);
    }
    else
    {
      printpoln(o2v(X.u));
      printpoln(o2v(X.v));
      printf("good\n");
    }
  }

  mktbl(X, ff);
  // while(I<P*P)
  {
    X = jac(JC + 1, ff);
    if (chkdiv(X, ff) == -1)
    {
      printf("bakayo\n");
      exit(1);
      // break;
    }
    I++;
  }

  return 0;
}