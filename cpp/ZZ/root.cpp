#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "pstruct.h"
#include <NTL/ZZ.h>

NTL_CLIENT

int bit(ZZ b, int i)
{
  int k = 1;

  if (((b & (1 << i)) >> i)%2 == 1)
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
    x = (ZZ)rand() % p;
    // y = rand() % p;
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



int main()
{
unsigned i;
ZZ P=to_ZZ("100000000000000003");
//ZZ xx=72976196454691585;
//printf("%llu\n",xx*xx%P);
ZZ xx;
PO x;

ZZ b = (P + 1) / 4;
ZZ c = pow_mod(to_ZZ("10"), b, P);
cout << c << endl;

//cout << root(to_ZZ("10"),P) << endl;
srand(clock());
//for(i=0;i<64;i++)
//printf("%d,",bit(P,i));
//printf("\n");
//exit(1);

//for(xx=0;xx<10000;xx++){
 x= tr1e(to_ZZ("1"),to_ZZ("2"),to_ZZ("3"),to_ZZ("4"),to_ZZ("5"),P);
 cout << x.x << ", " << x.y << endl;
//}

    return 0;
  }