#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <iostream>
#include <NTL/ZZ.h>


NTL_CLIENT

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


//jj=aa^bb mod oo
ZZ exp(ZZ aa, ZZ bb, ZZ oo)
{
  ZZ ii, jj, kk[8192];
  int j, c[8192], count = 0, i;
  ii = oo;
  j = 0;
  jj = 0;
  //  kk[4096]; //prime is 4096 bit table
  //  c[8192]  //mod is 8192 bit table
  count = 0;

  for (i = 0; i < 522; i++)
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

  //ex.1000=2**3+2**5+2**6+2**7+2**8+2**9 makes a array c=[3,5,6,7,8,9]
  for (i = 0; i < j + 1; i++)
  {
    if (bit(bb, i) != 0)
    { // testbit(bb,i)
      c[count] = i;
      count = count + 1;
    }
  }
  //    cout << bb << }l;
  //    cout << count << "\n";
  //exit(1);
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


ZZ pow_mod(ZZ x, ZZ n, ZZ p)
{
  if (n == 0)
    return to_ZZ("1");
  if ((n & 1)==1)
    return (pow_mod(x, n - 1, p) * x) % p;
  x = pow_mod(x, n / 2, p);
  return (ZZ)((x * x) % p);
}


int main(){
ZZ K,L;



}