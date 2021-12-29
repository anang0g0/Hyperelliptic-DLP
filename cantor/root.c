#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>


int bit(unsigned long long b, int i)
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

//jj=aa^bb mod oo
unsigned long long  eXp(unsigned long long  aa, unsigned long long  bb, unsigned long long  oo)
{
  unsigned long long  ii, jj, kk[8192];
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

unsigned long long root(unsigned long long a, unsigned long long p)
{
  unsigned long long  c, b;

  printf("p mod = %llu == %llu , %llu\n",a, p % 4, p );
  if (p % 4 == 3 || p % 8 == 5)
  {
    if (p % 4 == 3)
    {
      b = (p + 1) / 4;
      c = eXp(a, b, p);
      {
        //printf("good\n");
        printf("good c= %llu\n", c);
        return c;
      }
    }
    if (p % 8 == 5)
    {
      c = eXp(a, (p + 3) / 8, p);
      if ((c * c) % p != a)
        printf("baka2\n");
        return -1;
    }
    if (c * c % p == a)
    {
      printf("good\n");
      printf("%llu\n", c);
      return c;
    }
  }
  if (p % 8 == 5)
  {
    c = 2 * a * eXp(4 * a, (p - 5) / 8, p);
    if (c * c % p != a)
    {
      printf("dangerous\n");
      return -1;
    }
    if (c * c % p == a)
    {
      printf("good\n");
      printf("%llu\n", c);
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

int main()
{
unsigned i;
unsigned long long P=100000000000000003;
unsigned long long xx=72976196454691585;
//printf("%llu\n",xx*xx%P);

unsigned long long b = (P + 1) / 4;
unsigned long long c = pow_mod(10, b, P);
printf("%llu\n",c);

//printf("%llu\n",root(10,P));
//exit(1);

for(i=0;i<64;i++)
printf("%d,",bit(P,i));
printf("\n");
exit(1);

for(xx=0;xx<10000;xx++){
printf("%llu %llu\n", root(xx,P) , xx);
}


    return 0;
  }