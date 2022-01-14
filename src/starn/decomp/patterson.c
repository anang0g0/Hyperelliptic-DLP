#include "common.c"

//パターソンアルゴリズムでバイナリGoppa符号を復号する
vec patterson(OP w, OP f)
{
  OP g1 = {0}, ll = {0}, op = {0};
  int i, j, k, l, c, n, count = 0, flg2 = 0;
  int flg, o1 = 0;
  OP tt = {0}, ff = {0}, h = {0};
  EX hh = {0}, opp = {0};
  vec v;
  oterm rr = {0};
  OP r2 = {0}, b2 = {
                   0};
  //unsigned short g[K+1]={2,2,12,1,2,8,4,13,5,10,8,2,15,10,7,3,5};

  tt.t[0].n = 1;
  tt.t[0].a = 1;
  o1 = 0;

  ff = inv(f, w);
  //  //printpol (o2v (ff));
  b2 = omod(omul(ff, f), w);
  if (odeg((b2)) > 0)
  {
    printf("逆元が計算できません。error\n");
    //wait ();
    exit(1);
  }
  printf("locater==========\n");
  //exit(1);
  r2 = oadd(ff, tt);
  printpol(o2v(r2));
  printf(" h+x==============\n");
  //wait();
  //  exit(1);
  g1 = osqrt(r2, w);
  printpol(o2v(g1));
  printf(" sqrt(h+x)=======\n");
  b2 = omod(omul(g1, g1), w);
  printpol(o2v(b2));
  printf(" g1*g1%w===========\n");
  if (deg(o2v(b2)) == 0)
  {
    printf("b2 is 0\n");
    exit(1);
  }
  if (LT2(b2).a != LT2(r2).a)
  {
    printpol(o2v(w));
    printf(" w============\n");
    printpol(o2v(r2));
    printf(" r2============\n");
    printpol(o2v(g1));
    printf(" g1============\n");
    printf(" g1は平方ではありません。error");
    //wait ();
    exit(1);
  }
  printpol(o2v(w));
  printf(" ========goppa\n");
  printpol(o2v(g1));
  printf(" g1!=========\n");
  if (LT(g1).n == 0 && LT(g1).a == 0)
  {
    //printpol (o2v (w));
    printf(" badkey=========\n");
    //printpol (o2v (g1));
    printf(" g1============\n");
    printf("平方が０になりました。 error\n");
    //wait ();
    exit(1);
  }
  h = gcd(f, w);
  if (odeg((h)) > 0)
  {
    printsage(o2v(w));
    printf(" =goppa\n");
    printsage(o2v(f));
    printf(" =syn\n");
    printpol(o2v(h));
    printf(" =gcd\n");

    //printf("f-ben_or=%d\n",ben_or(f));
    //printf("w-ben_or=%d\n",ben_or(w));
    printf(" s,wは互いに素じゃありません。\n");
    flg2 = 1;
    //wait ();
    exit(1);
    //goto label;
  }

  //exit(1);
  OP ppo = {0};

  //hh = xgcd(w, g1);
  ppo = zgcd(w, g1,T);
  printpol(o2v(ppo));
  printf(" =========ppo\n");
  printpol(o2v(ppo));
  printf(" =========ppo\n");
  //if(LT(ppo).a!=LT(ppo).a)
  //exit(1);
  flg = 0;
  printpol(o2v(ppo));
  printf(" =========ppo\n");
  ff = omod(omul(ppo, g1), w);
  printpol(o2v(ppo));
  printf(" alpha!=========\n");
  printpol(o2v(ff));
  printf(" beta!=========\n");
  printpol(o2v(w));
  printf(" goppa=========\n");
  printpol(o2v(f));
  printf(" syn=========\n");

  printpol(o2v(ppo));
  printf(" alpha!=========\n");
  //hh = xgcd(w, g1);
  ppo = zgcd(w, g1,T);
  //  h=zgcd(w,g1);
  ff = omod(omul(ppo, g1), w);
  ll = oadd(omul(ff, ff), omul(tt, omul(ppo, ppo)));
  v = chen(ll);
  if (v.x[K - 1] > 0)
  {
    //wait ();
    return v;
  }

  ppo = zgcd(w, g1,T);
  ff = omod(omul(ppo, g1), w);
  if (odeg((ff)) != K / 2)
  {

    printf("\nbefore h.d\n");
    ff = omod(omul(ppo, g1), w);
    flg = 1;
    printpol(o2v(ff));
    printf(" ==========beta!\n");
    printpol(o2v(ppo));
    printf(" alpha!=========\n");
    ll = oadd(omul(ff, ff), omul(tt, omul(ppo, ppo)));
    v = chen(ll);
    if (v.x[K - 1] > 0)
    {
      wait();
      return v;
    }
  }

  if (odeg((ff)) == 1)
  {
    ll = oadd(omul(ff, ff), omul(tt, omul(ppo, ppo))); //ff;
    printf("deg==1\n");
    exit(1);
  }

  printf("あっ、でる・・・！\n");
  count = 0;
  printpol(o2v(ll));
  printf(" ll=================\n");
  printpol(o2v(f));
  printf(" syn=================\n");
  v = chen(ll);

  if (v.x[K - 1] > 0)
  {
    return v;
  }

  printf("flg2=%d\n", flg2);
  //exit(1);
  //ero2(v);
  return v;
}



//512bitの秘密鍵を暗号化
void encrypt(char buf[], unsigned char sk[64])
{
  const uint8_t *hash = {0};
  sha3_context c = {0};
  int image_size = 512, i;
  FILE *fp;
  //  unsigned short dd=0;

  printf("plain text=");
  for (i = 0; i < 64; i++)
    printf("%u,", sk[i]);
  printf("\n");
  //  puts(buf);
  //printf("\n");
  //exit(1);

  //scanf("%s",buf);
  sha3_Init256(&c);
  sha3_Update(&c, (char *)buf, strlen(buf));
  hash = sha3_Finalize(&c);

  //j=0;

  for (i = 0; i < 64; i++)
  {
    printf("%d", hash[i]);
    //char s[3];
    //byte_to_hex(hash[i],s);

    sk[i] ^= hash[i];
  }
  printf("\nencrypt sk=");
  for (i = 0; i < 64; i++)
    printf("%d,", sk[i]);
  printf("\n");

  fp = fopen("enc.sk", "wb");
  fwrite(sy, 2, K, fp);
  fwrite(sk, 1, 64, fp);
  fclose(fp);
}
OP synd(unsigned short zz[])
{
  unsigned short syn[K] = {0};
  unsigned short s = 0;
  int i, j, t1;
  OP f = {0};

  printf("in synd2\n");

  for (i = 0; i < K; i++)
  {
    syn[i] = 0;
    s = 0;
    //#pragma omp parallel num_threads(16)
    for (j = 0; j < N; j++)
    {
      s ^= gf[mlt(fg[zz[j]], fg[mat[j][i]])];
    }
    syn[i] = s;
    //printf ("syn%d,", syn[i]);
  }
  //printf ("\n");

  f = setpol(syn, K);
  printpol(o2v(f));
  printf(" syn=============\n");
  //  exit(1);

  return f;
}


void decrypt(OP w)
{
  FILE *fp;
  int i, j;
  unsigned char sk[64] = {0}, err[N] = {
                                  0};
  unsigned short buf[K] = {0}, tmp[K] = {
                                   0};
  OP f = {0}, r = {0};
  vec v = {0};
  const uint8_t *hash = {0};
  sha3_context c = {0};
  int image_size = 512;

  j = 0;
  fp = fopen("enc.sk", "rb");

  fread(tmp, 2, K, fp);
  fread(sk, 1, 64, fp);
  fclose(fp);

  for (i = 0; i < K; i++)
    buf[i] = tmp[K - i - 1];

  printf("in decrypt\n");
  f = setpol(buf, K);
  v = patterson(w, f);

  // elo(r);
  //exit(1);
  //v=o2v(r);

  j = 0;
  if (v.x[1] > 0 && v.x[0] == 0)
  {
    err[0] = 1;
    j++;
  }

  printf("j=%d\n", j);
  printf("after j\n");
  for (i = j; i < 2 * T; i++)
  {
    if (v.x[i] > 0 && v.x[i] < N)
    {
      err[v.x[i]] = 1;
    }
  }

  char buf0[8192] = {0}, buf1[10] = {
                             0};

  //#pragma omp parallel for
  for (i = 0; i < N; i++)
  {
    snprintf(buf1, 10, "%d", err[i]);
    strcat(buf0, buf1);
  }
  //puts (buf0);
  printf("vector=%d\n", strlen(buf0));
  //exit(1);
  printf("cipher sk2=");
  for (i = 0; i < 64; i++)
    printf("%u,", sk[i]);
  printf("\n");

  sha3_Init256(&c);
  sha3_Update(&c, (char *)buf0, strlen(buf0));
  hash = sha3_Finalize(&c);

  j = 0;
  printf("hash=");
  for (i = 0; i < 64; i++)
  {
    printf("%d", hash[i]);
    //char s[3];
    //byte_to_hex(hash[i],s);

    sk[i] ^= hash[i];
  }
  printf("\ndecript sk=");
  for (i = 0; i < 64; i++)
    printf("%u,", sk[i]);
  printf("\n");
  //  exit(1);

  return;
}


//64バイト秘密鍵の暗号化と復号のテスト
void test(OP w, unsigned short zz[])
{
  int i;
  vec v = {0};
  const uint8_t *hash;
  sha3_context c;
  //int image_size=512;
  OP f = {0};
  FILE *fp;

  fp = fopen("aes.key", "rb");
  /*
     static char base64[] = {
     'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H',
     'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P',
     'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X',
     'Y', 'Z', 'a', 'b', 'c', 'd', 'e', 'f',
     'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n',
     'o', 'p', 'q', 'r', 's', 't', 'u', 'v',
     'w', 'x', 'y', 'z', '0', '1', '2', '3',
     '4', '5', '6', '7', '8', '9', '+', '/',
     };
   */
  char buf[8192] = {0}, buf1[10] = {
                            0};
  unsigned char sk[64] = {0};
  // unsigned short s[K]={0};
  //fread(sk,1,32,fp);
  for (i = 0; i < 64; i++)
    sk[i] = i + 1;

  for (i = 0; i < N; i++)
  {
    snprintf(buf1, 10, "%d", zz[i]);
    strcat(buf, buf1);
  }
  //puts (buf);
  printf("vector=%u\n", strlen(buf));
  //exit(1);

  printf("sk0=");
  for (i = 0; i < 64; i++)
    printf("%u,", sk[i]);
  printf("\n");
  //exit(1);

  f = synd(zz);
  v = o2v(f);
  //printf("v=");
  for (i = 0; i < K; i++)
  {
    sy[i] = v.x[i];
    printf("%d,", sy[i]);
  }
  printf("\n");

  encrypt(buf, sk);
  decrypt(w);

  sha3_Init256(&c);
  sha3_Update(&c, (char *)buf, strlen(buf));
  hash = sha3_Finalize(&c);
}


