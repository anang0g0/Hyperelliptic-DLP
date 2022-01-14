#include "decode.c"
//Xeuc

//言わずもがな
int main(void)
{
  unsigned short z1[N] = {0}; //{1,0,1,1,1,0,0,0,0,0,1,1,1,0,0,1};
  OP f = {0}, r = {0}, w = {0};
  vec v = {0};
  unsigned short ss[K] = {0}, zz[N] = {0};

  if (K > N)
  {
    printf("configuration error! K is too big K\n");
    exit(1);
  }

  //kabatiansky example
  unsigned short s[K + 1] = {0, 15, 1, 9, 13, 1, 14};
  //Berlekamp-Massey法（UnderConstruction）
  //bms(s);
  //exit(1);

  //公開鍵を生成する(H'=SHP)
  w = pubkeygen();

  //while(1){

  //エラーベクトルを生成する
  memset(z1, 0, sizeof(z1));
  mkerr(z1, T);
  //exit(1);

  //encryotion test
  //test (w, z1);
  //exit(1);

  //  mkd(w); ??

  memset(zz, 0, sizeof(zz));
  memset(ss, 0, sizeof(ss));

  int j = 0, count = 0;
  //decode開始
  int k = 0;
  while (1)
  {

    memset(z1, 0, sizeof(z1));
    mkerr(z1, T);

    for (int i = 0; i < N; i++)
    {
      if (z1[i] > 0)
        printf("la=%d %d\n", i, z1[i]);
    }
    //暗号化(v=eH')
    v = sin2(z1);
    //f=v2o(v);
    //復号(S^-1)
    f = dec(v.x);
    r = decode(w, f);

    //m=m'P^-1
    //count = elo2(r);
    elo2(r);
    //exit(1);

    for (int i = 0; i < N; i++)
    {
      //検算
      if (z1[i] > 0)
        printf("error position= %d %d\n", i, z1[i]);
    }
    if (count < 0)
    {
      printf("baka-@\n");
      exit(1);
    }
    j++;
    printf("err=%dっ！！\n", count);
    for (int i = 0; i < N; i++)
      printf("%d,", z1[i]);
    printf("\n");
    exit(1);

    if (j == 10000)
      exit(1);
  }

  return 0;
}

