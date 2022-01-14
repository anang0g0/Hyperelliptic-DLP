#include "oplib.c"




void chu(void)
{
  int i, j, k, l, ii = 0;
  OP w = {0};
  unsigned short tr[N] = {0};
  unsigned short ta[N] = {0};

  j=0;
aa:

  //printf("\n");

  //既約性判定のためのBen-Orアルゴリズム。拡大体にも対応している。デフォルトでGF(8192)
  //既約多項式しか使わない。
  l = -1;
  ii = 0;


  while (l == -1)
  {
    w = mkpol();
    l = ben_or(w);
    printf("irr=%d\n", l);
    if (ii > 300)
    {
      printf("too many tryal\n");
      goto aa;
      //exit(1);
    }
    printf("ben=%d\n",ii);
    ii++;
    //
  } 
  
  //w = mkpol();

  //多項式の値が0でないことを確認
  for (i = 0; i < N; i++)
  {
    ta[i] = trace(w, i);
    if (ta[i] == 0)
    {
      printf("trace 0 @ %d\n", i);
      //fail = 1;
      exit(1);
    }
  }
  for (i = 0; i < N; i++)
  {
    tr[i] = oinv(ta[i]);
    //printf("%d,", tr[i]);
  }

  printpol(o2v(w));
  printf(" =irreducible\n");
  printsage(o2v(w));
  printf("\n");
  j++;
  if(j>100)
  exit(1);
    goto aa;

}


int main(void)
{
  int i, j = 0, count = 0;
  unsigned short zz[N] = {0}, z1[N] = {0};
  OP f = {0}, r = {0}, w = {0};
  vec v, x = {0};
  MTX R = {0}, O = {0}, Q = {0};
  unsigned short s[K + 1] = {0};
  arrayul ary={0};
  //srand(clock());
/*
//while(1){
chash(ary.d);
printf("%d\n",ary.t[1]);
//}
exit(1);
*/
f=mkpol();
ben_or(f);
//chu();

  Pgen();
  do
  {
    memset(S.x, 0, sizeof(S.x));
    memset(inv_S.x, 0, sizeof(inv_S.x));
    for (i = 0; i < K * E; i++)
    {
      for (j = 0; j < K * E; j++)
        S.x[i][j] = xor128() % 2;
    }
  } while (is_reg(S, inv_S.x) == -1);
  //exit(1);
  
  if (K > N)
    printf("configuration error! K is bigger than N\n");

  // （謎）
  memset(mat, 0, sizeof(mat));

  // 公開鍵を生成する(Berlekamp-Massey用)
  R = pk_gen();
  while(1){
  // エラーベクトルの初期化
  memset(zz, 0, sizeof(zz));
  //重みTのエラーベクトルを生成する
  mkerr(zz, T);
  // 暗号文の生成(s=eH)
  x = sin2(zz, R);
  // 復号化１(m'=sS^{-1})
  r = dec(x.x);
  v = o2v(r);
  for (i = 0; i < K; i++)
    s[i + 1] = v.x[i];

  // Berlekamp-Massey Algorithm
  f = bma(s, K);
  x = chen(f);
  // 平文の表示(m=m'P^{-1})
  ero2(x);
  //wait();
  break;
  }
exit(1);
  
  //公開鍵を生成する(G=SHP : Niederreiter , Patterson共用)
  O = pubkeygen(w);

  //decode開始
  while (1)
  {
    memset(z1, 0, sizeof(z1));
    //T重みエラーベクトルの生成
    mkerr(z1, T);
    //公開鍵を使ったシンドロームの計算(s=eG)
    v = sin2(z1, O);
    //シンドロームの復号(s'=sS^{-1})
    f = dec(v.x);
    //復号(m'=D(s'))
    r = decode(O.f, f);
    //エラー（平文）の表示(m=m'P^{-1})
    count = elo2(r);
    for (int i = 0; i < N; i++)
    {
      //検算
      if (z1[i] > 0)
        printf("error position= %d %d\n", i, z1[i]);
    }
    if (count < T)
    {
      printpol(o2v(w));
      printf(" == Goppa polynomial\n");
      exit(1);
    }
    j++;
    if (count == T)
      printf("err=%dっ！！\n", count);

    //if (j == 10000)
    break;
  }
  wait();

  //  O = pubkeygen(w);
  while (1)
  {

    //エラーベクトルを生成する
    memset(z1, 0, sizeof(z1));
    mkerr(z1, T * 2);
    //exit(1);

    //encryotion
    test(O.f, z1);

    //シンドロームを計算する
    x = sin2(z1, O);
    f = dec(x.x);
        //f=v2o(x);
        printpol(o2v(f));
    printf(" ==syndrome\n");

    //復号化の本体
    v = patterson(O.f, f);
    //エラー表示
    ero(v);

    break;
  }
  wait();

 do
  {
    memset(Q.x, 0, sizeof(Q.x));
    memset(O.x, 0, sizeof(O.x));
    for (i = 0; i < (K/2+1) * E; i++)
    {
      for (j = 0; j < (K/2+1) * E; j++)
        S.x[i][j] = xor128() % 2;
    }
  } while (mkRS(S, inv_S.x) == -1);
// exit(1);
  O = mk_pub();
  memset(zz, 0, sizeof(zz));
  mkerr(zz, T);
  r = sendrier2(zz, K, O);
  x = o2v(r);
  for (i = 0; i < K; i++)
    s[i + 1] = x.x[i];
  //for (i = 0; i < K; i++)
  //    printf("%d,", s[i]);
  //printf("\n");
  f = bma(s, K);
  x=chen(f);
  ero2(x);
  for (i = 0; i < N; i++)
    if (zz[i] > 0)
      printf("%d,", i);
  printf("\n");

  if (odeg(f) < T)
  {
    printpol(o2v(r));
    printf("==goppa\n");
    for (i = 0; i < N; i++)
      printf("%d,", zz[i]);
    printf("\n");
    exit(1);
  }
  wait();
  

  
while(1){
  memset(O.x,0,sizeof(O.x));
  O.row=64;
  O.col=64;
  for(i=0;i<64;i++){
    for(j=0;j<64;j++)
    O.x[i][j]=rand()%2;
  }
  if(binv(O,Q.x,64)==0)
  break;
}
Q.col=64;
Q.row=64;
printf("aaaaaaa\n");
//exit(1);
mmul(O,Q);
//matinv(K);

  return 0;
}
