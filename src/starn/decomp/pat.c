//patx
#include "patterson.c"
//#include "8192.h"
//#include "4096.h"



//言わずもがな
int main(void)
{
  unsigned short z1[N] = {0}; //{1,0,1,1,1,0,0,0,0,0,1,1,1,0,0,1};
  OP f = {0}, r = {0}, w = {0};
  vec v;
  unsigned short ss[K] = {0};

srand(clock());
if (K > N){
  printf("configuration error! K is too big K\n");
  exit(1);
}

//kabatiansky example
unsigned short s[K+1]={0,15,1,9,13,1,14};
//Berlekamp-Massey法（UnderConstruction）
//bms(s);
//exit(1);

int j=0;

  //chu();


  //公開鍵を生成する
 w = pubkeygen();
 
while(1){

//エラーベクトルを生成する
  memset(z1, 0, sizeof(z1));
  mkerr(z1, T * 2);
  //exit(1);

  //encryotion
  test (w, z1);

  //シンドロームを計算する
  v=sin2(z1);
  printpol(o2v(f));
  printf(" ==syndrome\n");
  f=dec(v.x);
  //復号化の本体
  v=patterson(w, f);
  //エラー表示
  ero(v);

  break;
}

  return 0;
}
