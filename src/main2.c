
int main()
{
  unsigned int i, count = 0;
  unsigned long long  aaa[O] = {0};
  

  unsigned long long  f[K + 1] = {1, 7, 6, 2, 8, 2};
/*
  unsigned long long  u2[K + 1] = {0, 0, 0, 1, 21, 16};
  unsigned long long  u1[K + 1] = {0, 0, 0, 1, 19, 20};
  unsigned long long  v2[K + 1] = {0, 0, 0, 0, 21, 21};
  unsigned long long  v1[K + 1] = {0, 0, 0, 0, 12, 8};
*/
// unsigned long long  f[K+1]={1 ,0, 2,  30,  5,  1};
/*
  //unsigned long long  f[K+1]= {1, 1597 , 1041 ,5503 , 6101 , 1887 };
//f1 = x + 28555025517563816 and f2 = x + 74658844563359755 ;
unsigned long long  u2[K+1]={0,0,0,1,1571353025997967 , 12198441063534328};
unsigned long long  v2[K+1]={0,0,0,0,32227723250469108 , 68133247565452990};
unsigned long long  u1[K+1]={0,0,0,1, 70887725815800572 , 94321182398888258};
unsigned long long  v1[K+1]={0,0,0,0, 42016761890161508 , 3182371156137467 };
*/
unsigned long long u2[K+1]={0,0,0,1,26,20};
unsigned long long v2[K+1]={0,0,0,0,29,26};
unsigned long long u1[K+1]={0,0,0,1,9,27};
unsigned long long v1[K+1]={0,0,0,0,29,16};

//  unsigned long long  f[K + 1] = {1, 0, 3, 7, 1, 2};
/*
  unsigned long long  u2[K + 1] = {0, 0, 0, 1, 7, 10};
  unsigned long long  u1[K + 1] = {0, 0, 0, 1, 0, 10};
  unsigned long long  v2[K + 1] = {0, 0, 0, 0, 1, 9};
  unsigned long long  v1[K + 1] = {0, 0, 0, 0, 7, 9};
*/
  OP ff, k, uu1, uu2, vv1, vv2, s, l, u3, v3, u, ll, t, m, o, d, c;
  unsigned long long  tst1[K + 1] = {0, 0, 0, 0, 8, 7};
  unsigned long long  tst2[K + 1] = {0, 0, 0, 0, 0, 10};
  unsigned tmp[2][3]={0};
  EX V;
  oterm a;
  OP b={0};
//unsigned long long  q1 = to_unsigned long long ("1208925819614629174708801");
int a1 = 1331;
//J1 =to_unsigned long long ("1461501637326815988079848163961117521046955445901");
//e y2 = x5+a, a ∈ Fp

//unsigned long long  q2 = to_unsigned long long ("1208925819614629174709941");
int a2 = 2;
//J2 = to_unsigned long long ("1461501637331762771847359428275278989652932675771");

vec vx={0},xv={0};
Div G0,G1,X;
  ff = setpol(f, K + 1);


  uu1 = setpol(u1, K + 1);
  uu2 = setpol(u2, K + 1);
  vv1 = setpol(v1, K + 1);
  vv2 = setpol(v2, K + 1);
  o=setpol(tst1,K+1);
  m=setpol(tst2,K+1);

//V=xgcd(uu1,uu2);
G0=cadd(ff,uu1,uu2,vv1,vv2);
printf("%d\n",chkdiv(G0,ff));
//exit(1);

G1.u=uu1;
G1.v=vv1;
X.u=uu2;
X.v=vv2;
if(chkdiv(G1,ff)==-1 || chkdiv(X,ff)==-1){
    printf("erro!\n");
exit(1);
}
G0=cadd(ff,G1.u,X.u,G1.v,X.v);
if(chkdiv(G0,ff)==-1)
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
srand(clock());
/*
G0=cadd(ff,uu1,uu2,vv1,vv2);
if(chkdiv(G0,ff)==-1){
    printf("bug\n");
    exit(1);
}
*/
/*
for(rr=0;rr<1000;rr++){
  printf("%llu %llu\n",root(rr,P),rr);
}
exit(1);
vx=o2v(ff);
xx = tr1e(vx.x[4], vx.x[3], vx.x[2], vx.x[1], vx.x[0], P); // cofficient of function
printf("%llu, %llu\n",xx.x, xx.y);
exit(1);
*/
/*
for(i=0;i<1000;i++){
  x=i;
  
}
*/

PO xx;
unsigned long long rr=0;
count=0;
int xount=0;

G1=gendiv(ff);
X=gendiv(ff);
uu1=G1.u;
vv1=G1.v;
srand(clock());
while(1){
G1=gendiv(ff);
X=gendiv(ff);
G0=cadd(ff,G1.u,X.u,G1.v,X.v);

if(chkdiv(G,ff)==-1)
{
printf("baka\n");
 count++;
  //exit(1);
  V=xgcd(X.u,G1.u);
  if(LT(V.d).n>0){
  printf("gcd!\n");

  //exit(1);
  }
  }else if(oequ(G1.u,uu1)==0){
    printf("order #J= %d\n",xount);
    exit(1);
  }else if(chkdiv(G1,ff)!=-1){
  printf("ウホッ！いい因子。\n");
  xount++;
  //exit(1);
}else{
    printpoln(o2v(G0.u));
    printpoln(o2v(G0.v));
    printpoln(o2v(G1.u));
    printpoln(o2v(G1.v));
    printpoln(o2v(X.u));
    printpoln(o2v(X.v));
  printf("why?\n");
  //count++;
  exit(1);
  }
if(count>100)
break;
printf("%u xount=%u\n",count,xount);
}
printf("%u %u\n",count,xount);
exit(1);

/*
G0=cdbl(X,ff);
if(chkdiv(G0,ff)==-1)
{
  printf("naze?\n");
}else{
  printf("イイっ！この因子すげえいいっ！\n");
}
*/
//exit(1);
/*
G0=g2add(ff,uu1,uu2,vv1,vv2);
if(chkdiv(G0,ff)==-1){
  printf("nani?\n");
}else{
printf("だろうな\n");
}
exit(1);
*/
//b=odiv(o,m);
//exit(1);
/*
X=g2add(ff,uu1,uu2,vv1,vv2);
printpol(o2v(X.u));
printf("  Xu\n");
printpol(o2v(X.v));
printf("  Xv\n");
//printf("%llu\n",chkdiv(X,ff));
exit(1);
*/


/*
EX F;
Div U={0},U0={0};
while(1)
{
U=gendiv(ff);
U0=gendiv(ff);
printf("1111111111111111111111111111\n");
//F=xgcd(U.u,U0.u);
if(chkdiv(U,ff)==-1 || chkdiv(U0,ff)==-1)
{
  printf("ee?!\n");
  exit(1);
}
X=cadd(ff,U.u,U0.u,U.v,U0.v);
if(chkdiv(X,ff)==-1)
{
  printpoln(o2v(U.u));
  printpoln(o2v(U0.u));
  printpoln(o2v(U.v));
  printpoln(o2v(U0.v));
  printf("baka^^\n");
  exit(1);
}else{
  printf("happy!\n");
  //exit(1);
}
}
printf("isdiv=%llu\n",chkdiv(X,ff));
printpol(o2v(F.u));
printf("  Fuwwwwwwww\n");
printpol(o2v(F.v));
printf("  Fvwwwwwwww\n");
printpol(o2v(F.d));
printf("  Fdwwwwwwww\n");
printpol(o2v(F.h));
printf("  Fhwwwwwwww\n");
printpol(o2v(X.u));
printf("  Xuwwwwwwww\n");
printpol(o2v(X.v));
printf("  Xvwwwwwwww\n");
printpol(o2v(U.u));
printf("  Uu wwwwwwww\n");
printpol(o2v(U0.u));
printf("  U0u vwwwwwwww\n");
printpol(o2v(U.v));
printf("  Uv wwwwwwww\n");
printpol(o2v(U0.v));
printf("  U0v vwwwwwwww\n");
*/
//exit(1);
/*
//printf("isdiv=%llu\n",chkdiv(U,ff));
//printf("isdiv=%llu\n",chkdiv(U0,ff));
X=g2add(ff,U.u,U0.u,U.v,U0.v);
printpol(o2v(U.u));
printf("  isUwwwwwwww\n");
printpol(o2v(U0.u));
printf("  isU0wwwwwwww\n");
  printf("isdiv=%llu\n",chkdiv(X,ff));
  exit(1);
  */
OP oi,rem;
Div D1,D2;
/*
while(1){
  D1=gendiv(ff);
  D2=gendiv(ff);
  printpol(o2v(D1.u));
  printf(" ===D1.u\n");
  printpol(o2v(D1.v));
  printf(" ===D1.v\n");
  printpol(o2v(D2.u));
  printf(" ===D2.u\n");
  printpol(o2v(D2.v));
  printf(" ===D2.v\n");
  oi=qinv(D2.u,D1.u);
  printpol(o2v(oi));
  printf(" ===Doi\n");
  rem=omod(omul(oi,D2.u),D1.u);
  printpol(o2v(rem));
  printf(" ==D rem\n");
  exit(1);
  printf("isdiv=%llu D\n",chkdiv(D1,ff));
  printf("isdiv=%llu D\n",chkdiv(D2,ff));
  U = g2add(ff, D1.u, D2.u, D1.v, D2.v);
  printf("isdiv=%llu D\n",chkdiv(U,ff));
  //exit(1);
  
  if(chkdiv(U,ff)==-1){
    printf("buggy\n");
  exit(1);
  }
}
  //U = g2add(ff, uu1, uu2, vv1, vv2);
  //printf("isdiv=%llu\n",chkdiv(U,ff));
//    printf("isdiv=%llu\n",chkdiv(D2,ff));
//  exit(1);
  //V=xgcd(uu1,uu2,2);
  printpol(o2v(V.u));
  printf(" =====u3\n");
  printpol(o2v(V.v));
  printf(" =====v3\n");
*/

D1.u=uu1;
D1.v=vv1;
D2.u=uu2;
D2.v=vv2;

printf("isdiv=%llu\n",chkdiv(D1,ff));
printf("isdiv=%llu\n",chkdiv(D2,ff));
//exit(1);

D1=gendiv(ff);
D2=gendiv(ff);

printf("isdiv=%llu\n",chkdiv(D1,ff));
printf("isdiv=%llu\n",chkdiv(D2,ff));
//exit(1);
//U=g2add(ff,D1.u,D2.u,D1.v,D2.v);

  o=oadd(vv1,vv2);
  printpol(o2v(o));
  printf("\n");

  Div W;
  // below undercondtruction
  k = odiv(oadd(ff, minus(omul(D1.v, D1.v))), D1.u);
  s = omod(odiv(k, scr(2, D1.v)), D1.u);
  l = omul(s, uu1);
  u3 = omod(osub(omul(s, s), (osub(scr(2, omul(D1.v, s)), (k)))), D1.u);
  v3 = omod(minus(oadd(l, D1.v)), u3);
  W.u=u3;
  W.v=v3;
  printpol(o2v(u3));
  printf("======du3\n");
  printpol(o2v(v3));
  printf("======dv3\n");
  printf("W's isdiv=%llu\n",chkdiv(W,ff));
//  exit(1);


  return 0;
}