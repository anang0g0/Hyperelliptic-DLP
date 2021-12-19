#include <NTL/ZZ.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define N 256

NTL_CLIENT

ZZ le_u[N];
ZZ le_v[N];
ZZ uu[N];
ZZ vv[N];
ZZ U[N];
ZZ V[N];

ZZ ue1,ue0,ve0,ve1,u31,u30,v31,v30;

/* definition of a curve */
void HEC(){
  //find by Gaudry
  //fg=x^5+2682810822839355644900736*x^3+226591355295993102902116*x^2+2547674715952929717899918*x+4797309959708489673059350
  ZZ P;
  P = to_ZZ("5000000000000000008503491");
  ZZ FF[5];// ={0};
  char* t[5]={"0", "2682810822839355644900736", "226591355295993102902116", "2547674715952929717899918", "4797309959708489673059350"};
  ZZ Jga;
  Jga = to_ZZ("24999999999994130438600999402209463966197516075699");
  char* t2[2]={"1", "-1713538969626908355896596"}; 
  char* t3[2]={"0", "138905579055173741542118"};
  ZZ ug0[2];
  ZZ vg0[2];
  char* t4[2]={"1738366273424896804842766", "3184841659043138633535652"};
  ZZ ug1[2];
char* t5[2]={"2931056213155642836850986", "402980510097415333052905"};
  ZZ vg1[2];
for(int i=0;i<5;i++)
  FF[i] = to_ZZ(t[i]);
for(int i=0;i<2;i++){
  ug0[i]=to_ZZ(t2[i]);
  vg0[i]=to_ZZ(t3[i]);
  ug1[i]=to_ZZ(t4[i]);
  vg1[i]=to_ZZ(t5[i]);
}
}

// invert of integer
ZZ inv(ZZ a,ZZ n){
ZZ  d;
ZZ q,t,r,x,s,gcd;

 x = to_ZZ("0");
 s = to_ZZ("1");

d = n;
  while (a != 0){
    q = d / a;
    r = d % a;
    d = a;
    a = r;
    t = x - q * s;
    x = s;
    s = t;
  }
  gcd = d;

  return ((x + n) % (n / d));
}

ZZ g2add(ZZ u11,ZZ u10,ZZ v11,ZZ v10,ZZ u21,ZZ u20,ZZ v21,ZZ v20,ZZ f4,ZZ f3,ZZ f2,ZZ f1,ZZ f0,ZZ p){
ZZ a0, a1, ue1, ul, ul0, uu0, uu1, u30, u31, b0, b1, e0, e1, c0, c1, d0, d1, ss2, ss3, s1, s2, s3, v30, v31;
ZZ iv, l2, l3, l1;
  a0 = u11 * u10; a1 = u21 * u20; d0 = u10 - u20; d1 = u11 - u21;
  b0 = u11 * u11; b1 = u21 * u21; c0 = v20 - v10; c1 = v21 - v11;
  e0 = -u20 + u10; s1 = a1 - a0; s2 = b1 - b0; s3 = s2 + e0;
  iv = inv(d0 * s3 - d1 * s1, p);
  ss2 = u20 + u11 * u21 + u10; ss3 = u21 + u11;
  l3 = inv(d0 * c1 - d1 * c0, p); l2 = inv(c0 * s3 - c1 * s1, p);
  l1 = u11 * l2 + v11 + (b0 - u10) * l3;
  u31 = 2 * l2 * l3 + 1 - ss3;
  u30 = 2 * l1 * l3 + l2 * l2 + 1 - f4 - u31 * ss3 - ss2;
  ul0 = u30 * l3; ul = -u31 * l2 + l1;
  v31 = -(u31 * u31 - ul0 + ul); v30 = -(u31 * ul0 + ul);
 
}

void g2dbl(ZZ u1,ZZ u0,ZZ v1,ZZ v0,ZZ f3,ZZ f2,ZZ f1,ZZ f0,ZZ p){
ZZ uu1, uu0, uv01, uv00, uv10, uv11, d0, d1, d2, d3, d4, d5, uue0, uue1, s1, s2, s3, ss2, ss3, v10, v11, v20, v21, v30, v31, ve0, ve1, l0, l1, l2, l3, iv, e0, e1, m0, m1, m3, m4, ue0, ue1;
  uu1 = u1 * u1; uu0 = u1 * u0; uv01 = u0 * v1; uv10 = u1 * v0; uv11 = u1 * v1;
  uv00 = u0 * v0;
  d0 = 6 * v1 * uu1 - uv10; d1 = -4 * uv11 + 4 * v0; d2 = 2 * v1;
  d3 = 6 * v1 * uu0 - 6 * uv00; d4 = -4 * uv01; d5 = 2 * v0;
  e0 = 5 * (-u1 * uu1 + 2 * uu0) - 3 * f3 * u1 + 2 * f2;
  e1 = 5 * (-u0 * uu0 + u0 * u0) - 3 * f3 * u0 + f1;
  m0 = d3 - d5 * (uu1 - u0); m1 = d4 - d5 * (-u1);
  m3 = d0 - d2 * (uu1 - u0); m4 = d1 - d2 * (-u1);
  s1 = e1 - d5 * v1; s2 = e0 - d2 * v1;
  iv = inv(m0 * m4 - m1 * m3, p);
  l3 = inv(m4 * s1 - m1 * s2, p); l2 = inv(m0 * s2 - m3 * s1, p);
  l1 = v1 + u1 * l2 - (uu1 - u0) * l3; l0 = v0 + u0 * l2 - (uu0) * l3;
  ue1 = 2 * l3 * l2 - 2 * u1 - 1;
  ue0 = 2 * l3 * l1 + l2 * l2 - 2 * u0 - uu0 - 2 * ue1 * u1;
  uue1 = ue1 * ue1; uue0 = ue1 * ue0;
  ve1 = (uu1 - ue0) * l3 - ue1 * l2 + l1;
  ve0 = uue0 * l3 - ue0 * l2 + l0;

//  return ue1, ue0, ve1, ve0
}

void mktable(ZZ u,ZZ v,ZZ q){
  // print ZZ q ,"\n"
  ZZ uu = u;
  ZZ vv = v;

  // enzan table
  le_u[0] = uu[0];
  le_v[0] = vv[0];
//  print ZZ uu, " l673\n"
//  print ZZ vv, " l674\n"

  for(int i = 1;i<N;i++){ //begin Pub_key at plain
    g2dbl(uu[1], uu[0], vv[1], vv[0], FF[0], FF[1], FF[2], FF[3], q);
    le_u[i] = uu[i];
    le_v[i] = vv[i];

  } //of for
//  print i, "i\n"
}


//D'=mD
void jac(ZZ kk,ZZ p){
  ZZ j = to_ZZ("0");
  Z i;
  ZZ ki[N]={0};
  ZZ U,V;

  for(j=0;j<N+1;j++){
    ki[j] = to_ZZ("0");
  }
  j = to_ZZ("0");
  for (i=0;i<N+1;i++){
    if (((kk ^ (1 << i)) >> i) % 2 == 0){ //testbit(kk,i)
      ki[j] = i;
      j = j + 1
      //print i, "L704\n"
    }
  }

   U[0] = le_u[ki[0]];
   V[0] = le_v[ki[0]];

  //print ZZ le_u[ki[0]],"\n"
  //print ZZ le_u[ki[2]],"\n"
  //exit();
  for (i = 1;i<j;i++){
    if (U != le_u[ki[i]]){
      g2add(ZZ U[0], ZZ U[1], ZZ le_u[ki[i]][0], ZZ le_u[ki[i]][1], ZZ V[0], ZZ V[1], ZZ le_v[ki[i]][0], ZZ le_v[ki[i]][1], 0, 0, 0, ZZ FF[0], ZZ FF[1], ZZ FF[2], ZZ FF[3], ZZ FF[4], p);
      //print "i=",i," ",ZZ le_u[ki[i]][0],"\n"
    }

    if (ZZ U == le_u[ki[i]]){
      cout << le_u[ki[i]] << ", =" << U << endl;
      cout << i << "de('A`)\n";
      exit(1);
    }
  } //of for
}

int main(){

HEC();

mktable(ug0, vg0, P);
jac(Jga, P);

return 0;
}