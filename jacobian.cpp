#include <iostream>
#include <NTL/ZZ.h>
#include <stdio.h>
#include <stdlib.h>

NTL_CLIENT

int main(){


char *a[1]={"4856350849"};
char *b[1]={"170765737"};
char *c[2]={"8402287361","293048739936698113783226667142139285987"};

char *d[1]={"3081439808"};
char *e[1]={"1398375548"};
char *f[2]={"6183226511","87584922866592430046826643578455405177"};

ZZ p=to_ZZ("340282366920938463463374607431768211283");
ZZ q=to_ZZ("170141183460469231731687303715884105727");
ZZ a1[2],a2[2],a3[2];
ZZ b1[2],b2[2],b3[2];

a1[0]=to_ZZ(a[0]);
a1[1]=to_ZZ("1");
a2[0]=to_ZZ(b[0]);
a2[1]=to_ZZ("1");
b1[0]=to_ZZ(d[0]);
b1[1]=to_ZZ("1");
b2[0]=to_ZZ(e[0]);
b2[1]=to_ZZ("1");

for(int i=0;i<2;i++){
    a3[i]=to_ZZ(c[i]);
    b3[i]=to_ZZ(f[i]);
    }

ZZ tmp1[3];
ZZ tmp2[3];
for(int i=0;i<3;i++){
    tmp1[i]=to_ZZ("0");
    tmp2[i]=to_ZZ("0");
}
for(int i=0;i<2;i++){
    for(int j=0;j<2;j++){
        tmp1[i+j]+=(a1[i]*a2[j])%p;
        tmp2[i+j]+=(b1[i]*b2[j])%q;
    }
}
for(int i=0;i<3;i++){
    tmp1[i]=tmp1[i]%p;
    tmp2[i]=tmp2[i]%q;
}

for(int i=0;i<3;i++)
cout << tmp1[i] << endl;
for(int i=0;i<3;i++)
cout << tmp2[i] << endl;
/*
cout << (b2[0]*b1[0])%q << endl;
cout << (b2[1]*b1[0]+b1[1]*b2[0])%q << endl;
cout << (a2[0]*a1[0])%p << endl;
cout << (a2[1]*a1[0]+a1[1]*a2[0])%p << endl;
*/

return 0;
}
