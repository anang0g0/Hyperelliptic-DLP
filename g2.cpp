#include <NTL/ZZ.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define N 256

NTL_CLIENT

ZZ le_u[2][N];
ZZ le_v[2][N];
ZZ uu[2];
ZZ vv[2];
ZZ U[2];
ZZ V[2];

ZZ FF[5];// ={0};
ZZ ue1,ue0,ve0,ve1,u31,u30,v31,v30;
ZZ ug0[2];
ZZ vg0[2];

ZZ P;
ZZ Jga;


/* definition of a curve */
void HEC(){

    #find by Harley
  @q = 10 ** 19 + 51
  @a = [3141592653589793238, 4626433832795028841, 9716939937510582097, 4944592307816406286, 2089986280348253421]
  #@u1_=[13131182302866750318,6953593084278582387]
  #@v1_=[@q,0] #infinity
  @uh = [8940387226809150403, 3838225076702837943]
  @vh = [8035450087728851271, 1893861348804881148]
  # ff=x^5+314159265358979338*x^4+4626433832795028841*x^3+9716939937510582097*x^2+4944592307816406286*x+2089986280348253421
  @u1 = [10027301878627002813, 9681764174062850433]
  #616419646419685014=a*3542790122851877922+b
  #0=a*6484511755775124891+b
  #616419646419685014=a*7058278367076753082
  @v1 = [9406915506559133975, 920961725690419616]
  @u2 = [15109848135481867673, 5563304430399854240]
  #2935061693073737419=a*5239897978117534135+b
  #3524464046627319761=a*9869950157364333538+b
  #589402353553582342=a*4630052179246799403
  @v2 = [7250939689363649434, 6461431514924022130]
  @J = 99999999982871020671452277000281660080
  #? 7054215880371151972602291562049
  s1 = 1712898036
  s2 = 11452277089352355350

  @b = [1597, 1041, 5503, 6101, 1887]
  @u = [1571353025997967, 12198441063534328]
  @v = [32227723250469108, 67133247565452990]
  @u0 = [70887725815800572, 94321182398888258]
  @v0 = [42016761890161508, 3182371156137467]

  #find by Lange
  # f=xx^5+153834295433461683634059*xx^3+1503542947764347319629935*xx^2+1930714025804554453580068*xx+790992824799875905266969
  @p3 = 1932005208863265003490787
  @F1 = [0, 153834295433461683634059, 1503542947764347319629935, 1930714025804554453580068, 790992824799875905266969]
  @uu0 = [1594018975878036024296315, 52552598504459997856285]
  @uu1 = [1791061143796384566472590, 160038959612724914387201]
  #504894935863953268767725=a*106028591185649525291891+b
  #704210062398295465154981=a*1487990384692386499004424+b
  @vv0 = [1288294670775269897135356, 942599769284250370891960]
  #1=a*704665787761008893641614+b
  #824513992484349685277891=a*1086395356035375672830976+b
  @vv1 = [325694573428709528176000, 1410279347067078324080410]
  @p4 = 3713820117856140824697372689
  # ff=x^5+241216435998068557682742515*x^3+553011586465186980114036462*x^2+1456621446251091989731057514*x+3440013483680364963850133535
  @F0 = [0, 241216435998068557682742515, 553011586465186980114036462, 1456621446251091989731057514, 3440013483680364963850133535]
  @uua = [3090907731099435637713212933, 3430279740253146837327450789]
  #1090190095529845640563737560=a*901573235033529767913345809+b
  #1607209110223949051233778495=a*2189334496065905869799867124+b
  @vva = [1516936926385660062377077660, 177161035616247877668423903]
  #1044171804289905858226438503=a*10486172923382208811538764+b
  #928996521279782754969642877=a*1874123356916514825274712274+b
  #3598644834846017721440577063=a*1863637183993132616463173510
  @uub = [1884609529839897034086251038, 265492627929763696542013047]
  @vvb = [515973747200726989346030903, 2866067948417124660466300132]

  //find by Gaudry
  //fg=x^5+2682810822839355644900736*x^3+226591355295993102902116*x^2+2547674715952929717899918*x+4797309959708489673059350

  P = to_ZZ("5000000000000000008503491");

  char* t[5]={"0", "2682810822839355644900736", "226591355295993102902116", "2547674715952929717899918", "4797309959708489673059350"};

  Jga = to_ZZ("24999999999994130438600999402209463966197516075699");
  char* t2[2]={"1", "-1713538969626908355896596"}; 
  char* t3[2]={"0", "138905579055173741542118"};
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

void g2add(ZZ u11,ZZ u10,ZZ v11,ZZ v10,ZZ u21,ZZ u20,ZZ v21,ZZ v20,ZZ f4,ZZ f3,ZZ f2,ZZ f1,ZZ f0,ZZ P){
ZZ a0, a1, ue1, ul, ul0, uu0, uu1, u30, u31, b0, b1, e0, e1, c0, c1, d0, d1, ss2, ss3, s1, s2, s3, v30, v31;
ZZ iv, l2, l3, l1;
  a0 = (u11 * u10)%P; a1 = (u21 * u20)%P; d0 = (u10 - u20)%P; d1 = (u11 - u21)%P;
  b0 = (u11 * u11)%P; b1 = (u21 * u21)%P; c0 = (v20 - v10)%P; c1 = (v21 - v11)%P;
  e0 = (-u20 + u10)%P; s1 = (a1 - a0)%P; s2 = (b1 - b0)%P; s3 = (s2 + e0)%P;
  iv = inv(d0 * s3 - d1 * s1, P);
  ss2 = (u20 + u11 * u21 + u10)%P; ss3 = (u21 + u11)%P;
  l3 = inv(d0 * c1 - d1 * c0, P); l2 = inv(c0 * s3 - c1 * s1, P);
  l1 = (u11 * l2 + v11 + (b0 - u10) * l3)%P;
  u31 = (2 * l2 * l3 + 1 - ss3)%P;
  u30 = (2 * l1 * l3 + l2 * l2 + 1 - f4 - u31 * ss3 - ss2)%P;
  ul0 = (u30 * l3)%P; ul = (-u31 * l2 + l1)%P;
  v31 = -(u31 * u31 - ul0 + ul)%P; v30 = -(u31 * ul0 + ul)%P;

  U[0]=u30;
  U[1]=u31;
  V[0]=v30;
  V[1]=v31;
}

void g2dbl(ZZ u1,ZZ u0,ZZ v1,ZZ v0,ZZ f3,ZZ f2,ZZ f1,ZZ f0,ZZ p){
ZZ uu1, uu0, uv01, uv00, uv10, uv11, d0, d1, d2, d3, d4, d5, uue0, uue1, s1, s2, s3, ss2, ss3, v10, v11, v20, v21, v30, v31, ve0, ve1, l0, l1, l2, l3, iv, e0, e1, m0, m1, m3, m4, ue0, ue1;
  uu1 = (u1 * u1)%P; uu0 = (u1 * u0)%P; uv01 = (u0 * v1)%P; uv10 = (u1 * v0)%P; uv11 = (u1 * v1)%P;
  uv00 = (u0 * v0)%P;
  d0 = (6 * v1 * uu1 - uv10)%P; d1 = (-4 * uv11 + 4 * v0)%P; d2 = (2 * v1)%P;
  d3 = (6 * v1 * uu0 - 6 * uv00)%P; d4 = (-4 * uv01)%P; d5 = (2 * v0)%P;
  e0 = (5 * (-u1 * uu1 + 2 * uu0) - 3 * f3 * u1 + 2 * f2)%P;
  e1 = (5 * (-u0 * uu0 + u0 * u0) - 3 * f3 * u0 + f1)%P;
  m0 = (d3 - d5 * (uu1 - u0))%P; m1 = (d4 - d5 * (-u1))%P;
  m3 = (d0 - d2 * (uu1 - u0))%P; m4 = (d1 - d2 * (-u1))%P;
  s1 = (e1 - d5 * v1)%P; s2 = (e0 - d2 * v1)%P;
  iv = inv(m0 * m4 - m1 * m3, P);
  l3 = inv(m4 * s1 - m1 * s2, P); l2 = inv(m0 * s2 - m3 * s1, P);
  l1 = (v1 + u1 * l2 - (uu1 - u0) * l3)%P; l0 = (v0 + u0 * l2 - (uu0) * l3)%P;
  ue1 = (2 * l3 * l2 - 2 * u1 - 1)%P;
  ue0 = (2 * l3 * l1 + l2 * l2 - 2 * u0 - uu0 - 2 * ue1 * u1)%P;
  uue1 = (ue1 * ue1)%P; uue0 = (ue1 * ue0)%P;
  ve1 = ((uu1 - ue0) * l3 - ue1 * l2 + l1)%P;
  ve0 = (uue0 * l3 - ue0 * l2 + l0)%P;

  uu[1]=ue1;
  uu[0]=ue0;
  vv[1]=ve1;
  vv[0]=ve0;
//  return ue1, ue0, ve1, ve0
}

void mktable(ZZ u[],ZZ v[],ZZ q){
  // print ZZ q ,"\n"
   uu[1] = u[0];
   uu[0] = u[1];
   vv[1] = v[0];
   vv[0] = v[1];

  // enzan table
  le_u[0][0] = uu[0];
  le_v[0][0] = vv[0];
  le_u[1][0] = uu[1];
  le_v[1][0] = vv[1];
  cout << uu[0] << " l673\n";
  cout << vv[0] << " l674\n";

  for(int i = 1;i<N;i++){ //begin Pub_key at plain
    g2dbl(uu[1], uu[0], vv[1], vv[0], FF[0], FF[1], FF[2], FF[3], q);
    le_u[0][i] = uu[0];
    le_u[1][i] = uu[1];
    le_v[0][i] = vv[0];
    le_v[1][i] = vv[1];
  cout << i << "i\n";
  } //of for

}


//D'=mD
void jac(ZZ kk,ZZ p){
  int j = 0;
  int i;
  unsigned int ki[N];


  for(j=0;j<N+1;j++){
    ki[j] = 0;
  }
  j = 0;
  for (i=0;i<N+1;i++){
    if (((kk ^ (1 << i)) >> i) % 2 == 0){ //testbit(kk,i)
      ki[j] = i;
      j = j + 1;
      cout << i << "L704\n";
    }
  }

   U[0] = le_u[0][ki[0]];
   U[1] = le_u[1][ki[0]];
   V[0] = le_v[0][ki[0]];
   V[1] = le_v[1][ki[0]];

  cout << le_u[0][ki[0]] << "\n";
  cout << le_u[1][ki[0]] << "\n";
//  exit(1);

  for (i = 1;i<j;i++){
    if (U[0] != le_u[0][ki[i]] && U[1] != le_u[1][ki[i]]){
      g2add(U[0], U[1], le_u[0][ki[i]], le_u[1][ki[i]], V[0], V[1], le_v[0][ki[i]], le_v[1][ki[i]], FF[0], FF[1], FF[2], FF[3], FF[4], p);
      cout << "i=" << i << " " << le_u[0][ki[i]] << "," << U[0] << "\n";
    }

    if (U[0] == le_u[0][ki[i]]){
      cout << le_u[0][ki[i]] << ", =" << U[0] << endl;
      cout << i << "de('A`)\n";
      exit(1);
    }
  } //of for
}



//jj=aa^bb mod oo
ZZ exp(ZZ aa,ZZ bb,ZZ oo){
ZZ ii,jj,kk[8192];
int j,c[8192],count=0,i;
  ii=oo;
  j=0;
  jj=0;
//  kk[4096]; //prime is 4096 bit table
//  c[8192]  //mod is 8192 bit table
  count=0;

  for(i=0;i<8192;i++){
    kk[i]=0;
    }
  while(ii>0){
    ii = (ii>>1);
    j=j+1;
  }


  kk[0]=aa;

//  cout << j << "\n";
  
//ex.1000=2**3+2**5+2**6+2**7+2**8+2**9 makes a array c=[3,5,6,7,8,9]
  for(i=0;i<j+1;i++){
      if(bit(bb,i) != 0){ // testbit(bb,i)
	c[count]=i;
	count=count+1;
      }
    }
//    cout << bb << endl;
//    cout << count << "\n";
//exit(1);
    for(i=1;i<c[count-1]+1;i++){
      kk[i] = kk[i-1]*kk[i-1]%oo;
    }

    jj=1;
    for(i=0;i<count;i++){
      jj=kk[c[i]]*jj%oo;
      if (jj==0){
//	print i,"\n"
      }
    }

    return jj;
}


ZZ root(ZZ a, ZZ P){
  ZZ b,c,d,e;

  b=(P+1)/to_ZZ("4");
  c=exp(a,b,P);

  cout << a << endl;
  cout << c << endl;

  return c;
}


int main(){

HEC();

mktable(ug0, vg0, P);
jac(Jga, P);

return 0;
}