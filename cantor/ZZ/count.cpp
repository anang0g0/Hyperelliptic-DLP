#include <NTL/ZZ.h>
#include <stdio.h>

NTL_CLIENT


  // ZZ f[K+1] = {to_ZZ("1"),to_ZZ("0"), to_ZZ("1184"), to_ZZ("1846"), to_ZZ("956"),  to_ZZ("560")};
/*
  ZZ J = to_ZZ("115792089237316195429342203801033554170931615651881657307308068079702089951781");
  ZZ f[K + 1] = {to_ZZ("1"), to_ZZ("0"), to_ZZ("318258242717201709453901384328569236653"), to_ZZ("75380722035796344355219475510170298006"), to_ZZ("129416082603460579272847694630998099237"), to_ZZ("143864072772599444046778416709082679388")};
*/
  /*
    ZZ f[K+1] = {to_ZZ("1"),to_ZZ("3141592653589793238"), to_ZZ("4626433832795028841"), to_ZZ("9716939937510582097"), to_ZZ("4944592307816406286"), to_ZZ("2089986280348253421")};
    //#@u1_=[13131182302866750318,6953593084278582387]
    //#@v1_=[@q,0] #infinity
    ZZ u1[K+1] = {to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("1"),to_ZZ("8940387226809150403"), to_ZZ("3838225076702837943")};
    ZZ v1[K+1] = {to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("8035450087728851271"), to_ZZ("1893861348804881148")};
    // ff=x^5+314159265358979338*x^4+4626433832795028841*x^3+9716939937510582097*x^2+4944592307816406286*x+2089986280348253421
   ZZ J = to_ZZ("99999999982871020671452277000281660081");
   */
  /*
    ZZ uh[K+1] = {0,0,0,1,10027301878627002813, 9681764174062850433};
    //#616419646419685014=a*3542790122851877922+b
    //#0=a*6484511755775124891+b
    //#616419646419685014=a*7058278367076753082
    ZZ vh[K+1] = {0,0,0,0,9406915506559133975, 920961725690419616};
    ZZ u2[K+1] = {0,0,0,1,15109848135481867673, 5563304430399854240};
    //#2935061693073737419=a*5239897978117534135+b
    //#3524464046627319761=a*9869950157364333538+b
    //#589402353553582342=a*4630052179246799403
    ZZ v2[K+1] = {0,0,0,0,7250939689363649434, 6461431514924022130};
  */
  //  #find by Gaudry
  //  #fg=x^5+2682810822839355644900736*x^3+226591355295993102902116*x^2+2547674715952929717899918*x+4797309959708489673059350
  //  @P = 5 * 10 ** 24 + 8503491
  /*
     ZZ f[K+1] = {to_ZZ("1"),to_ZZ("0"), to_ZZ("2682810822839355644900736"), to_ZZ("226591355295993102902116"), to_ZZ("2547674715952929717899918"), to_ZZ("4797309959708489673059350")};
     ZZ Jga = to_ZZ("24999999999994130438600999402209463966197516075699");
     ZZ Jgb = to_ZZ("25000000000005869731468829402229428962794965968171");
      //#50724386855111482309402*2779199501981512279739817%P
      //#2055622596816515886446193*1553122609714208136553134%P
      //#1520505942073936921231867
      ZZ u1[K+1] = {to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("1"),to_ZZ("1738366273424896804842766"), to_ZZ("3184841659043138633535652")};
      ZZ v1[K+1] = {to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("0"),to_ZZ("2931056213155642836850986"), to_ZZ("402980510097415333052905")};
    */

  //  ZZ f[K + 1] = {to_ZZ("1"), to_ZZ("7"), to_ZZ("6"), to_ZZ("2"), to_ZZ("8"), to_ZZ("2")};
  /*
      char  u2[K + 1] = {0, 0, 0, 1, 21, 16};
      char  u1[K + 1] = {0, 0, 0, 1, 19, 20};
      char  v2[K + 1] = {0, 0, 0, 0, 21, 21};
      char  v1[K + 1] = {0, 0, 0, 0, 12, 8};
    */
  // ZZ f[K + 1] = {to_ZZ("1"), to_ZZ("0"), to_ZZ("2"), to_ZZ("30"), to_ZZ("5"), to_ZZ("1")};
  /*

    unsigned long long u2[K + 1] = {0, 0, 0, 1, 26, 20};
    unsigned long long v2[K + 1] = {0, 0, 0, 0, 29, 26};
    unsigned long long u1[K + 1] = {0, 0, 0, 1, 9, 27};
    unsigned long long v1[K + 1] = {0, 0, 0, 0, 29, 16};

    unsigned long long u2[K + 1] = {0, 0, 0, 1, 30, 3};
    unsigned long long v2[K + 1] = {0, 0, 0, 0, 12, 8};
    unsigned long long u1[K + 1] = {0, 0, 0, 1, 0, 12};
    unsigned long long v1[K + 1] = {0, 0, 0, 0, 10, 4};
  */
  //  unsigned long long  f[K + 1] = {1, 0, 3, 7, 1, 2};
  /*
    unsigned long long  u2[K + 1] = {0, 0, 0, 1, 7, 10};
    unsigned long long  u1[K + 1] = {0, 0, 0, 1, 0, 10};
    unsigned long long  v2[K + 1] = {0, 0, 0, 0, 1, 9};
    unsigned long long  v1[K + 1] = {0, 0, 0, 0, 7, 9};
  */
//ZZ JC =to_ZZ("115792089237316195413084848286254559191634403377276171765138133738613798321603");

//ZZ aa = to_ZZ("1331");
//ZZ J = to_ZZ("1461501637326815988079848163961117521046955445901");
//ZZ P = to_ZZ("1208925819614629174709941");
//ZZ a = to_ZZ("2");
//ZZ JC = to_ZZ("1461501637331762771847359428275278989652932675771");



/* definition of a curve */

void HEC()
{
  /*
  //p2 = 2 ** 128 - 173
  ZZ P2=to_ZZ("340282366920938463463374607431768211283");
  ZZ r2 = to_ZZ("115792089237316195429342203801033554170931615651881657307308068079702089951781");
  ZZ ef[K+1] = {1,0,318258242717201709453901384328569236653, 75380722035796344355219475510170298006, 129416082603460579272847694630998099237, 143864072772599444046778416709082679388};

  //2**127-1
  P3=to_ZZ("170141183460469231731687303715884105727");
  ZZ p=to_ZZ("170141183460469231731687303715884105727");
  ZZ r = to_ZZ("28948022309329048848169239995659025138451177973091551374101475732892580332259");
  //which is 254 bits. A possible degree 5 model is C : y2 = x**5 + f3*x**3 + f2*x**2 + f1x + f0, where
  ZZ f3[K+1] = {1,0,34744234758245218589390329770704207149,  132713617209345335075125059444256188021, 90907655901711006083734360528442376758,  6667986622173728337823560857179992816};


  // p127m = (2**63 - 27443) · 2**64 + 1.
  ZZ P4=to_ZZ("170141183460468725497689688904659107841");
  //f:y^2 = x^5 + 17
  fa[K+1]={1,0,0,0,0,17};
  ZZ r2 = to_ZZ("28948022309328876595115567994214488524823328209723866335483563634241778912751");


  ZZ P5=to_ZZ("28=340282366920938463463374607431768186521");
  // y2 = x5+ 37
  fb[K+1]={1,0,0,0,0,37};
  ZZ r3 = to_ZZ("115792089237316195401210495125503591471546519982099914586091636775415022457661");
  */

  /*
    //find by Harley
  //  @q = 10 ** 19 + 51
  q=10000000000000000051;
    ZZ a[5] = [3141592653589793238, 4626433832795028841, 9716939937510582097, 4944592307816406286, 2089986280348253421]
    //#@u1_=[13131182302866750318,6953593084278582387]
    //#@v1_=[@q,0] #infinity
    ZZ uh = [8940387226809150403, 3838225076702837943]
    ZZ vh = [8035450087728851271, 1893861348804881148]
    // ff=x^5+314159265358979338*x^4+4626433832795028841*x^3+9716939937510582097*x^2+4944592307816406286*x+2089986280348253421
    ZZ u1 = [10027301878627002813, 9681764174062850433]
    //#616419646419685014=a*3542790122851877922+b
    //#0=a*6484511755775124891+b
    //#616419646419685014=a*7058278367076753082
    ZZ v1 = [9406915506559133975, 920961725690419616]
    ZZ u2 = [15109848135481867673, 5563304430399854240]
    //#2935061693073737419=a*5239897978117534135+b
    //#3524464046627319761=a*9869950157364333538+b
    //#589402353553582342=a*4630052179246799403
    ZZ v2 = [7250939689363649434, 6461431514924022130]
    ZZ J = 99999999982871020671452277000281660080
    //#? 7054215880371151972602291562049


    #find by Lange
    //# f=xx^5+153834295433461683634059*xx^3+1503542947764347319629935*xx^2+1930714025804554453580068*xx+790992824799875905266969
    ZZ p3 = 1932005208863265003490787
    ZZ F1 = [0, 153834295433461683634059, 1503542947764347319629935, 1930714025804554453580068, 790992824799875905266969]
    ZZ uu0 = [1594018975878036024296315, 52552598504459997856285]
    ZZ uu1 = [1791061143796384566472590, 160038959612724914387201]
    //#504894935863953268767725=a*106028591185649525291891+b
    //#704210062398295465154981=a*1487990384692386499004424+b
    ZZ vv0 = [1288294670775269897135356, 942599769284250370891960]
    //#1=a*704665787761008893641614+b
    //#824513992484349685277891=a*1086395356035375672830976+b
    ZZ vv1 = [325694573428709528176000, 1410279347067078324080410]
    ZZ p4 = 3713820117856140824697372689
    //# ff=x^5+241216435998068557682742515*x^3+553011586465186980114036462*x^2+1456621446251091989731057514*x+3440013483680364963850133535
    ZZ F0 = [0, 241216435998068557682742515, 553011586465186980114036462, 1456621446251091989731057514, 3440013483680364963850133535]
    ZZ uua = [3090907731099435637713212933, 3430279740253146837327450789]
    //#1090190095529845640563737560=a*901573235033529767913345809+b
    //#1607209110223949051233778495=a*2189334496065905869799867124+b
    ZZ vva = [1516936926385660062377077660, 177161035616247877668423903]
    //#1044171804289905858226438503=a*10486172923382208811538764+b
    //#928996521279782754969642877=a*1874123356916514825274712274+b
    //#3598644834846017721440577063=a*1863637183993132616463173510
    ZZ uub = [1884609529839897034086251038, 265492627929763696542013047]
    ZZ vvb = [515973747200726989346030903, 2866067948417124660466300132]
  */
  //  #find by Gaudry
  //  #fg=x^5+2682810822839355644900736*x^3+226591355295993102902116*x^2+2547674715952929717899918*x+4797309959708489673059350
  //  @P = 5 * 10 ** 24 + 8503491
  ZZ P6 = to_ZZ("5000000000000000008503491");
  ZZ FF[K + 1] = {to_ZZ("1"), to_ZZ("0"), to_ZZ("2682810822839355644900736"), to_ZZ("226591355295993102902116"), to_ZZ("2547674715952929717899918"), to_ZZ("4797309959708489673059350")};
  ZZ Jga = to_ZZ("24999999999994130438600999402209463966197516075699");
  ZZ Jgb = to_ZZ("25000000000005869731468829402229428962794965968171");
  //#50724386855111482309402*2779199501981512279739817%P
  //#2055622596816515886446193*1553122609714208136553134%P
  //#1520505942073936921231867
  ZZ ug0[K + 1] = {to_ZZ("0"), to_ZZ("0"), to_ZZ("0"), to_ZZ("1"), to_ZZ("1"), to_ZZ("-1713538969626908355896596")};
  ZZ vg0[K + 1] = {to_ZZ("0"), to_ZZ("0"), to_ZZ("0"), to_ZZ("0"), to_ZZ("0"), to_ZZ("138905579055173741542118")};
  ZZ ug1[K + 1] = {to_ZZ("0"), to_ZZ("0"), to_ZZ("0"), to_ZZ("1"), to_ZZ("1738366273424896804842766"), to_ZZ("3184841659043138633535652")};
  ZZ vg1[K + 1] = {to_ZZ("0"), to_ZZ("0"), to_ZZ("0"), to_ZZ("0"), to_ZZ("2931056213155642836850986"), to_ZZ("402980510097415333052905")};
}


int main(){
int i=0;
ZZ T=to_ZZ("187145766228145807790746688818282675710128013948861177145226708619702335205202779181217");

while(T>0){
    T=(T>>1);
    i++;
}

printf("%d\n",i);

return 0;
}

