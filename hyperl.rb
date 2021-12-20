#/* HyperElliptic Curve Domain of Parametor */
# Copyright(C) 2005 by tcshacina (20050314)
# This is a example of a HECDLP imprementation
# hyperl.rb ver0.1
# ?????Y????(????)
# 1.?}???t?H?[?h?\????????????q???o??????????B
# 2.???R?r?A?????????????????????(v=0)????????????B
# 3.????????????????????????????????????B
# 4.?????v?Z????????????????H
# g=3,4,5?????????????????�??A?b?v?B
# ???B


N=256
@le_u=[N,N,N,N]
@le_v=[N,N,N,N]
@uu=[N,N,N,N]
@vv=[N,N,N,N]
@U=[N,N,N,N]
@V=[N,N,N,N]

@uu4=[N,N,N,N]
@vv4=[N,N,N,N]
@U4=[N,N,N,N]
@V4=[N,N,N,N]

# invert of integer
def inv(a, n)
  
  d = n
  x = 0
  s = 1
  while (a != 0)
    q = d / a
    r = d % a
    d = a
    a = r
    t = x - q * s
    x = s
    s = t
  end
  gcd = d  # $\gcd(a, n)$ 

  return ((x + n) % (n / d))
end



# /* definition of a curve */
def HEC()
  #find by Smart
  #y^2+xy=x^3+b1B^2+b2B^4+b3B^8
  @pa=1208925819614311295169073
  @b1=1127280
  @b2=171398
  @b3=1370436
  
  #y^2+(x^4+a1x^3+a2x^2+a3x+a4)y+(x^9+a5x^6+a6x^4+a7x^3+a8x^2+a9x+e)
  @pb=1208925819614311295169073
  @a1=624429
  @a2=1248858
  @a3=1442662
  @a4=386860
  @a5=1859582
  @a6=293124
  @a7=1783647
  @a8=1541982
  @a9=1370912
  @e=1888298
  
  
  @pc=584600649323611672814739995379292203636332479268
  #y^2+(x^2+c1x)y=x^5+x^4+c2x^2+c3x+c4
  c1=2012013793551629036365609
  c2=1586464037343056940725724
  c3=43334222987849600951547
  c4=774788345987798314632240
  
  
  #find by Harley
  @q=10**19+51
  @a=[3141592653589793238,4626433832795028841,9716939937510582097,4944592307816406286,2089986280348253421]
  #@u1_=[13131182302866750318,6953593084278582387]
  #@v1_=[@q,0] #infinity
  @uh=[8940387226809150403, 3838225076702837943]
  @vh=[8035450087728851271, 1893861348804881148]
  # ff=x^5+314159265358979338*x^4+4626433832795028841*x^3+9716939937510582097*x^2+4944592307816406286*x+2089986280348253421
  @u1=[10027301878627002813,9681764174062850433]
  #616419646419685014=a*3542790122851877922+b
  #0=a*6484511755775124891+b
  #616419646419685014=a*7058278367076753082
  @v1=[9406915506559133975,920961725690419616]
  @u2=[15109848135481867673,5563304430399854240]
  #2935061693073737419=a*5239897978117534135+b
  #3524464046627319761=a*9869950157364333538+b
  #589402353553582342=a*4630052179246799403
  @v2=[7250939689363649434,6461431514924022130]
  @J=99999999982871020671452277000281660080
  #? 7054215880371151972602291562049
  s1=1712898036
  s2=11452277089352355350
  
  @b=[1597,1041,5503,6101,1887]
  @u=[1571353025997967,12198441063534328]
  @v=[32227723250469108,67133247565452990]
  @u0=[70887725815800572,94321182398888258]
  @v0=[42016761890161508,3182371156137467]
  
  #find by Lange
  # f=xx^5+153834295433461683634059*xx^3+1503542947764347319629935*xx^2+1930714025804554453580068*xx+790992824799875905266969
  @p3=1932005208863265003490787
  @F1=[0,153834295433461683634059,1503542947764347319629935,1930714025804554453580068,790992824799875905266969]
  @uu0=[1594018975878036024296315,52552598504459997856285]
  @uu1=[1791061143796384566472590,160038959612724914387201]
  #504894935863953268767725=a*106028591185649525291891+b
  #704210062398295465154981=a*1487990384692386499004424+b
  @vv0=[1288294670775269897135356,942599769284250370891960]
  #1=a*704665787761008893641614+b
  #824513992484349685277891=a*1086395356035375672830976+b
  @vv1=[325694573428709528176000,1410279347067078324080410]
  @p4=3713820117856140824697372689
  # ff=x^5+241216435998068557682742515*x^3+553011586465186980114036462*x^2+1456621446251091989731057514*x+3440013483680364963850133535
  @F0=[0,241216435998068557682742515,553011586465186980114036462,1456621446251091989731057514,3440013483680364963850133535]
  @uua=[3090907731099435637713212933,3430279740253146837327450789]
  #1090190095529845640563737560=a*901573235033529767913345809+b
  #1607209110223949051233778495=a*2189334496065905869799867124+b
  @vva=[1516936926385660062377077660,177161035616247877668423903]
  #1044171804289905858226438503=a*10486172923382208811538764+b
  #928996521279782754969642877=a*1874123356916514825274712274+b
  #3598644834846017721440577063=a*1863637183993132616463173510
  @uub=[1884609529839897034086251038,265492627929763696542013047]
  @vvb=[515973747200726989346030903,2866067948417124660466300132]
  
  
  #find by Gaudry
  #fg=x^5+2682810822839355644900736*x^3+226591355295993102902116*x^2+2547674715952929717899918*x+4797309959708489673059350
  @P=5*10**24+8503491
  @FF=[0,2682810822839355644900736,226591355295993102902116,2547674715952929717899918,4797309959708489673059350]
  @Jga=24999999999994130438600999402209463966197516075699
  #@Jgb=25000000000005869731468829402229428962794965968171
  #50724386855111482309402*2779199501981512279739817%P
  #2055622596816515886446193*1553122609714208136553134%P
  #1520505942073936921231867
  @ug0=[1,-1713538969626908355896596]
  @vg0=[0,138905579055173741542118]
  @ug1=[1738366273424896804842766,3184841659043138633535652]
  @vg1=[2931056213155642836850986,402980510097415333052905]
  
#find by takahashi (g=4)
#y2 = x^9 + 29x, 
@ft=[0,0,0,0,0,0,0,29,0]
@g4p = 1759218504481 #(41-bit)
@Jta=4789034620376653463540859489797855263219497047089 #(162-bit)
@ut1=[1583083862007,1691057975742,1024593740896,782130844060]
@vt1=[1558543996506,719475088010,1650902695840,1703257876129]
#y2 = x^9 + 1953125x, 
@ft1=[0,0,0,0,0,0,0,1953125,0]
@g4p2 = 2199023315233 #(41-bit)
@Jtb=11691986799636433497742258013292719544703684675777 #(163-bit)
# g=2
@fp=1208925819614629175095961
@ff = [1,0,0,0,3,0]
@fu = [1,1];
@fv = [0,2];
@order = 1461501637332961738997140052922587693332824903682

#find by waseda
@p2=9671406556917033397660861
@J2=93536104789297918652114038095154103630376284213875
@p3=631153760340591307
@J3=251423300188980936808533376314868064530443303970434811
@wu=[377166208405888276, 65838350729795034, 258349841662486010]
@wv=[82817318614465843, 235991131611388836, 245978366599867468]
@p6=269049691
@J6=379316622800922233782741202725478330656627788904081

#test by me
@Q=31
@Fe=[1,1,1,1,1]
@ug=[18,5]
@vg=[10,23]
#g=3 y^2=x^5+2988249846838211
@ord=1707775189051778189231079919462521312159759509959
@Pa=11952999387352843
@then=478 
#g=5 over GF(2^29)
#C:110=101000000001
@ord1=1498314861672848797366241841475989276779
#g=4 over GF(2^31)
#C:110=1111100001
@ord2=924789699702887396540027854335006631
#g=3 over GF(2^61)
#C:110=10000011
@orz=401900157961648362039235675748707122986758302471449
#g=6 over GF(2^31)
#C:110=11111010000001
@org=3382155885826628920558461386553295741968057149350507061



=begin
a=100003
p=a**4+a**3+a**2+a+1
print p,"\n"

p=88741
for x in 0..p-1
 for y in 0..p-1
  if((y*y+y)%p==x**5%p)
   print "(",x,",",y,")\n"
  end
 end
end


a=100012
p=a**4+a**3+a*a+a+1
print p,"\n"

a=100018
p=a**4+a**3+a*a+a*1
print p,"\n"

=end

end



def cantora(a,b)

d1=gcd(a,b);
d2=e1*a+e2*b;

d=gcd(d1,b1+b2+h);
dd=c1*d1+c2*(b1+b2+h);
s1=c1*e1;
s2=c1*e2;
s3=c2;

d=s1*a1+s2*a2+s3*(b1+b2+h);
a=a1*a2/(d*d);

b=inv(s1a1b2+s2a2b1+s3(b1b2+f,d,a));

end


def cantord(a,b)

d=gcd(a,2*b+h) #=s1*a+s3*(2*b+h)
a0=a*a/d**2
b0=((s1*a*b+s3*(b*b+f))/d)%a0

end


def g2add(u11,u10,u21,u20,v11,v10,v21,v20,h2,h1,h0,f4,f3,f2,f1,f0,p)
z1=(u11-u21)%p; z2=(u20-u10)%p; z3=(u11*z1+z2)%p;
r=(z2*z3+(z1**2)*u10)%p

inv1=z1; inv0=z2;
w1=(v10-v20)%p; w2=(v11-v21)%p; w3=inv0*w1%p; w4=inv1*w2%p;
s1_=((inv0+inv1)*(w1+w2)-w3-w4*(1+u11))%p; s0_=(w3-u10*w4)%p;

=begin
if s1==0
inb=inv(r,p)

s0=s0_*inb%p
l1=s0*u21%p
l0=s0*u20%p

u0_=(f4-u21-u11-s0**2-s0*h2)%p

w1=(l1+v21+h1)%p
w2=(s0+v20+h0)%p
w1=(w1-u0_*(h2+s0))%p
v0_=(u0_*w1-w2)%p

=end

w1=inv(r*s1_,p)
w2=r*w1%p
w3=(s1_**2)*w1%p
w4=r*w2%p
w5=(w4**2)%p
s0__=s0_*w2%p

l2_=(u21+s0__)%p
l1_=(u21*s0__+u20)%p
l0_=u20*s0__%p

u0_=((s0__-u11)*(s0__-z1+h2*w4)-u10+l1_+(h1+2*v21)*w4+(2*u21-z1-f4)*w5)%p
u1_=(2*s0__-z1-w5+h2*w4)%p
@U=[u1_,u0_]

w1=(l2_-u1_)%p
w2=(u1_*w1+u0_-l1_)%p
v1_=(w2*w3-v21-h1+h2*u1_)%p
w2=(u0_*w1-l0_)%p
v0_=(w2*w3-v20-h0+h2*u0_)%p
@V=[v1_,v0_]
print "diva (s^2+",u1_,"*s+",u0_,",",v1_,"*s+",v0_,")\n"

end


def g2ads(u11,u10,u21,u20,v11,v10,v21,v20,h2,h1,h0,f4,f3,f2,f1,f0,p)

r=(u21-(u21*u10)*u10)%p

inb=inv(r,@p3)

s0=inb*(v10-v20-v21*u10)%p

l1=s0*u21%p
l0=s0*u20%p

k2=(f4-u21)%p
k1=(f3-(f4-u21)*u21-u20-v21*h2)%p

u1_=(k2-s0**2-u10-s0*h2)%p
u0_=(k1-s0*(l1+2*v21+h1)-u10*u1_)%p

v1_=((h2+s0)*u1_-(h1+l1+v21))%p
v0_=((h2+s0)*u0_-(h0+l0+v20))%p

end


def g2dbl(u1,u0,v1,v0,h2,h1,h0,f4,f3,f2,f1,f0,p)

v1goal=(h1+2*v1-h2*u1)%p
v0goal=(h0+2*v0-h2*u0)%p

w0=v1**2%p
w1=u1**2%p
w2=v1goal**2%p
w3=u1*v1goal%p
r=(u0*w2+v0goal*(v0goal-w3))%p

inv1_=-v1goal%p
inv0_=(v0goal-w3)%p

w3=(f3+w1)%p
w4=2*u0%p
k1_=(2*(w1-f4*u1)+w3-w4-v1*h2)%p
k0_=(u1*(2*w4-w3+f4*u1+v1*h2)+f2-w0-2*f4*u0-v1*h1-v0*h2)%p

w0=k0_*inv0_%p
w1=k1_*inv1_%p
s1_=((inv0_+inv1_)*(k0_+k1_)-w0-w1*(1+u1))%p
s0_=(w0-u0*w1)%p

=begin
if s1==0
w1=inv(r,p)
s0=s0_*w1%p
w1=(u1*s0+v1+h1)%p
w2=(u0*s0+v0+h0)%p

u0_=(f4-s0**2-2*u1-s0*h2)%p

w1=(w1-u0_*(s0+h2))%p
v0_=(u0_*w1-w2)%p

=end

w1=inv(r*s1_,p)
w2=r*w1%p
w3=s1_**2*w1%p
w4=r*w2%p
w5=w4**2%p
s0__=s0_*w2%p

l2_=(u1+s0__)%p
l1_=(u1*s0__+u0)%p
l0_=u0*s0__%p

u0_=(s0__**2+w4*(2*v1+h1+h2*(s0__+u1))+w5*(2*u1-f4))%p
u1_=(2*s0__-w5+w4*h2)%p
@uu=[u1_,u0_]

w1=(l2_-u1_)%p
w2=(u1_*w1+u0_-l1_)%p
v1_=(w2*w3-v1-h1-u1_*h2)%p
w2=(u0_*w1-l0_)%p
v0_=(w2*w3-v0-h0)%p
@vv=[v1_,v0_]
print "div (s^2+",u1_,"*s+",u0_,",",v1_,"*s+",v0_,")\n"

end


def chao(u11,u10,v11,v10,f4,f3,f2,p)

w1=v11**2%p
w2=u11*v11%p
r=(u10*w1+v10*(v10-w2))%p

if r==0
 print "r=0\n"
 u11=(v10*inv(v11,p)-u11)%p
 v11=-v10*inv(v11,p)%p
 w1=v11**2%p
 w2=u11*v11%p
 r=(u10*w1+v10*(v10-w2))%p
end

w3=inv(2*r,p)
i11=-v11*w3%p
i10=(v10-w2)*w3%p
#print "I=",i11,"x+",i10,"\n"

w2=(u11-f4)%p; w3=2*u10%p;
t10=(u11*(2*w3-u11*w2-f3)-f4*w3+f2-w1)%p;
t11=(u11*(2*w2+u11)+f3-w3)%p;

w1=i10*t10%p; w2=i11*t11%p;
w3=((i10+i11)*(t10+t11)-w1-w2)%p;
s1=(w3-u11*w2)%p
s0=(-u10*w2+w1)%p

if s1==0
 print "s1=0\n"
 exit();
end

w1=inv(s1,p)
u20=(w1*(w1*(s0**2+2*u11*-f4)+2*v11))%p;
u21=w1*(2*s0-w1)%p; u22=1;
@uu=[u21,u20]
print "factormod(s^2+",u21,"*s+",u20,",",p,")\n"

w1=(u11-u21)%p
v20=(u20*(s1*w1+s0)-s0*u10-v10)%p;
v21=(s1*(u21*w1+u20-u10)-s0*w1-v11)%p;
@vv=[v21,v20]
print "v'=",v21,"s",v20,"\n";

if(v20==0 || v21==0)
 print "infinity!\n"
 exit()
end

end



def kotehan(u21,u20,u11,u10,v21,v20,v11,v10,f4,p)

w1=(u21-u11)%p; w2=(u21*w1+u10-u20)%p;
r=(u10*(w2-u20)+u20*(u20-u11*w1))%p
if r==0
 print "r=0!\n"
 exit();
end

w3=inv(r,p); i11=w1*w3%p; i10=w2*w3%p;

w3=(v20-v10)%p; w4=(v21-v11)%p;
w5=i10*w3%p; w6=i11*w4%p;
w7=((i10+i11)*(w3+w4)-w5-w6)%p;
s1=(w7-u21*w6)%p; s0=(-u20*w6+w5)%p;

if s1==0
 print "s1=0!\n"
 exit();
end

w3=inv(s1,p)
u30=(w3*(w3*(s0**2+u11+u21-f4)+2*(v11-s0*w1))+w2)%p;
u31=(w3*(2*s0-w3)-w1)%p; u32=1;

w1=(u30-u10)%p; w2=(u11-u31)%p;
v30=(s1*u30*w2+s0*w1-v10)%p;
v31=(s1*(u31*w2+w1)-s0*w2-v11)%p;
@U=[u31,u30]
@V=[v31,v30]
print "divk (s^2+",u31,"*s+",u30,",",v31,"*s+",v30,")\n"

end


def g3add(u12,u11,u10,u22,u21,u20,v12,v11,v10,v22,v21,v20,f5,f4,f3,f2,f1,f0,p)

t1=(u11*u20-u10*u21)%p; t2=(u12*u20-u10*u22)%p; t3=(u20-u10)%p; t4=(u21-u11)%p; t5=(u22-u12)%p; t6=t4**2%p;
t7=t3*t4%p; t8=(u12*u21-u11*u22+t3)%p; t9=(t3**2-t1*t5)%p; t10=(t2*t5-t7)%p; r=(t8*t9+t2*(t10-t7)+t1*t6)%p;

if r==0
print "r=0!\n"
exit();
end

i2=(t5*t8-t6)%p; i1=(u22*i2-t10)%p; i0=(u21*i2-(u22*t10+t9))%p;

t1=(v10-v20)%p; t2=(v11-v21)%p; t3=(v12-v22)%p; t4=t2*i1%p; t5=t1*i0%p; t6=t3*i2%p; t7=u22*t6%p;
t8=(t4+t6+t7-(t2+t3)*(i1+i2))%p; t9=(u20+u22)%p; t10=(t9+u21)*(t8-t6)%p;
t9=(t9-u21)*(t8+t6)%p; s0_=-(u20*t8+t5)%p; s2_=(t6-(s0_+t4+(t1+t3)*(i0+i2)+(t10+t9)*inv(2,p)))%p;
s1_=(t4+t5+(t9-t10)*inv(2,p)-(t7+(t1+t2)*(i0+i1)))%p;

if s2_==0
print "s2_=0\n"
exit();
end

t1=inv(r*s2_,p); t2=r*t1%p; w=t1*s2_**2%p; wi=r*t2%p; s0=t2*s0_%p; s1=t2*s1_%p;

t6=(s0+s1)%p; t1=(u10+u12)%p; t2=t6*(t1+u11)%p; t3=(t1-u11)*(s0-s1)%p; t4=u12*s1%p;
z0=u10*s0%p; z1=((t2-t3)*inv(2,p)-t4)%p; z2=((t2+t3)*inv(2,p)-z0+u10)%p; z3=(u11+s0+t4)%p; z4=(u12+s1)%p;

ut3=(z4+s1-u22)%p; 
t1=s0*z3%p; t5=(s1*z4-u22*ut3)%p; ut2=(z3+s0+t5-u21)%p; t3=u21*ut2%p; t4=(t1-t3)%p; 
t2=(u22+u21)*(ut3+ut2)%p; ut1=(z2+t6*(z4+z3)+wi*(2*v12-wi)-(t5+t2+t4+u20))%p;
ut0=(z1+t4+s1*z2+wi*(2*(v11+s1*v12)+wi*u12)-(u22*ut1+u20*ut3))%p;

t1=(ut3-z4)%p; vt0=(w*(t1*ut0+z0)+v10)%p; vt1=(w*(t1*ut1+z1-ut0)+v11)%p;
vt2=(w*(t1*ut2+z2-ut1)+v12)%p; vt3=w*(t1*ut3+z3-ut2)%p;

t1=2*vt3%p; u32=-(ut3+vt3**2)%p; u31=(f5-(ut2+u32*ut3+t1*vt2))%p;
u30=(f4-(ut1+vt2**2+u32*ut2+u31*ut3+t1*vt1))%p;
@U=[u32,u31,u30]
print "s^3+",u32,"*s^2+",u31,"*s+",u30,"\n"

v32=(vt2-u32*vt3)%p; v31=(vt1-u31*vt3)%p; v30=(vt0-u30*vt3)%p;
@V=[v32,v31,v30]

end


def g3dbl(u12,u11,u10,v12,v11,v10,f5,f4,f3,f2,f1,f0,p)

t1=(u11*v10-u10*v11)%p; t2=(u12*v10-u10*v12)%p; t3=v11**2%p; t4=v11*v10%p; t5=(v10+u12*v11-u11*v12)%p;
t6=(v10**2-v12*t1)%p; t7=(v12*t2-t4)%p; r=(t5*t6+t2*(t7-t4)+t1*t3)%p;

if r==0
print "r=0!\n"
exit();
end

i2=(t3-v12*t5)%p; i1=(u12*i2+t7)%p; i0=(u11*i2+u12*t7+t6)%p;

t1=2*u10%p; t2=2*u11%p; t3=u12**2%p; t4=(f4-(t1+v12**2))%p; t5=(f5+t3-t2)%p; t10=2*v12%p; z2=(t5+2*t3)%p;
z1=(u12*(t2-t5)+t4)%p; z0=(f3+t3*(t5-u11)+u12*(t1-t4)+u11*(u11-f5)-t10*v11)%p;

t1=i1*z1%p; t2=i0*z0%p; t3=i2*z2%p; t4=u12*t3%p; t5=((i2+i1)*(z2+z1)-(t1+t3+t4))%p; t6=u10*t5%p;
t7=(u10+u12)%p; t8=(t7+u11)%p; t9=(t7-u11)%p; t7=t8*(t3+t5)%p; t11=t9*(t5-t3)%p;
s2_=(t1+t6+(i2+i0)*(z2+z0)-(t2+t3+(t7+t11)*inv(2,p)))%p;
s1_=(t4+(i0+i1)*(z1+z0)+(t11-t7)*inv(2,p)-(t1+t2))%p; s0_=(t2-t6)%p;

if s2_==0 
print "s2_=0!\n"
exit();
end

t1=2*r%p; t2=inv(t1*s2_,p); t3=t1*t2%p; w=t2*s2_**2%p; wi=t1*t3%p; s0=t3*s0_%p; s1=t3*s1_%p;

t1=t8*(s1+s0)%p; t2=t9*(s0-s1)%p; t3=u12*s1%p;
g0=u10*s0%p; g1=((t1-t2)*inv(2,p)-t3)%p; g2=(u10+(t1+t2)*inv(2,p)-g0)%p; g3=(t3+u11+s0)%p; g4=(u12+s1)%p;

ut3=2*s1%p; ut2=(s1**2+2*s0)%p; ut1=(ut3*s0+wi*(t10-wi))%p; ut0=(s0**2+2*wi*((s1-u12)*v12+v11+wi*u12))%p;

t1=(ut3-g4)%p; vt0=(w*(t1*ut0+g0)+v10)%p; vt1=(w*(t1*ut1+g1-ut0)+v11)%p;
vt2=(w*(t1*ut2+g2-ut1)+v12)%p; vt3=(w*(t1*ut3+g3-ut2))%p;

t1=2*vt3%p; u22=-(ut3+vt3**2)%p; u21=(f5-(ut2+u22*ut3+t1*vt2))%p;
u20=(f4-(ut1+vt2**2+u22*ut2+u21*ut3+t1*vt1))%p;
@uu=[u22,u21,u20]
print "s^3+",u22,"*u^2+",u21,"*u+",u20,"\n"

v22=(vt2-u22*vt3)%p; v21=(vt1-u21*vt3)%p; v20=(vt0-u20*vt3)%p;
@vv=[v22,v21,v20]

end

def j2add(u11, u10, v11, v10, u21, u20, v21, v20, f4, f3, f2, f1, f0, p)
  a0 = u11 * u10; a1 = u21 * u20; d0 = u10 - u20; d1 = u11 - u21
  b0 = u11 * u11; b1 = u21 * u21; c0 = v20 - v10; c1 = v21 - v11
  e0 = -u20 + u10; s1 = a1 - a0; s2 = b1 - b0; s3 = s2 + e0
  iv = inv(d0 * s3 - d1 * s1, p)
  ss2 = u20 + u11 * u21 + u10; ss3 = u21 + u11
  l3 = inv(d0 * c1 - d1 * c0, p); l2 = inv(c0 * s3 - c1 * s1, p)
  l1 = u11 * l2 + v11 + (b0 - u10) * l3
  u31 = 2 * l2 * l3 + 1 - ss3
  u30 = 2 * l1 * l3 + l2 * l2 + 1 - f4 - u31 * ss3 - ss2
  ul0 = u30 * l3; ul = -u31 * l2 + l1
  v31 = -(u31 * u31 - ul0 + ul); v30 = -(u31 * ul0 + ul)
  @U=[u31,u30]
  @V=[v31,v30]
  #return u31, u30, v31, v30
end

def j2dbl(u1, u0, v1, v0, f3, f2, f1, f0, p)
  uu1 = u1 * u1; uu0 = u1 * u0; uv01 = u0 * v1; uv10 = u1 * v0; uv11 = u1 * v1
  uv00 = u0 * v0
  d0 = 6 * v1 * uu1 - uv10; d1 = -4 * uv11 + 4 * v0; d2 = 2 * v1
  d3 = 6 * v1 * uu0 - 6 * uv00; d4 = -4 * uv01; d5 = 2 * v0
  e0 = 5 * (-u1 * uu1 + 2 * uu0) - 3 * f3 * u1 + 2 * f2
  e1 = 5 * (-u0 * uu0 + u0 * u0) - 3 * f3 * u0 + f1
  m0 = d3 - d5 * (uu1 - u0); m1 = d4 - d5 * (-u1)
  m3 = d0 - d2 * (uu1 - u0); m4 = d1 - d2 * (-u1)
  s1 = e1 - d5 * v1; s2 = e0 - d2 * v1
  iv = inv(m0 * m4 - m1 * m3, p)
  l3 = inv(m4 * s1 - m1 * s2, p); l2 = inv(m0 * s2 - m3 * s1, p)
  l1 = v1 + u1 * l2 - (uu1 - u0) * l3; l0 = v0 + u0 * l2 - (uu0) * l3
  ue1 = 2 * l3 * l2 - 2 * u1 - 1
  ue0 = 2 * l3 * l1 + l2 * l2 - 2 * u0 - uu0 - 2 * ue1 * u1
  uue1 = ue1 * ue1; uue0 = ue1 * ue0
  ve1 = (uu1 - ue0) * l3 - ue1 * l2 + l1
  ve0 = uue0 * l3 - ue0 * l2 + l0
  @uu=[ue1,ue0];
  @vv=[ve1,ve0];
  #return ue1, ue0, ve1, ve0
end

def g4add(a1,a2,b1,b2)

r23=e+a; r22=f+b; r21+g+c; r20=h+d; t20=1; r33=r23*a+r22;
r32=r23*b+r21; r31=r23*c+r20; r30=r23*d; t31=1; t30=r23;
r42=r33*r22+r23*r32; r41=r33*r21+r23*r31; r40=r33*r20+r23*r30;
t41=r23; t40=r33+r23**2; r52=r42+r33*t40; t50=r42*r23;
r50=r42*r30; t52=r33*t41; t51=r42+r33*t40; t50=r42*r23;
r61=r52*r41+r42*r51; r60=r52*r40+r42*r50; t62=r42*t52;
t61=r52*t41+r42*t51; t60=r52*t40+r42*t50; r71=r61*r51+r52*r60;
r70=r61*r50; t73=r52*t62; t72=r61*t52+r52*t61; t71=r61*t51+r52*t60;
t70=r61*t50;rt70 r80=r71*r60+r61;
inv1=r71*t61+r61+t61*t71; inv0=r71=r71*t60+r61*t70;

ta=m+i; tb=n+j; tc=o+k; td=p+l; te=inv3; tf=inv2; tg=inv1;
th=inv0; t0=te*tg; t1=tb*tf; t2=ta*te; t3=tb*tg; t4=tc*tf; t10=td*th
t11=(tc+td)*(tg+th)+t0+t10; t12=(tb+t-d)*(tf+th)+t10+t1+t0;
t13=(ta+td)*(te+th)+t10+t2+t3+t4; t14=(ta+tc)*(te+tg)+t2+t0+t1;
t15=(ta+tb)*(te+tf)+t2+t1; t16=t2; t17=t15+e*t16; t18=e*t17+t16*f+t14;
s3_=e*t18+f*t17+g*t16+t13; s2_=f*t18+g*t17+h*t16+t12;
s1_=g*t18+h*t17+t11; s0_=h*t18+t10;

t1=r80*s3_; w6=inv(t1,@p3); w7=r80*w6; w4=r80*w7; w3=s3_**2*w6;
w5=w4**2; s0=s0_*w7; s1=s1_*w7; s2=s2_*w7;

t0=c*s1; t1=b*s2; z0=s0*d; z1=(c+d)*(s1+s0)+z0+t0; 
z2=(b+d)*(s2+s0)+z0+t1+t0; z3=(a+d)*(1+s0)+z0+a+b*s1+c*s2;
z4=(a+c)*(1+s1)+a+t0+t1+s0; z5=(a+b)*(1+s2)+a+t1+s1; z6=a+s2;
z7=1;

t1=s2*w4; t2=s1*w4; diff4=s2*z6+z5+s1+f; diff3=g+z4+s0+s2*z5+s1*z6;
diff2=h+z3+s2*z4+s1*z5+s0*z4+w5*a+z1; diff1=s2*z3+s1*z4+s0*z5+z2+w5;
diff0=w4+s2*z2+s1*z3+s0*z4+w5*a+z1; u5_=z6+s2+e; u4_=diff4+e*u5_;
u3_=diff3+e*u4_+f*u5_; u2_=diff2+e*u3_+f*u4_+g*u5_;
u1_=diff1+e*u2_+f*u3_+g*u4_+h*u5_; u0_=diff0+e*u1_+f*u2_+g*u3_+h*u4_;

t1=u5_+z6; v5_=w3*(u5_*t1+u4_+z5); v4_=w3*(u4_*t1+u3_+z4);
v3_=w3*(u3_*t1+u2_+z3)+i; v2_=w3*(u2_*t1+u1_+z2)+j;
v1_=w3*(u1_*t1+u0_+z1)+1+k; v0_=w3*(u0_*t1+z0)+l;

diff3=1; diff2=v4_**2; diff1=f7; diff0=f6+v5_+v3_**2; u34=v5_**2;
u33=diff3+u34*u5_; u32=diff2+u34*u4_+u33*u5_;
u31=diff1+u34*u3_+u33*u4_+u32*u5_; 
u30=diff0+u34*u2_+u33*u3_+u32*u4_+u31*u5_;

t0=inv(u34,@p3); agoal=u33*t0; bgoal=u32*t0; cgoal=u31*t0; dgoal=u30*t0;
@U=[agoal,bgoal,cgoal,dgola]

t0=v4_+v5_*a3; igoal=agoal*t0+bgoal*v5_+v3_; jgoal=bgoal*t0+cgoal*v5_+v2_;
kgoal=cgoal*t0+dgoal*v5_+v1_+1; lgoal=d*t0+v0_;
@V=[igoal,jgoal,kgoal,lgoal]

end


def g4dbl(a,b,c,d,i,j,k,l,f7,f6,f5,f4,f3,f2,f1,f0,p)

inv3=1; inv2=a; inv1=b; inv0=c;

diff2=(i*i+c)%p; diff0=(j*j+i)%p; z4_=a; z3_=(b+a*a+f7)%p; z2_=(diff2+a*(z3_+b)+f6)%p;
z1_=(d+a*(z2_+c)+b*z3_+f5)%p; z0_=(diff0+a*(z1_+d)+b*z2_+c*z3_+f4)%p;
z3=(b+z3_)%p; z2=(c+z2_)%p; z1=(d+z1_)%p; z0=z0_;

t1=(z1*inv1)%p; t2=(z2*inv2)%p; t10=(z0*inv0)%p; t16=(z3*inv3)%p;
t15=((z3+z2)*(inv2+inv3)+t16+t2)%p; t14=((z3+z1)*(inv3+inv1)+t16+t1+t2)%p;
t13=((z3+z0)*(inv3+inv0)+t10+t16+z2*inv1+z1*inv2)%p; 
t12=((z2+z0)*(inv2+inv0)+t10+t2+t1)%p; t11=((z1+z0)*(inv0+inv1)+t1+t10)%p;
t3=(t15+a*t16)%p; t4=(a*t3+b*t16+t14)%p; s0_=(d*t4+t10)%p; s1_=(c*t4+d*t3+t11)%p;
s2_=(b*t4+c*t3+d*t16+t12)%p; s3_=(a*t4+b*t3+c*t16+t13)%p;

t1=d*s3_%p; w6=inv(t1,p); w7=d*w6%p; w4=d*w7%p; w3=s3_*s3_*w6%p;
w5=w4*w4%p; s0=s0_*w7%p; s1=s1_*w7%p; s2=s2_*w7%p;

t0=c*s1%p; t1=b*s2%p; g0=s0*d%p; g1=((c+d)*(s1+s0)+t0+g0)%p;
g2=((b+d)*(s2+s0)+g0+t1+t0)%p; g3=((a+d)*(1+s0)+g0+a+b*s1+c*s2)%p;
g4=((a+c)*(1+s1)+a+t0+t1+s0)%p; g5=((a+b)*(1+s2)+a+t1+s1)%p; g6=(a+s2)%p;

t0=a*a%p; t1=b*b%p; t2=c*c%p; t3=g6*g6%p; t4=g5*g5%p; u5_=0; u4_=(t3+t0)%p; u3_=0;
u2_=(t4+t1+t0*(t3+t0))%p; u1_=w5; u0_=(g4*g4+w4+t2+t0*u2_+t1*(t0+t3))%p; 

v5_=w3*(u4_+g5)%p; v4_=w3*(u4_*g6+g4)%p; v3_=(w3*(u2_+g3)+i)%p;
v2_=(w3*(u2_*g6+u1_+g2)+j)%p; v1_=(w3*(u1_*g6+u0_+g1)+1+k)%p;
v0_=(w3*(u0_*g6+g0)+l)%p;

diff3=1; diff2=v4_*v4_%p; diff1=f7; diff0=(f6+v5_+v3_*v3_)%p; u24=v5_*v5_%p;
u23=diff3; u22=(diff2+u24*u4_)%p; u21=(diff1+u23*u4_)%p;
u20=(diff0+u24*u2_+u22*u4_)%p;

t0=inv(u24,p); agoal=u23*t0%p; bgoal=u22*t0%p; cgoal=u21*t0%p; dgoal=u20*t0%p;
@uu=[agoal,bgoal,cgoal,dgoal]
print "s^4+",agoal,"s^3+",bgoal,"s^2+",cgoal,"s+",dgoal,"\n"

t0=(v4_+v5_*agoal)%p; igoal=(agoal*t0+bgoal*v5_+v3_)%p; jgoal=(bgoal*t0+cgoal*v5_+v2_)%p;
kgoal=(cgoal*t0+dgoal*v5_+v1_+1)%p; lgoal=(dgoal*t0+v0_)%p;
@vv=[igoal,jgoal,kgoal,lgoal]

end


def mktable(u,v)
# print @q ,"\n"
 pp= 31
 print  pp ,"\n"
 @uu= u
 print  @uu ,"\n"
 @vv= v
 print  @vv ,"\n"
 a= [1,1,1,1,1]
# print  a ,"\n"
# print  @CRV_b ,"\n"


 # enzan table
   @le_u[0]=@uu
   @le_v[0]=@vv
print @uu," l673\n"
print @vv," l674\n"

for i in 1..N-1 #begin Pub_key at plain

   #chao(@uu[0],@uu[1],@vv[0],@vv[1],@FF[0],@FF[1],@FF[2],@P)
   #j2dbl(@uu[0],@uu[1],@vv[0],@vv[1],@FF[1],@FF[2],@FF[3],@FF[4],@P)
   #g4dbl(@uu[0],@uu[1],@uu[2],@uu[3],@vv[0],@vv[1],@vv[2],@vv[3],0,0,0,0,0,0,29,0,@g4p)
    g3dbl(@uu[0],@uu[1],@uu[2],@vv[0],@vv[1],@vv[2],0,0,0,0,0,157788440085147827,@p3)
   @le_u[i]=@uu
   @le_v[i]=@vv
#print i,"=",@le_u[i],"\n"

#   @uu=@le_u[i]
#   @vv=@le_v[i]
 end #of for
print i,"i\n"

end



#D'=mD
def jac(kk,p)

j=0
ki=[N]

   for j in 0..N
     ki[j]=0
   end
  j=0
   for i in 0..N
     if(((kk^(1<<i))>>i)%2 == 0) #testbit(kk,i)
       ki[j]=i
	j=j+1
	print i,"L704\n" 
     end
   end

   @U=@le_u[ki[0]]
   @V=@le_v[ki[0]]

#print @le_u[ki[0]],"\n"
#print @le_u[ki[2]],"\n"
#exit();
for i in 1..j-1
print j,",",i,"\n"
 if(@U != @le_u[ki[i]])
    #j2add(@U[0],@U[1],@le_u[ki[i]][0],@le_u[ki[i]][1],@V[0],@V[1],@le_v[ki[i]][0],@le_v[ki[i]][1],@FF[0],@FF[1],@FF[2],@FF[3],@FF[4],p)
#print "i=",i," ",@le_u[ki[i]][0],"\n"
    #kotehan(@U[0],@U[1],@le_u[ki[i]][0],@le_u[ki[i]][1],@V[0],@V[1],@le_v[ki[i]][0],@le_v[ki[i]][1],@FF[0],p)
	g3add(@U[0],@U[1],@U[2],@le_u[ki[i]][0],@le_u[ki[i]][1],@le_u[ki[i]][2],@V[0],@V[1],@V[2],@le_v[ki[i]][0],@le_v[ki[i]][1],@le_v[ki[i]][2],0,0,0,0,0,157788440085147827,p)
end

 if(@U==@le_u[ki[i]])
  print @le_u[ki[i]],"=",@U,"\n"
  print i,"de('A`)\n"
  exit();
end
end #of for

end


def zeta()



end

#/* frobenius map */
def fro()


p=10**17+3
x=41899918437810026
y=68564906982200716
print (x**5+@b[0]*x**4+@b[1]*x**3+@b[2]*x*x+@b[3]*x+@b[4])%p,"\n"
print y*y%p,"\n"

x=73271434588187944
y=91647395753297815
print (x**5+@b[0]*x**4+@b[1]*x**3+@b[2]*x*x+@b[3]*x+@b[4])%p,"\n"
print y*y%p,"\n"
x=79038146144619618
y=78852114870530227
print (x**5+@b[0]*x**4+@b[1]*x**3+@b[2]*x*x+@b[3]*x+@b[4])%p,"\n"
print y*y%p,"\n"

x=90305483630909247
y=98869188447477807
print (x**5+@b[0]*x**4+@b[1]*x**3+@b[2]*x*x+@b[3]*x+@b[4])%p,"\n"
print y*y%p,"\n"

q=10**19+51
aa=(inv(4630052179246799403,q)*589402353553582342)%q
print "a=",aa,"\n"
print inv(4630052179246799403,q)*4630052179246799403%q,"\n"
print "b=",(2935061693073737419-aa*5239897978117534135)%q,"\n"
print "b'=",(3524464046627319761-aa*9869950157364333538)%q,"\n"

aa=(616419646419685014*inv(7058278367076753082,q))%q
print "a=",aa,"\n"
print "b=",(616419646419685014-aa*3542790122851877922)%q,"\n"
print "b'=",(-1*aa*6484511755775124891)%q,"\n"

#504894935863953268767725=a*106028591185649525291891+b
#704210062398295465154981=a*1487990384692386499004424+b
#199315126534342196387256=a*1381961793506736973712533
aa=199315126534342196387256*inv(1381961793506736973712533,@p3)%@p3
print "a=",aa,"\n"
print "b=",(704210062398295465154981-aa*1487990384692386499004424)%@p3,"\n"
print "b'=",(504894935863953268767725-aa*106028591185649525291891)%@p3,"\n"

#1=a*704665787761008893641614+b
#824513992484349685277891=a*1086395356035375672830976+b
#824513992484349685277890=a*381729568274366779189362
aa=824513992484349685277890*inv(381729568274366779189362,@p3)%@p3
print "a=",aa,"\n"
print "b=",(1-aa*704665787761008893641614)%@p3,"\n"
print "b'=",(824513992484349685277891-aa*1086395356035375672830976)%@p3,"\n"

#1090190095529845640563737560=a*901573235033529767913345809+b
#1607209110223949051233778495=a*2189334496065905869799867124+b
#517019014694103410670040935=a*1287761261032376101886521315
aa=517019014694103410670040935*inv(1287761261032376101886521315,@p4)%@p4
print "a=",aa,"\n"
print "b=",(1090190095529845640563737560-aa*901573235033529767913345809)%@p4,"\n"
print "b'=",(1607209110223949051233778495-aa*2189334496065905869799867124)%@p4,"\n"

#1044171804289905858226438503=a*10486172923382208811538764+b
#928996521279782754969642877=a*1874123356916514825274712274+b
#262950140661994395088022052=a*1863637183993132616463173510
aa=262950140661994395088022052*inv(1863637183993132616463173510,@p4)%@p4
print "a=",aa,"\n"
print "b=",(928996521279782754969642877-aa*1874123356916514825274712274)%@p4,"\n"
print "b'=",(1044171804289905858226438503-aa*10486172923382208811538764)%@p4,"\n"



p=5000000000000000008503491
#2408698334445693709646495=a*2779199501981512279739817+b
#1802105445885051028400128=a*50724386855111482309402+b
#606592888560642681246367=a*2728475115126400797430415
aa=606592888560642681246367*inv(2728475115126400797430415,p)%p
print "a=",aa,"\n"
print "b=",2408698334445693709646495-aa*2779199501981512279739817,"\n"
print "b'=",1802105445885051028400128-aa*50724386855111482309402,"\n"

#1=a*2055622596816515886446193+b
#3=a*1553122609714208136553134+b
#-2=a*502499987102307749893059
aa=-2*inv(502499987102307749893059,p)%p
print "a=",aa,"\n"
print "b=",1-aa*2055622596816515886446193,"\n"
print "b'=",3-aa*1553122609714208136553134,"\n"

end



HEC()



mktable(@wu,@wv)
#exit()
jac(@J3,@p3)
#jac(@Jga,@P)
print "debug\n"
exit()
jac(6,@g4p)
print "debug\n"
jac(100,@g4p)
print "debug\n"
jac(@Jta,@g4p)

#jac(3,@p3)
