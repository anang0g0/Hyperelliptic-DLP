p3=631153760340591307n
J3=251423300188980936808533376314868064530443303970434811n
wu=[377166208405888276n, 65838350729795034n, 258349841662486010n]
wv=[82817318614465843n, 235991131611388836n, 245978366599867468n]

N=256n
le_u=[N,N,N,N]
le_v=[N,N,N,N]
uu=[N,N,N,N]
vv=[N,N,N,N]
U=[N,N,N,N]
V=[N,N,N,N]

//# invert of integer
function inv(a, n){
  
    d = n
    x = 0n
    s = 1n
    while (a != 0n){
      q = BigInt(d) / BigInt(a)
      r = BigInt(d) % BigInt(a)
      d = a
      a = r
      t = x - q * s
      x = s
      s = t
    }
    gcd = d  
  
    return ((x + n) % (n / d))
  }
  
function g3add(u12,u11,u10,u22,u21,u20,v12,v11,v10,v22,v21,v20,f5,f4,f3,f2,f1,f0,p){

    t1=(u11*u20-u10*u21)%p; t2=(u12*u20-u10*u22)%p; t3=(u20-u10)%p; t4=(u21-u11)%p; t5=(u22-u12)%p; t6=t4*t4%p;
    t7=t3*t4%p; t8=(u12*u21-u11*u22+t3)%p; t9=(t3*t3-t1*t5)%p; t10=(t2*t5-t7)%p; r=(t8*t9+t2*(t10-t7)+t1*t6)%p;
    
    if (r==0n){
    console.log("r=0!<br>")
    exit();
    }
    
    i2=(t5*t8-t6)%p; i1=(u22*i2-t10)%p; i0=(u21*i2-(u22*t10+t9))%p;
    
    t1=(v10-v20)%p; t2=(v11-v21)%p; t3=(v12-v22)%p; t4=t2*i1%p; t5=t1*i0%p; t6=t3*i2%p; t7=u22*t6%p;
    t8=(t4+t6+t7-(t2+t3)*(i1+i2))%p; t9=(u20+u22)%p; t10=(t9+u21)*(t8-t6)%p;
    t9=(t9-u21)*(t8+t6)%p; s0_=-(u20*t8+t5)%p; s2_=(t6-(s0_+t4+(t1+t3)*(i0+i2)+(t10+t9)*inv(2n,p)))%p;
    s1_=(t4+t5+(t9-t10)*inv(2n,p)-(t7+(t1+t2)*(i0+i1)))%p;
    
    if(s2_==0n){
    console.log("s2_=0\n")
    exit();
  }
    
    t1=inv(r*s2_,p); t2=r*t1%p; w=t1*s2_*s2_%p; wi=r*t2%p; s0=t2*s0_%p; s1=t2*s1_%p;
    
    t6=(s0+s1)%p; t1=(u10+u12)%p; t2=t6*(t1+u11)%p; t3=(t1-u11)*(s0-s1)%p; t4=u12*s1%p;
    z0=u10*s0%p; z1=((t2-t3)*inv(2,p)-t4)%p; z2=((t2+t3)*inv(2,p)-z0+u10)%p; z3=(u11+s0+t4)%p; z4=(u12+s1)%p;
    
    ut3=(z4+s1-u22)%p; 
    t1=s0*z3%p; t5=(s1*z4-u22*ut3)%p; ut2=(z3+s0+t5-u21)%p; t3=u21*ut2%p; t4=(t1-t3)%p; 
    t2=(u22+u21)*(ut3+ut2)%p; ut1=(z2+t6*(z4+z3)+wi*(2*v12-wi)-(t5+t2+t4+u20))%p;
    ut0=(z1+t4+s1*z2+wi*(2*(v11+s1*v12)+wi*u12)-(u22*ut1+u20*ut3))%p;
    
    t1=(ut3-z4)%p; vt0=(w*(t1*ut0+z0)+v10)%p; vt1=(w*(t1*ut1+z1-ut0)+v11)%p;
    vt2=(w*(t1*ut2+z2-ut1)+v12)%p; vt3=w*(t1*ut3+z3-ut2)%p;
    
    t1=2*vt3%p; u32=-(ut3+vt3*vt3)%p; u31=(f5-(ut2+u32*ut3+t1*vt2))%p;
    u30=(f4-(ut1+vt2*vt2+u32*ut2+u31*ut3+t1*vt1))%p;
    U=[u32,u31,u30]
    console.log( "s^3+",u32,"*s^2+",u31,"*s+",u30,"\n")
    
    v32=(vt2-u32*vt3)%p; v31=(vt1-u31*vt3)%p; v30=(vt0-u30*vt3)%p;
    V=[v32,v31,v30]
    
}
    
    
    function g3dbl(u12,u11,u10,v12,v11,v10,f5,f4,f3,f2,f1,f0,p){
    
    t1=(u11*v10-u10*v11)%p; t2=(u12*v10-u10*v12)%p; t3=(v11*v11)%p; t4=v11*v10%p; t5=(v10+u12*v11-u11*v12)%p;
    t6=(v10*v10-v12*t1)%p; t7=(v12*t2-t4)%p; r=(t5*t6+t2*(t7-t4)+t1*t3)%p;
    
    if (r==0){
    console.log("r=0!\n")
    exit();
    }
    
    i2=(t3-v12*t5)%p; i1=(u12*i2+t7)%p; i0=(u11*i2+u12*t7+t6)%p;
    
    t1=2n*u10%p; t2=2n*u11%p; t3=u12*u12%p; t4=(BigInt(f4)-BigInt(t1+v12*v12))%p; t5=(BigInt(f5)+BigInt(t3)-BigInt(t2))%p; t10=2n*v12%p; z2=(t5+2n*t3)%p;
    z1=(u12*(t2-t5)+t4)%p; z0=(BigInt(f3)+t3*(t5-u11)+u12*(t1-t4)+(u11*(BigInt(u11)-BigInt(f5)))-BigInt(t10*v11))%p;
    
    t1=i1*z1%p; t2=i0*z0%p; t3=i2*z2%p; t4=u12*t3%p; t5=((i2+i1)*(z2+z1)-(t1+t3+t4))%p; t6=u10*t5%p;
    t7=(u10+u12)%p; t8=(t7+u11)%p; t9=(t7-u11)%p; t7=t8*(t3+t5)%p; t11=t9*(t5-t3)%p;
    s2_=(t1+t6+(i2+i0)*(z2+z0)-(t2+t3+(t7+t11)*inv(2,p)))%p;
    s1_=(t4+(i0+i1)*(z1+z0)+(t11-t7)*inv(2,p)-(t1+t2))%p; s0_=(t2-t6)%p;
    
    if(s2_==0n){
    console.log("s2_=0!\n")
    exit();
  }
    
    t1=2n*r%p; t2=inv(t1*s2_,p); t3=t1*t2%p; w=t2*s2_*s2_%p; wi=t1*t3%p; s0=t3*s0_%p; s1=t3*s1_%p;
    
    t1=t8*(s1+s0)%p; t2=t9*(s0-s1)%p; t3=u12*s1%p;
    g0=u10*s0%p; g1=((t1-t2)*inv(2n,p)-t3)%p; g2=(u10+(t1+t2)*inv(2,p)-g0)%p; g3=(t3+u11+s0)%p; g4=(u12+s1)%p;
    
    ut3=2n*s1%p; ut2=(s1*s1+2n*s0)%p; ut1=(ut3*s0+wi*(t10-wi))%p; ut0=(s0*s0+2n*wi*((s1-u12)*v12+v11+wi*u12))%p;
    
    t1=(ut3-g4)%p; vt0=(w*(t1*ut0+g0)+v10)%p; vt1=(w*(t1*ut1+g1-ut0)+v11)%p;
    vt2=(w*(t1*ut2+g2-ut1)+v12)%p; vt3=(w*(t1*ut3+g3-ut2))%p;
    
    t1=2n*vt3%p; u22=-(ut3+vt3*vt3)%p; u21=(BigInt(f5)-(ut2+u22*ut3+t1*vt2))%p;
    u20=(BigInt(f4)-(ut1+vt2*vt2+u22*ut2+u21*ut3+t1*vt1))%p;
    uu=[u22,u21,u20]
    console.log("s^3+",u22,"*u^2+",u21,"*u+",u20,"\n")
    
    v22=(vt2-u22*vt3)%p; v21=(vt1-u21*vt3)%p; v20=(vt0-u20*vt3)%p;
    vv=[v22,v21,v20]
    
}

    function mktable(u,v){
         uu= u
         //print  uu ,"\n"
         vv= v
         //print  vv ,"\n"
         //a= [1,1,1,1,1]
                
         //# enzan table
           le_u[0]=uu
           le_v[0]=vv
        //print uu," l673\n"
        //print vv," l674\n"
        
        for (i=1 ;i<N;i++){// #begin Pub_key at plain
        
         g3dbl(uu[0],uu[1],uu[2],vv[0],vv[1],vv[2],0,0,0,0,0,157788440085147827,p3)
           le_u[i]=uu
           le_v[i]=vv
        //#print i,"=",le_u[i],"\n"
        
         } //#of for
        //print i,"i\n"
        
    }
        
        
        
        //D'=mD
        function jac(kk,p){
        
        j=0n
        ki=[N]
        
           for(j = 0n;j<N+1n;j++){
             ki[j]=0n
           }
          j=0n
           for(i= 0n;i<N+1n;i++){
             if(((kk^(1n<<i))>>i)%2n == 0n){ //#testbit(kk,i)
               ki[j]=i
            j=j+1n
            //print i,"L704\n" 
           }
          }
        
           U=le_u[ki[0n]]
           V=le_v[ki[0n]]
        
        for(i = 1n;i<j;i++){
         if(U != le_u[ki[i]]){
            g3add(U[0n],U[1n],U[2n],le_u[ki[i]][0n],le_u[ki[i]][1n],le_u[ki[i]][2n],V[0],V[1n],V[2n],le_v[ki[i]][0n],le_v[ki[i]][1n],le_v[ki[i]][2n],0n,0n,0n,0n,0n,157788440085147827n,p)
         }
        
         if(U==le_u[ki[i]]){
          console.log(le_u[ki[i]],"=",U,"\n")
          console.log(i,"de('A`)\n")
          exit();
         }
         }// #of for
        
        }

mktable(wu,wv)
jac(J3,p3)
        