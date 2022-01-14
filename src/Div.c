Div g2dbl(Div D,OP ff){
__int128_t  u1,u0,v1,v0,vu,vv,mm1,mm2,mm3,mm4,z1,z2,t1,t2,t3,t4,l2,l1,l3,d,aa,bb,cc,dd,ee,ue1,ue0,uue1,ve1,ve0,uue0,uu1,uu0;
Div D2={0};
vec vx;
unsigned f2,f3;
OP u,v;

u=D.u;
v=D.v;
vx=o2v(u);
u1=vx.x[1];
u0=vx.x[0];
vx=o2v(v);
v1=vx.x[1];
v0=vx.x[0];

vx=o2v(ff);
f2=vx.x[2];
f3=vx.x[3];

uu1=(u1*u1)%P; uu0=(u0*u1)%P;

vv=(v1*v1)%P; vu=(((v1+u1)*(v1+u1))-vv-uu1)%P; mm1=(2*v0-2*vu)%P; mm2=(2*v1*(u0+2*uu1))%P;
mm3=(-2*v1)%P; mm4=(vu+2*v0)%P; z1=(f2+2*uu1*u1+2*uu0-vv)%P; z2=(f3-2*u0+3*uu1)%P;
t1=((mm2-z1)*(z2-mm1))%P; t2=((-z1-mm2)*(z2+mm1))%P; 
t3=((mm4-z1)*(z2-mm3))%P; t4=((-z1-mm4)*(z2+mm3))%P;
l2=(t1-t2)%P; l3=(t3-t4)%P; d=(t3+t4-t1-t2-2*(mm2-mm4)*(mm1+mm3))%P;
aa=inv(d*l3,P); bb=(d*aa)%P; cc=(d*bb)%P; cc=(d*bb)%P; dd=(l2*bb)%P; ee=(l3*l3*aa)%P;
ue1=(2*dd-cc*cc-2*u1)%P; ve0=((dd-u1)*(dd-u1)+2*cc*(v1+cc*u1))%P; uue1=(ue1*ue1)%P; uue0=(ue1*ue0)%P;
ve1=(ee*ve1+v1)%P; ve0=(ee*ve0+v0)%P;
memset(vx.x,0,sizeof(vx.x));
vx.x[2]=1;
vx.x[1]=ue1;
vx.x[0]=ue0;
D2.u=v2o(vx);

memset(vx.x,0,sizeof(vx.x));
vx.x[1]=ve1;
vx.x[0]=ve0;
D2.v=v2o(vx);


return D2;
}

Div cdbl(Div D,OP f){
OP k,s,l,u,v,u1,u2,v2,ut,vt,vc;
Div D2;
EX F;
OP ud,vd;

u=D.u;
v=D.v;
if(chkdiv(D,f)==-1){
  printf("bakayo\n");
  exit(1);
}

OP t={0};
F=xgcd(u,scr(2,v));
//F=monic(F);
printpoln(o2v(F.u));
printpoln(o2v(F.v));
//exit(1);
u2=odiv(omul(u,u),omul(F.d,F.d));
v2=omod(scr(inv(LT(omul(F.d,F.d)).a,P),oadd(omul(F.v,omul(u,v)),omul(F.u,oadd(omul(v,v),f)))),u2);
printpoln(o2v(u2));
printpoln(o2v(v2));
//exit(1);

while(odeg(u2)>2){
ud=odiv(osub(f,omul(v2,v2)),u2);
vd=omod(minus(v2),ud);
printf("7777777777777777777777777777777777\n");
u2=ud;
v2=vd;
}

ud=u2;
vd=v2;
ud=monique(ud);


/*
k=odiv(osub(f,omul(v,v)),u);
s=omod(omul(k,qinv(scr(2,v),u)),u);
l=omul(s,u);
u1=osub(omul(s,s),omod(osub(omul(scr(2,v),s),k),u));
u2=monique(u1);
v2=omod(minus(oadd(l,v)),u2);
*/
printpoln(o2v(F.d));
printpoln(o2v(ud));
printpoln(o2v(vd));
//exit(1);
D2.u=ud;
D2.v=vd;

if(chkdiv(D2,f)==-1){
printf("dame\n");
exit(1);
}

printf("ii!\n");
return D2;
}
Div dobule(Div D, OP ff){
vec vx=o2v(ff);
__int128_t  d0,uv010,uu1,uu0,uuv11,d1,d2,d3,d4,d5,f2,e0,e1,f1,f3,mm0,mm1,mm3,mm4,s1,s2,l1,l0,l2,l3,u31,u30,v31,v30,uue0,uue1,iv,uuv01,uuv10,uuv00,u1,u0,v1,v0;
OP u,v;
Div D2;
f3=vx.x[3];
f2=vx.x[2];
f1=vx.x[1];

u=D.u;
v=D.v;
memset(vx.x,0,sizeof(vx.x));
vx=o2v(u);
u1=vx.x[1];
u0=vx.x[0];

memset(vx.x,0,sizeof(vx.x));
vx=o2v(v);
v1=vx.x[1];
v0=vx.x[0];

uu1=(u1*u1)%P; uu0=(u1*u0)%P; uuv01=(u0*v1)%P; uuv10=(u1*v0)%P; uuv11=(u1*v1)%P;
uuv00=(u0*v0)%P;

  d0=(6*v1*uu1-uuv10)%P; d1= (-4*uuv11+4*v0)%P; d2= (2*v1)%P;
  d3=(6*v1*uu0-6*uuv00)%P; d4= (-4*uuv01)%P; d5=(2*v0)%P;
  e0=(5*(-u1*uu1+2*uu0)-3*f3*u1+2*f2)%P;
  e1=(5*(-u0+uu0+u0*u0)-3*f3*u0+f1)%P;
  mm0=(d3-d5*(uu1-u0))%P; mm1=(d4-d5*(-u1))%P;
  mm3=(d0-d2*(uu1-u0))%P; mm4=(d1-d2*(-u1))%P;
  s1=(e1-d5*v1)%P; s2=(e0-d2*v1)%P;
  iv=inv((mm0*mm4-mm1*mm3),P);
  l3= (iv*(mm4*s1-mm1*s2))%P; l2=(iv*(mm0*s2-mm3*s1))%P;
  l1=(v1+u1*l2-(uu1-u0)*l3)%P; l0=(v0+u0*l2-(uu0)*l3)%P;
  u31=(2*l3*l2-2*u1-1)%P;
  u30=(2*l3*l1+l2*l2-2*u0-uu1-2*u31*u1)%P;
  uue1=(u31*u31)%P; uue0=(u31*u30)%P;
  v31=((uu1-u30)*l3-u31*l2+l1)%P;
  v30=(uue0*l3-u30*l2+l0)%P;

memset(vx.x,0,sizeof(vx.x));
vx.x[1]=u31;
vx.x[0]=u30;
vx.x[2]=1;
u=v2o(vx);
memset(vx.x,0,sizeof(vx.x));
vx.x[1]=v31;
vx.x[0]=v30;
vx.x[2]=1;
v=v2o(vx);

D2.u=u;
D2.v=v;

return D2;
}

Div g2add(OP ff, OP uu1, OP uu2, OP vv1, OP vv2)
{
  OP ll, u;
  OP v, s, l, k, v3, u3;
  Div X;

  ll = (omul(vv2, vv2));
  printpol(o2v(ll));
  printf("\n");
  ll = (osub(ff, (ll)));
  printpol(o2v(ll));
  printf("\n");
  k = (odiv(ll, uu2));
  printpol(o2v(k));
  printf(" ('A`)\n");
  //exit(1);
  ll=qinv(uu2,uu1);
  s=omod(omul(osub(vv1,vv2),ll),uu1);
  printpol(o2v(ll));
  printf("========inv\n");
l=omul(s,uu2);
u=odiv(osub(k,omul(s,oadd(l,omul(s,vv2)))),uu1);
printpol(o2v(u));
printf(" ===u\n");
u3=monique(u);
printpol(o2v(u3));
printf(" ===u3\n");
v3=omod(minus(oadd(l,vv2)),u3);
printpol(o2v(v3));
printf(" ===v3\n");
exit(1);


  X.u = u3;
  X.v = v3;

  return X;
}

