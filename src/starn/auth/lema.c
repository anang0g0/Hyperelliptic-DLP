///////////////////////////////////////////////////////////////////////////
//
// (lema-scheme)
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <stdint.h>


#define N 64 //order of GF(q)


 void rp(unsigned char* a) {
        int i, j, x;


        for (i = 0; i < N; i++) {
            a[i] = i;
        }
        for (i = 0; i < N - 2; i++) {
            // rand from i+1 to F-1
            j = (rand() % (N - 1 - i)) + i + 1;

            // swap a[i] and a[j]
            x = a[j];
            a[j] = a[i];
            a[i] = x;
        }
        if (a[N - 1] == N - 1) {
            a[N - 1] = a[N - 2];
            a[N - 2] = N - 1;
        }

 }



int main(){
  
    int i,j;
    unsigned char t[N]={0},t2[N]={0}; //public
    unsigned char u[N]; //public
    unsigned char ss[N]={0},r[N]={0},z[N]={0},e[N],e2[N],xx[N]={0};
    unsigned char v[N]={0};
    uint64_t l[N]={0},n[N]={0};
    time_t tt;

    srand(clock()+time(&tt));

    //頂点の値ssはすべて異なるものとする
    for(i=0;i<N;i++){
      e[i]=ss[i]=rand()%2; //secret
      e2[i]=xx[i]=rand()%2;
      printf("%d,%d %d,%d\n",e[i],ss[i],e2[i],xx[i]);
    }


      rp(u); //public

for(i=0;i<N;i++)
v[i]=e[i]^e2[i];


  for(j=0;j<N;j++){

    for(i=0;i<N;i++)
      t[i]=(ss[i]^ss[(i+1)%N]);
    for(i=0;i<N;i++)
      r[i]=ss[u[i]];
    memcpy(ss,r,sizeof(ss));
printf("\n");


    for(i=0;i<N;i++)
      t2[i]=(xx[i]^xx[(i+1)%N]);
    for(i=0;i<N;i++){ 
      r[i]=0;
      r[i]=xx[u[i]];
    }
      memcpy(xx,r,sizeof(xx));
printf("\n");


    for(i=0;i<N;i++)
      z[i]=(v[i]^v[(i+1)%N]);
    for(i=0;i<N;i++){
      r[i]=0;
      r[i]=v[u[i]];
    }
      memcpy(v,r,sizeof(v));
printf("\n\n");

/*
    for(i=0;i<N;i++)
      z[j]+=(e[i]^e[(i+1)%N])+(e2[i]^e2[(i+1)%N]);  
    

    for(i=0;i<N;i++)
      r[i]=e[u[i]];
      memcpy(e,r,sizeof(v));

    for(i=0;i<N;i++)
      r[i]=e2[u[i]];
      memcpy(e2,r,sizeof(v));
*/

for(i=0;i<N;i++){
  //printf("%d,%d\n",t[i]^t2[i],z[i]);
  l[j]+=t[i]^t2[i];
  n[j]+=z[i];
  l[j]=(l[j]<<1);
  n[j]=(n[j]<<1);
    printf("%llu,%llu\n",l[i],n[i]);
  }
  
l[j]=(l[j]>>1);
n[j]=(n[j]>>1);
/*
ll+=l[j];
nn+=n[j];
*/
  }
printf("\n\n");
for(i=0;i<N;i++)
printf("%llu,%llu\n",l[i],n[i]);
printf("\n\n");


//printf("%u,%u\n",ll,nn);
//printf("\n\n");



return 0;
}