#include <stdio.h>

int main(void){

int a1,b1,c1,a2,b2,c2; //方程式の値を格納する変数
unsigned x, y; //方程式の解

printf("一つ目の方程式ax + by = c :");
scanf("%d %d %d",&a1,&b1,&c1);

printf("二つ目の方程式ax + by = c :");
scanf("%d %d %d",&a2,&b2,&c2);
putchar('\n');

printf("%dx + %dy = %d\n",a1,b1,c1); //方程式を出力
printf("%dx + %dy = %d\n",a2,b2,c2);
putchar('\n');
if((a1*b2-a2*b1)%31==0){ //連立方程式が解けない場合

printf("解くことができません\n");

}else{ //解ける場合

x = (c1*b2-b1*c2)%31/(a1*b2-a2*b1)%31;
y = ((a1*c2-a2*c1)%31/(a1*b2-a2*b1)%31);

printf("解x = %u, y = %u\n", x, y);
}
return 0;

} //ここまでがプログラム

