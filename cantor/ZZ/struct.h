#include <NTL/ZZ.h>

#define DEG 512

/* -*- mode: C; coding:utf-8 -*- */
#define N 100
#define E 4

NTL_CLIENT

//monomial
typedef struct{
  int n; //単項式の次数
  ZZ a; //単項式の係数
} oterm;

//polynomial
typedef struct{
  oterm t[DEG]; //単項式の配列として多項式を表現する
} OP;

typedef struct {
  ZZ x[DEG]; //配列の添字を次数に、配列の値を係数に持つ多項式の表現
} vec;


//extra gcd
typedef struct{
  OP u; //inverse of polynomial?
  OP v; //error locater
  OP d; //gcd
  OP h;
} EX;




typedef struct  {
  ZZ x;
  ZZ y;
} PO;

typedef struct {
  OP u;
  OP v;
} Div;


  struct BF {
    unsigned int a : 10;
    unsigned int b : 8;
    unsigned int c : 12;
    unsigned int d : 2;
  };

  union {
    struct BF f;
    unsigned int u;
  } dd;
