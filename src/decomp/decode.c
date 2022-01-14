#include "common.c"

// GCD for decode
OP ogcd(OP xx, OP yy)
{
  OP tt;

  while (odeg(yy) > T - 1)
  {
    tt = omod(xx, yy);
    xx = yy;
    yy = tt;
  }

  printpol(o2v(yy));
  printf(" =========yy\n");
  printpol(o2v(tt));
  printf(" =========tt\n");

  return tt;
}


OP ww[T] = {0};

OP bib(int i, OP d)
{
  int id, j;

  OP t[T] = {0};
  //omp_set_num_threads(omp_get_max_threads());
  id = omp_get_thread_num();
  t[id] = d;

  //#pragma omp parallel for
  for (j = 0; j < T; j++)
  {
    // #pragma omp critical
    if (i != j)
    {
      t[id] = omul(t[id], ww[j]);
    }
  }

  return t[id];
}

//多項式の形式的微分
OP bibun(vec a)
{
  OP w[T * 2] = {0};
  OP l = {0}, t[T] = {0}, d = {0};
  int i, j, n, id;
  vec tmp = {0};

  n = deg(a);
  printf("n=%d\n", n);
  if (n == 0)
  {
    printf("baka8\n");
    //  exit(1);
  }
  memset(ww, 0, sizeof(ww));
  // #pragma omp parallel num_threads(8)
  for (i = 0; i < T; i++)
  {
    ww[i].t[0].a = a.x[i];
    ww[i].t[0].n = 0;
    ww[i].t[1].a = 1;
    ww[i].t[1].n = 1;

    ////printpol(o2v(w[i]));
  }
  //  exit(1);

  tmp.x[0] = 1;
  //
  d = v2o(tmp);

// omp_set_num_threads(omp_get_max_threads());
//#pragma omp parallel num_threads(omp_get_max_threads()) //num_threads(TH)
  {
    //#pragma omp parallel for
//#pragma omp for schedule(static)
    for (i = 0; i < T; i++)
    {
      t[i] = bib(i, d);
    }
  }

  for (i = 0; i < T; i++)
    l = oadd(l, t[i]);

  return l;
}



//ユークリッドアルゴリズムによる復号関数
OP sendrier(OP f, OP s)
{
  int i, j, k, count = 0;
  OP r = {0}, w = {0}, e = {0}, l = {0};
  oterm t1, t2, d1, a, b;
  vec x = {0};
  unsigned short d = 0;
  OP h = {0};
  EX hh = {0};

  printf("in decode\n");
  printpol(o2v(s));
  printf("\nsyn===========\n");
  r = vx(f, s);
  //h=ogcd(f,s);

  if (odeg((r)) == 0)
  {
    printf("baka12\n");
    exit(1);
  }
  k = 0;
  // exit(1);
  x = chen(r);
  // exit(1);

  for (i = 0; i < T*2; i++)
  {
    printf("x[%d]=1\n", x.x[i]);
    if (x.x[i] == 0)
      k++;
    if (k > 1)
    {
      printf("baka0\n");
      printvec(o2v(f));
      //for (i = 0; i < N; i++)
      //printf("%d,", zz[i]);
      exit(1);
      //return f;
    }
  }
  //exit(1);

  //  printf("\n");

  printf("あっ、でる！\n");
  //  exit(1);

  if (odeg((r)) < T*2)
  {
    printpol(o2v(r));
    printf("baka5 deg(r)<T\n");
    exit(1);
  }

  w = bibun(x);
  //exit(1);
  //  w=oterml(w,d1);
  printpol(o2v(w));
  printf("@@@@@@@@@\n");
  //exit(1);

  //hh = xgcd(f, s);
  h = zgcd(f, s,T-1);
  //printpol(o2v(hh.d));
  printpol(o2v(h));
  //wait();

  //  exit(1);
  t1 = LT(r);

  t2.a = t1.a;
  t2.n = 0;

  if (odeg((w)) == 0)
  {
    printpol(o2v(w));
  }
  l = oterml(w, t2);

  j = deg(x) + 1;
  printf("%d\n", j);

  //    exit(1);

  for (i = 0; i < j; i++)
  {
    //if (x.x[i] > 0)
    {
      //e.t[i].a =
      //  gf[mlt(fg[trace(hh.d, x.x[i])], oinv(trace(l, x.x[i])))];
      e.t[i].a = gf[mlt(fg[trace(h, x.x[i])], oinv(trace(l, x.x[i])))];
      e.t[i].n = x.x[i];
    }
  }
  printpol(o2v(f));
  printf(" f============\n");
  printpol(o2v(l));
  printf(" l============\n");
  //  exit(1);

  for (i = 0; i < T; i++)
    if (gf[trace(h, x.x[i])] == 0)
      printf("h=0");
  //printf("\n");
  for (i = 0; i < T; i++)
    if (gf[oinv(trace(l, x.x[i]))] == 0)
      printf("l=0\n");
  //  printf("\n");

  return e;
}



//ユークリッドアルゴリズムによる復号関数
OP decode(OP f, OP s)
{
  int i, j, k, count = 0;
  OP r = {0}, w = {0}, e = {0}, l = {0};
  oterm t1, t2, d1, a, b;
  vec x = {0};
  unsigned short d = 0;
  OP h = {0};
  EX hh = {0};

  printf("in decode\n");
  printpol(o2v(s));
  printf("\nsyn===========\n");
  r = vx(f, s);

  if (odeg((r)) == 0)
  {
    printf("baka12\n");
    exit(1);
  }
  k = 0;
  // exit(1);
  x = chen(r);
  // exit(1);

  for (i = 0; i < T; i++)
  {
    printf("x[%d]=1\n", x.x[i]);
    if (x.x[i] == 0)
      k++;
    if (k > 1)
    {
      printf("baka0\n");
      printvec(o2v(f));
      //for (i = 0; i < N; i++)
      //printf("%d,", zz[i]);
      exit(1);
      //return f;
    }
  }
  //exit(1);

  //  printf("\n");

  printf("あっ、でる！\n");
  //  exit(1);

  if (odeg((r)) < T)
  {
    printpol(o2v(r));
    printf("baka5 deg(r)<T\n");
    exit(1);
  }

  w = bibun(x);
  //exit(1);
  //  w=oterml(w,d1);
  printpol(o2v(w));
  printf("@@@@@@@@@\n");
  //exit(1);

  //hh = xgcd(f, s);
  h = ogcd(f, s);
  //printpol(o2v(hh.d));
  printpol(o2v(h));
  //wait();

  //  exit(1);
  t1 = LT(r);

  t2.a = t1.a;
  t2.n = 0;

  if (odeg((w)) == 0)
  {
    printpol(o2v(w));
  }
  l = oterml(w, t2);

  j = deg(x) + 1;
  printf("%d\n", j);

  //    exit(1);

  for (i = 0; i < j; i++)
  {
    //if (x.x[i] > 0)
    {
      //e.t[i].a =
      //  gf[mlt(fg[trace(hh.d, x.x[i])], oinv(trace(l, x.x[i])))];
      e.t[i].a = gf[mlt(fg[trace(h, x.x[i])], oinv(trace(l, x.x[i])))];
      e.t[i].n = x.x[i];
    }
  }
  printpol(o2v(f));
  printf(" f============\n");
  printpol(o2v(l));
  printf(" l============\n");
  //  exit(1);

  for (i = 0; i < T; i++)
    if (gf[trace(h, x.x[i])] == 0)
      printf("h=0");
  //printf("\n");
  for (i = 0; i < T; i++)
    if (gf[oinv(trace(l, x.x[i]))] == 0)
      printf("l=0\n");
  //  printf("\n");

  return e;
}

