
OP diviser(PO o, PO m)
{
  ZZ t1[2][3]={to_ZZ("0")}, cc[2]={to_ZZ("0")};
  vec c1, c2;
  int i, j, k;
  ZZ a, b;
  OP f;

  t1[0][0] = o.x; // o.t[1].a;
  t1[1][0] = m.x; // t[0].a;
  t1[0][1] = 1;   // t[1].a;
  t1[1][1] = 1;   // t[0].a;
  t1[0][2] = o.y;
  t1[1][2] = m.y;

  if (t1[0][0] == 0)
  {
    for (i = 0; i < 2; i++)
      t1[0][i] = (t1[0][i] + t1[1][i]) % P;
  }
  if (t1[1][1] == 0)
  {
    for (i = 0; i < 2; i++)
      t1[1][i] = (t1[1][i] + t1[0][i]) % P;
  }

  cc[0] = inv(t1[0][0], P)%P;
  cout << cc[0] << endl;
  if(cc[0]>P)
  exit(1);
  ZZ z;
  z = t1[1][0]%P;
  cout << "z=" << z << endl;

  cc[0] = (cc[0] * z) % P;
  cout << "c0=" << cc[0] << endl;
  // exit(1);

  for (i = 0; i < 2; i++)
  {
    for (j = 0; j < 3; j++)
      cout << t1[i][j] << endl;
    printf("\n");
  }
  printf("\n");

  for (j = 0; j < 3; j++)
  {
    t1[1][j] = (P+((t1[0][j] * cc[0]) % P - t1[1][j])) % P;
    cout << "," << t1[1][j] << endl;
  }
  printf("\n");
  b = (t1[1][2] * inv(t1[1][1], P) % P);
  cout << "b="<<  b << endl;
  for (i = 0; i < 3; i++)
    a = (t1[0][2] - (t1[0][1] * b) % P);
  if (a < 0)
    a += P;
  a = (inv(t1[0][0], P) * a) % P;
  cout << "a=" << a << endl;

  c1.x[0] = b;
  c1.x[1] = a;
  cout << c1.x[0] << " , " << c1.x[1] << endl;
  // exit(1);
  f = v2o(c1);
  printpoln(o2v(f));
  // exit(1);

  return f;
}

// ランダムな因子の生成
Div gendiv(OP f)
{
  int count = 0;
  PO a, b, e;
  OP d1 = {0}, d2 = {0}, c = {0}, d = {0}, vv1 = {0}, vv2 = {0}, uu1, v;
  //  ZZ  x, y, i, j, k,

  Div D = {0};
  vec v1, v2, z1, z2, ff, vx;
  EX V;

  vx = o2v(f);
  /// srand(clock());
  v1.x[1] = to_ZZ("1");
  v2.x[1] = to_ZZ("1");
//  while(1)
  {
    do{
    a = tr1e(vx.x[4], vx.x[3], vx.x[2], vx.x[1], vx.x[0], P); // cofficient of function
    b = tr1e(vx.x[4], vx.x[3], vx.x[2], vx.x[1], vx.x[0], P); // cofficient of function
    }while((a.x == b.x) || (a.x ==0 || b.x ==0));
    cout << "P= " << P << endl;
    cout << "P-x= " << P - a.x << endl;
    cout << "x= " << a.x << endl;
    v1.x[0] = P - a.x;
    
    printpol(v1);
    printf("ppppppppppppp\n");
    c = v2o(v1);

    v2.x[0] = P - b.x;
    d = v2o(v2);
    d2 = diviser(a, b);
    printpoln(o2v(d2));
    // exit(1);

    // d2 = v2o(z1);
    d1 = omul(c, d);
    printpoln(o2v(d1));
    // exit(1);
    D.u = d1;
    D.v = d2;
    if (chkdiv(D, f) != -1)
    {
      printf("line\n");
      return D;
    }else{
      printf("just say buggy\n");
    exit(1);
    }
    // 
  } //while (chkdiv(D, f) == -1);

  if (chkdiv(D, f) == -1)
  {
    printf("so buggy!\n");
    exit(1);
  }
  return D;
}
