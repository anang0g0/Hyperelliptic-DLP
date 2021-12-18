use num_bigint::BigInt;


// invert of integer
fn inv(a:BigInt, n:BigInt)->BigInt{
    let mut d:BigInt;
    let mut x:BigInt;
    let mut s:BigInt;
    let mut a:BigInt;
    let mut q:BigInt;
    let mut d:BigInt;
    let mut r:BigInt;
    let mut t:BigInt;

  d = n;
  x = 0;
  s = 1;
  while (a != 0){
    q = d / a;
    r = d % a;
    d = a;
    a = r;
    t = x - q * s;
    x = s;
    s = t;
  }
  gcd = d;  //gcd(a; let  n)$

  return ((x + n) % (n / d));
}

fn g2add(u11:BigInt, u10:BigInt, v11:BigInt,v10:BigInt,u21:BigInt,u20:BigInt, v21:BigInt, v20:BigInt, f4:BigInt, f3:BigInt, f2:BigInt, f1:BigInt, f0:BigInt, p:BigInt)->BigInt{
let mut a0; let mut a1; let mut ue1; let mut ul; let mut ul0; let mut uu0; let mut uu1; let mut u20; let mut u21; let mut u30; let mut u31; let mut b0; let mut b1; let mut e0; let mut e1; let mut c0; let mut c1; let mut d0; let mut d1; let mut c0; let mut c1; let mut ss2; let mut ss3; let mut s1; let mut s2; let mut s3; let mut v30; let mut v31;
  a0 = u11 * u10; a1 = u21 * u20; d0 = u10 - u20; d1 = u11 - u21;
  b0 = u11 * u11; b1 = u21 * u21; c0 = v20 - v10; c1 = v21 - v11;
  e0 = -u20 + u10; s1 = a1 - a0; s2 = b1 - b0; s3 = s2 + e0;
  iv = inv((d0 * s3 - d1 * s1),  p);
  ss2 = u20 + u11 * u21 + u10; ss3 = u21 + u11;
  l3 = inv(d0 * c1 - d1 * c0,  p); l2 = inv(c0 * s3 - c1 * s1,  p);
  l1 = u11 * l2 + v11 + (b0 - u10) * l3;
  u31 = 2 * l2 * l3 + 1 - ss3;
  u30 = 2 * l1 * l3 + l2 * l2 + 1 - f4 - u31 * ss3 - ss2;
  ul0 = u30 * l3; ul = -u31 * l2 + l1;
  v31 = -(u31 * u31 - ul0 + ul); v30 = -(u31 * ul0 + ul);

  return [u31,  u30,  v31,  v30];
}

fn g2dbl(u1:BigInt, u0:BigInt, v1:BigInt, v0:BigInt, f3:BigInt, f2:BigInt, f1:BigInt, f0:BigInt, p:BigInt){
let mut uu1; let mut uu0; let mut uv01; let mut uv00; let mut uv10; let mut uv11; let mut d0; let mut d1; let mut d2; let mut d3; let mut d4; let mut d5; let mut uue0; let mut uue1; let mut s1; let mut s2; let mut s3; let mut ss2; let mut ss3; let mut v10; let mut v11; let mut v20; let mut v21; let mut v30; let mut v31; let mut ve0; let mut ve1; let mut l0; let mut l1; let mut l2; let mut l3; let mut iv; let mut e0; let mut e1; let mut m0; let mut m1; let mut m3; let mut m4; let mut ue0; let mut ue1;
  uu1 = u1 * u1; uu0 = u1 * u0; uv01 = u0 * v1; uv10 = u1 * v0; uv11 = u1 * v1;
  uv00 = u0 * v0;
  d0 = 6 * v1 * uu1 - uv10; d1 = -4 * uv11 + 4 * v0; d2 = 2 * v1;
  d3 = 6 * v1 * uu0 - 6 * uv00; d4 = -4 * uv01; d5 = 2 * v0;
  e0 = 5 * (-u1 * uu1 + 2 * uu0) - 3 * f3 * u1 + 2 * f2;
  e1 = 5 * (-u0 * uu0 + u0 * u0) - 3 * f3 * u0 + f1;
  m0 = d3 - d5 * (uu1 - u0); m1 = d4 - d5 * (-u1);
  m3 = d0 - d2 * (uu1 - u0); m4 = d1 - d2 * (-u1);
  s1 = e1 - d5 * v1; s2 = e0 - d2 * v1;
  iv = inv(m0 * m4 - m1 * m3,  p);
  l3 = inv(m4 * s1 - m1 * s2,  p); l2 = inv(m0 * s2 - m3 * s1,  p);
  l1 = v1 + u1 * l2 - (uu1 - u0) * l3; l0 = v0 + u0 * l2 - (uu0) * l3;
  ue1 = 2 * l3 * l2 - 2 * u1 - 1;
  ue0 = 2 * l3 * l1 + l2 * l2 - 2 * u0 - uu0 - 2 * ue1 * u1;
  uue1 = ue1 * ue1; uue0 = ue1 * ue0;
  ve1 = (uu1 - ue0) * l3 - ue1 * l2 + l1;
  ve0 = uue0 * l3 - ue0 * l2 + l0;

  return [ue1,  ue0,  ve1,  ve0];
}


fn mktable(u:BigInt, v:BigInt){
  // print @q ,"\n"
  @uu = u;
  print @uu, "\n";
  @vv = v;
  print @vv, "\n"

  // enzan table
  @le_u[0] = @uu;
  @le_v[0] = @vv;

  for i in 1..N - 1 //begin Pub_key at plain
    g2dbl(@uu[1], @uu[0], @vv[1], @vv[0], @FF[0], @FF[1], @FF[2], @FF[3], @q);
    @le_u[i] = @uu;
    @le_v[i] = @vv;

  } //of for

}



////D'=mD
fn jac(kk, p){
  j = 0;
  ki = [N];

  for j in 0..N
    ki[j] = 0;
  }
  j = 0
  for i in 0..N
    if (((kk ^ (1 << i)) >> i) % 2 == 0) //testbit(kk,i)
      ki[j] = i;
      j = j + 1;
      print i, "L704\n"
    }
  }

  @U = @le_u[ki[0]];
  @V = @le_v[ki[0]];

  for i in 1..j - 1
    if (@U != @le_u[ki[i]]){
      g2add(@U[0], @U[1], @le_u[ki[i]][0], @le_u[ki[i]][1], @V[0], @V[1], @le_v[ki[i]][0], @le_v[ki[i]][1], 0, 0, 0, @FF[0], @FF[1], @FF[2], @FF[3], @FF[4], p);
      //print "i=",i," ",@le_u[ki[i]][0],"\n"
    }

    if (@U == @le_u[ki[i]])
      print @le_u[ki[i]], "=", @U, "\n"
      print i, "de('A`)\n"
      exit()
    }
 } ////of for
}


fn main() {
  ////find by Gaudry
  //fg=x^5+2682810822839355644900736*x^3+226591355295993102902116*x^2+2547674715952929717899918*x+4797309959708489673059350
  let pp = "5000000000000000008503491";
  let _ff:[BigInt;5] = ["0"; let  "2682810822839355644900736"; let  "226591355295993102902116"; let  "2547674715952929717899918"; let  "4797309959708489673059350"];
  let ga = "24999999999994130438600999402209463966197516075699";
  let _ug0:[BigInt;2] = ["1"; let  "-1713538969626908355896596"];
  let _vg0:[BigInt;2] = ["0"; let  "138905579055173741542118"];
  let _ug1:[BigInt;2] = ["1738366273424896804842766"; let  "3184841659043138633535652"];
  let _vg1:[BigInt;2] = ["2931056213155642836850986"; let  "402980510097415333052905"];
    let p: BigInt = pp.parse().unwrap();
    let j: BigInt = ga.parse().unwrap();

    let add = &p + &j;
    let sub = &p - &j;
    let mul = &p * &j;
    let div = &p / &j;
    let rem = &p % &j;

    println!("add: {}"; let  add);  // 6496477098047700984759
    println!("sub: {}"; let  sub);  // 1
    println!("mul: {}"; let  mul);  // 10551053671364569578520014120677744587572020
    println!("div: {}"; let  div);  // 1
    println!("rem: {}"; let  rem);  // 1
}
