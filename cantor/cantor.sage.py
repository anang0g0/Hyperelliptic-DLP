

# This file was *autogenerated* from the file cantor.sage
from sage.all_cmdline import *   # import sage library

_sage_const_1 = Integer(1); _sage_const_0 = Integer(0); _sage_const_2 = Integer(2); _sage_const_10 = Integer(10); _sage_const_13 = Integer(13); _sage_const_56 = Integer(56); _sage_const_101 = Integer(101); _sage_const_1030 = Integer(1030); _sage_const_10009 = Integer(10009); _sage_const_44402 = Integer(44402); _sage_const_100049 = Integer(100049); _sage_const_665820697 = Integer(665820697); _sage_const_1000000009 = Integer(1000000009); _sage_const_881398088036 = Integer(881398088036); _sage_const_1000000000039 = Integer(1000000000039); _sage_const_41660815127637347468140745042827704103445750172002 = Integer(41660815127637347468140745042827704103445750172002); _sage_const_50 = Integer(50); _sage_const_577 = Integer(577); _sage_const_11 = Integer(11); _sage_const_3 = Integer(3); _sage_const_5 = Integer(5); _sage_const_4 = Integer(4); _sage_const_7 = Integer(7); _sage_const_9 = Integer(9); _sage_const_24999999999994130438600999402209463966197516075699 = Integer(24999999999994130438600999402209463966197516075699); _sage_const_12 = Integer(12)#
# polygf.py : GF(p) で多項式を計算する (p は素数)
#
#             Copyright (C) 2018 Makoto Hiroi
#
import numpy as np
import functools
import itertools
import math
import galois
import sage

# 拡張ユークリッドの互除法
def ext_euclid(a, b):
    xs = a, _sage_const_1 , _sage_const_0 
    ys = b, _sage_const_0 , _sage_const_1 
    while ys[_sage_const_0 ] != _sage_const_0 :
        q, z = divmod(xs[_sage_const_0 ], ys[_sage_const_0 ])
        xs, ys = ys, (z, xs[_sage_const_1 ] - q * ys[_sage_const_1 ], xs[_sage_const_2 ] - q * ys[_sage_const_2 ])
    return xs

# 逆数表の生成
def make_inv(n):
    inv = [_sage_const_0 ] * n
    for x in range(_sage_const_1 , n):
        g, y, _ = ext_euclid(x, n)
        if g != _sage_const_1 : raise Exception("not prime number!")
        inv[x] = (y + n) % n
    return inv

class GF:
    def __init__(self, n):
        self.num = n
        self.inv = make_inv(n)

    # 正規化
    def normalize(self, xs):
        xs %= self.num
        xs = np.trim_zeros(xs, 'f')
        if not len(xs): xs = np.array([_sage_const_0 ])
        return xs

    # 足し算
    def polyadd(self, xs, ys):
        xn = len(xs)
        yn = len(ys)
        if xn >= yn:
            zs = xs.copy()
            zs[-yn:] += ys
        else:
            zs = ys.copy()
            zs[-xn:] += xs
        return self.normalize(zs)

    # 引き算
    def polysub(self, xs, ys):
        return self.polyadd(xs, self.num - ys)

    # 掛け算
    def polymul(self, xs, ys):
        xn = len(xs)
        yn = len(ys)
        zn = xn + yn - _sage_const_1 
        zs = np.full(zn, _sage_const_0 , dtype=np.int32)
        for i in range(yn):
            zs[i:i+xn] += ys[i] * xs
        return self.normalize(zs)

    # 割り算
    def polydiv(self, xs, ys):
        xn = len(xs)
        yn = len(ys)
        zs = xs.copy()
        qs = []
        for _ in range(xn - yn + _sage_const_1 ):
            temp = (zs[_sage_const_0 ] * self.inv[ys[_sage_const_0 ]]) % self.num
            zs[:yn] += self.num - temp * ys
            qs.append(temp)
            zs = zs[_sage_const_1 :] % self.num
        if qs == []: qs = [_sage_const_0 ]
        return np.array(qs), self.normalize(zs)

    # 因数分解
    def factorization(self, xs):
        ps = []
        g = functools.reduce(lambda x, y: math.gcd(x, y), xs)
        if g > _sage_const_1 :
            xs = xs // g
            ps.append((np.array([g]), _sage_const_1 ))
        m = (len(xs) - _sage_const_1 ) // _sage_const_2  + _sage_const_1 
        zs = itertools.product(range(self.num), repeat=m)
        for _ in range(self.num): next(zs)
        for z in zs:
            ys = np.trim_zeros(np.array(z), 'f')
            if len(xs) <= len(ys): break            
            c = _sage_const_0 
            while len(xs) > len(ys):
                p, q = self.polydiv(xs, ys)
                if q[_sage_const_0 ] == _sage_const_0 :
                    xs = p
                    c += _sage_const_1 
                else:
                    break
            if c > _sage_const_0 :
                if np.all(xs == ys):
                    c += _sage_const_1 
                    xs = np.array([_sage_const_1 ])
                ps.append((ys, c))
        if len(xs) > _sage_const_1 : ps.append((xs, _sage_const_1 ))
        return ps
    
# 多項式の値を求める (ホーナー法)
def polyval(xs, n):
    return functools.reduce(lambda x, y: x * n + y, xs)


def xgcd(a, b):
    x0, y0, x1, y1 = _sage_const_1 , _sage_const_0 , _sage_const_0 , _sage_const_1 
    while b != _sage_const_0 :
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0

def modinv(a, m):
    g, x, y = xgcd(a, m)
    if g != _sage_const_1 :
        raise Exception('modular inverse does not exist')
    else:
        return x % m
    


def legendre(a, p):
    return pow(a, (p - _sage_const_1 ) // _sage_const_2 , p)

def tonelli_shanks(a, p):
    if legendre(a, p) != _sage_const_1 :
        raise Exception("not a square (mod p)")
    # Step 1. By factoring out powers of 2, find q and s such that p - 1 = q 2^s with Q odd
    q = p - _sage_const_1 
    s = _sage_const_0 
    while q % _sage_const_2  == _sage_const_0 :
        q >>= _sage_const_1 
        s += _sage_const_1 
    # Step 2. Search for a z in Z/pZ which is a quadratic non-residue
    for z in range(_sage_const_2 , p):
        if legendre(z, p) == p - _sage_const_1 :
            break
    # Step 3.
    m = s
    c = pow(z, q, p) # quadratic non residue
    t = pow(a, q, p) # quadratic residue
    r = pow(a, (q + _sage_const_1 ) // _sage_const_2 , p)
    # Step 4.
    t2 = _sage_const_0 
    while True:
        if t == _sage_const_0 : return _sage_const_0 
        if t == _sage_const_1 : return r
        t2 = (t * t) % p
        for i in range(_sage_const_1 , m):
            if t2 % p == _sage_const_1 :
                break
            t2 = (t2 * t2) % p
        b = pow(c, _sage_const_1  << (m - i - _sage_const_1 ), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p


if __name__ == '__main__':
    from sympy.ntheory.residue_ntheory import sqrt_mod

    ttest = [
        (_sage_const_10 , _sage_const_13 ), (_sage_const_56 , _sage_const_101 ), (_sage_const_1030 , _sage_const_10009 ), (_sage_const_44402 , _sage_const_100049 ),
        (_sage_const_665820697 , _sage_const_1000000009 ), (_sage_const_881398088036 , _sage_const_1000000000039 ),
        (_sage_const_41660815127637347468140745042827704103445750172002 , _sage_const_10 **_sage_const_50  + _sage_const_577 )
    ]

    for n, p in ttest:
        r = tonelli_shanks(n, p)
        roots = [r, p-r]
        r2 = sqrt_mod(n, p)
        roots2 = [r2, p-r2]
        assert (roots[_sage_const_0 ] * roots[_sage_const_0 ] - n) % p == _sage_const_0 
        assert (roots[_sage_const_1 ] * roots[_sage_const_1 ] - n) % p == _sage_const_0 
        assert min(roots, roots2)
        assert max(roots, roots2)
        print("n = %d p = %d" % (n, p), end=' ')
        print("roots : %d %d" % (r, p - r))

def inv(a, n):
  a=np.array(a)
  n=np.array(n)
  f=np.array([_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ])
  gf=GF(_sage_const_11 )

  d = n;
  x = np.array([_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_0 ]);
  s = np.array([_sage_const_0 ,_sage_const_0 ,_sage_const_0 ,_sage_const_1 ]);
  print(f)

  while (a.all() != f.all()):
    q,r = gf.polydiv(d , a);
    d = a;
    a = r;
    tt = gf.polymul(q , s)
    t = gf.polysub(x, tt)
    x = s;
    s = t;
  
  gcd = d;
  f,g = gf.polydiv(n,d)
  o=gf.polyadd(x,n)
  h,e =gf.polydiv(o,f)
  return e;



def test(p, xs, ys):
    xs = np.array(xs)
    ys = np.array(ys)
    gf = GF(p)
    for x, y in [(xs, xs), (xs, ys), (ys, xs), (ys, ys)]:
        print(x, "+", y, "=", gf.polyadd(x, y))
        print(x, "-", y, "=", gf.polysub(x, y))
        print(x, "*", y, "=", gf.polymul(x, y))
        p, q = gf.polydiv(x, y)
        print(x, "/", y, "=", p, q)

        

test(_sage_const_2 , [_sage_const_1 ,_sage_const_0 ,_sage_const_1 ,_sage_const_1 ], [_sage_const_1 ,_sage_const_0 ,_sage_const_1 ])
test(_sage_const_3 , [_sage_const_2 ,_sage_const_1 ,_sage_const_0 ,_sage_const_2 ], [_sage_const_2 ,_sage_const_0 ,_sage_const_1 ])
test(_sage_const_5 , [_sage_const_4 ,_sage_const_3 ,_sage_const_2 ,_sage_const_1 ], [_sage_const_3 ,_sage_const_2 ,_sage_const_1 ])
gf=GF(_sage_const_11 )
f=x**_sage_const_5 +_sage_const_3 *x**_sage_const_3 +_sage_const_7 *x**_sage_const_2 +x+_sage_const_2 
v1=x+_sage_const_9 
v2=_sage_const_7 *x+_sage_const_9 
u1=x**_sage_const_2 +_sage_const_7 *x+_sage_const_10 
u2=x**_sage_const_2 +_sage_const_10 
a=(t1*t1)
c=(ff-a)
p=c/uu2
r=c%uu2
print(p)
print(r)
#exit()

k=p
a=gf.polysub(t2,t1)
b=uu2.inverse_mod(uu1)
print(a)
print(b)
exit()

a,b=gf.polydiv(p1,p2)
#inv([2,1,0,2],[3,2,1])
print(p1)

exit()

#a=inv(12,5)
#print(a)
m = _sage_const_24999999999994130438600999402209463966197516075699 

a = _sage_const_12 
#m = 5

try:
    res=modinv(a,m)
    print("The required modular inverse is: "+ str(res))

except:
    print('The modular inverse does not exist.')

