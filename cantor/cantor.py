#
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
    xs = a, 1, 0
    ys = b, 0, 1
    while ys[0] != 0:
        q, z = divmod(xs[0], ys[0])
        xs, ys = ys, (z, xs[1] - q * ys[1], xs[2] - q * ys[2])
    return xs

# 逆数表の生成
def make_inv(n):
    inv = [0] * n
    for x in range(1, n):
        g, y, _ = ext_euclid(x, n)
        if g != 1: raise Exception("not prime number!")
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
        if not len(xs): xs = np.array([0])
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
        zn = xn + yn - 1
        zs = np.full(zn, 0, dtype=np.int32)
        for i in range(yn):
            zs[i:i+xn] += ys[i] * xs
        return self.normalize(zs)

    # 割り算
    def polydiv(self, xs, ys):
        xn = len(xs)
        yn = len(ys)
        zs = xs.copy()
        qs = []
        for _ in range(xn - yn + 1):
            temp = (zs[0] * self.inv[ys[0]]) % self.num
            zs[:yn] += self.num - temp * ys
            qs.append(temp)
            zs = zs[1:] % self.num
        if qs == []: qs = [0]
        return np.array(qs), self.normalize(zs)

    # 因数分解
    def factorization(self, xs):
        ps = []
        g = functools.reduce(lambda x, y: math.gcd(x, y), xs)
        if g > 1:
            xs = xs // g
            ps.append((np.array([g]), 1))
        m = (len(xs) - 1) // 2 + 1
        zs = itertools.product(range(self.num), repeat=m)
        for _ in range(self.num): next(zs)
        for z in zs:
            ys = np.trim_zeros(np.array(z), 'f')
            if len(xs) <= len(ys): break            
            c = 0
            while len(xs) > len(ys):
                p, q = self.polydiv(xs, ys)
                if q[0] == 0:
                    xs = p
                    c += 1
                else:
                    break
            if c > 0:
                if np.all(xs == ys):
                    c += 1
                    xs = np.array([1])
                ps.append((ys, c))
        if len(xs) > 1: ps.append((xs, 1))
        return ps
    
# 多項式の値を求める (ホーナー法)
def polyval(xs, n):
    return functools.reduce(lambda x, y: x * n + y, xs)


def xgcd(a, b):
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b != 0:
        q, a, b = a // b, b, a % b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    return a, x0, y0

def modinv(a, m):
    g, x, y = xgcd(a, m)
    if g != 1:
        raise Exception('modular inverse does not exist')
    else:
        return x % m
    


def legendre(a, p):
    return pow(a, (p - 1) // 2, p)

def tonelli_shanks(a, p):
    if legendre(a, p) != 1:
        raise Exception("not a square (mod p)")
    # Step 1. By factoring out powers of 2, find q and s such that p - 1 = q 2^s with Q odd
    q = p - 1
    s = 0
    while q % 2 == 0:
        q >>= 1
        s += 1
    # Step 2. Search for a z in Z/pZ which is a quadratic non-residue
    for z in range(2, p):
        if legendre(z, p) == p - 1:
            break
    # Step 3.
    m = s
    c = pow(z, q, p) # quadratic non residue
    t = pow(a, q, p) # quadratic residue
    r = pow(a, (q + 1) // 2, p)
    # Step 4.
    t2 = 0
    while True:
        if t == 0: return 0
        if t == 1: return r
        t2 = (t * t) % p
        for i in range(1, m):
            if t2 % p == 1:
                break
            t2 = (t2 * t2) % p
        b = pow(c, 1 << (m - i - 1), p)
        m = i
        c = (b * b) % p
        t = (t * c) % p
        r = (r * b) % p


if __name__ == '__main__':
    from sympy.ntheory.residue_ntheory import sqrt_mod

    ttest = [
        (10, 13), (56, 101), (1030, 10009), (44402, 100049),
        (665820697, 1000000009), (881398088036, 1000000000039),
        (41660815127637347468140745042827704103445750172002, 10**50 + 577)
    ]

    for n, p in ttest:
        r = tonelli_shanks(n, p)
        roots = [r, p-r]
        r2 = sqrt_mod(n, p)
        roots2 = [r2, p-r2]
        assert (roots[0] * roots[0] - n) % p == 0
        assert (roots[1] * roots[1] - n) % p == 0
        assert min(roots, roots2)
        assert max(roots, roots2)
        print("n = %d p = %d" % (n, p), end=' ')
        print("roots : %d %d" % (r, p - r))

# ここが動かない！
def inv(a, n):
  a=np.array(a)
  n=np.array(n)
  f=np.array([0,0,0,0])
  gf=GF(11)

  d = n;
  x = np.array([0,0,0,0]);
  s = np.array([0,0,0,1]);
  print(f)
  r=np.array([0,0,0,0])
  q=np.array([0,0,0,0])
  while (a.all() != 0):
    q, r = gf.polydiv(d , a);
    d = a;
    a = r;
    tt = gf.polymul(q , s)
    t = gf.polysub(x, tt)
    x = s;
    s = t;
  
  gcd = d;
  print(r)
  exit()
  f,g = gf.polydiv(n,d)
  o=gf.polyadd(x,n)
  h,e =gf.polydiv(o,f)
  return q, r



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

        

test(2, [1,0,1,1], [1,0,1])
test(3, [2,1,0,2], [2,0,1])
test(5, [4,3,2,1], [3,2,1])
gf=GF(11)
f=[1,0,3,7,1,2]
v1=[1,9]
u1=[1,7,10]
u2=[1,0,10]
v2=[7,9]
t1=np.array(v2)
t2=np.array(v1)
ff=np.array(f)
uu2=np.array(u2)
uu1=np.array(u1)
a=gf.polymul(t1,t1)
c=gf.polysub(ff,a)
p,q=gf.polydiv(c,uu2)
print(p)
print(q)
#exit()

k=p
a=gf.polysub(t2,t1)
print(inv(uu2,uu1))
print(a)
#print(b)
exit()

a,b=gf.polydiv(p1,p2)
#inv([2,1,0,2],[3,2,1])
print(p1)

exit()

#a=inv(12,5)
#print(a)
m = 24999999999994130438600999402209463966197516075699

a = 12
#m = 5

try:
    res=modinv(a,m)
    print("The required modular inverse is: "+ str(res))

except:
    print('The modular inverse does not exist.')

