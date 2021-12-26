def inv(a, n)

  d = n
  x = 0
  s = 1
  while (a != 0)
    q = d / a
    r = d % a
    d = a
    a = r
    t = x - q * s
    x = s
    s = t
  end
  gcd = d  # $\gcd(a, n)$ 

  return ((x + n) % (n / d))
end

n=inv(2,31)

print n,"\n"
