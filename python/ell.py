p = random_prime(100000000)
a = random_prime(10)
F = GF(p)
for b in range(1, p):
    E = EllipticCurve(F, [a, b])
    if is_prime(E.order()):
        print(E)
        break
    E = EllipticCurve(F, [a, -b])
    if is_prime(E.order()):
        print(E)
        break

print('p =', p)
print('a =', a)
print('b =', b)
print('#E =', E.order())
