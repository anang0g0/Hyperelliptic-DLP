from fracModulo import extEuclidPoly, modPoly, multPoly, reModulo
import numpy
import sys



N=11
p=31

f=[-1,1,1,0,-1,0,1,0,0,1,-1]
m=[1,0,1,0,1,1,1,1]

if (len(sys.argv)>1):
	N=int(sys.argv[1])
if (len(sys.argv)>2):
	p=int(sys.argv[2])
if (len(sys.argv)>3):
	f=eval("["+sys.argv[3]+"]")
if (len(sys.argv)>4):
	m=eval("["+sys.argv[4]+"]")



print("Values used:")
print(" N=",N)
print(" p=",p)
print("========")
print("\nBob picks a polynomials (f):")

print("f(x)= ",f)
print ("\n",numpy.poly1d(f[::-1]))

D=[0]*(N+1)
D[0]=-1
D[N]=1


print("\n====Now we determine F_p ===")
[gcd_f,s_f,t_f]=extEuclidPoly(f,D)

f_p=modPoly(s_f,p)

print("F_p:",f_p)
print ("\n",numpy.poly1d(f_p[::-1]))


x=multPoly(f,m)
enc=reModulo(x,D,p)

x=multPoly(enc,f_p)
dec=reModulo(x,D,p)[:len(m)]

print("\n====Now we determine F_p ===")
print("Alice's Message:\t",m)
print ("\n",numpy.poly1d(m[::-1]))

print("\n====Encrypted message ===")
print("Encrypted message:\t",enc)
print ("\n",numpy.poly1d(enc[::-1]))

print("\n====Decrypted message ===")
print("Decrypted message:\t",dec)
print ("\n",numpy.poly1d(dec[::-1]))
