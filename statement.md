# Bonjour!

Ce programme présente la factorisation par la méthode de Lenstra qui utilise des courbes elliptiques sur les entiers modulaires.
Cette version utilise des courbes de Montgomery, pour lesquelles les opérations arithmétiques peuvent être optimisées plus efficacement.
Les calculs n'utilisent que les coordonées x et z du point projectif.

```python runnable
from time import time
from random import randint
from math import gcd

def millerTest(a,d,n,r):
    # test de Miller pour un témoin a
    # Retourne faux si n est composé et vrai si n est probablement premier
    # d et r doivent vérifier n = 2^r * d + 1 avec d impair   
           
    x = pow(a, d, n) # Calcule a^d % n   
    if (x == 1  or x == n-1): 
       return True

    for _ in range(r):    
        x = (x * x) % n 
        if (x == 1):
            
            return False 
        if (x == n-1):
            return True    
    
    return False 

def isPrime(n, k=25): 
    # Test de primalité de Miller Rabin  
    # Si faux alors n est composé et si vrai alors n est probablement premier 
    # k determine le niveau de certitude : P(erreur) < 1/4**k
    
    if (n <= 1 or n == 4):
        return False 
    if (n <= 5):
        return True   
    
    # Trouver d et r tels que n = 2^r * d + 1 avec d impair 
    d = n - 1 
    r = 0
    while (d&1 == 0): 
        d  >>= 1 
        r += 1 
    
    # Effectuer k tests de Miller
    for i in range(k):
        a = randint(2,n-2) 
        if (not millerTest(a, d, n, r)):
              return False  
    return True

def nextPrime(n):
    # premier suivant n
    while not isPrime(n):
        n += 1
    return n

# Crible d'Eratosthenes
def sieve(n):
    b = [True] * n
    ps = []
    for p in range(2, n):
        if b[p]:
            ps.append(p)
            for i in range(p, n, p):
                b[i] = False
    return ps

# Addition de deux points d'une courbe de Montgomery
# Calcule P+Q sur base de P, Q et P-Q
def addPoints(xp, zp, xq, zq, xpq, zpq, n):
    
    u = (xp + zp) * (xq - zq) %n
    v = (xp - zp) * (xq + zq) %n
    w = (u + v) 
    t = (u - v) 
    w = (w * w) %n
    t = (t * t) %n
    X = (w * zpq) %n
    Z = (t * xpq) %n
    return X, Z

# Duplication d'un point d'une courbe de Montgomery
def duplicatePoint(x, z, a, n):
    
    x2 = x*x %n
    z2 = z*z %n
    xz = x*z %n
    u = (x2-z2)
    z = 4*xz*(x2+a*xz+z2)%n
    x = u*u %n
      
    return x, z

# Multiplication d'un point par un entier 
# Algorithme de l'échelle de Montgomery
def montgomeryLadder(k, px, pz, a, n):
    qx, qz = px, pz
    rx, rz = duplicatePoint(px, pz, a, n)
    for c in bin(k)[3:]:
        if c == '1':
            qx, qz = addPoints(rx, rz, qx, qz, px, pz, n)
            rx, rz = duplicatePoint(rx, rz, a, n)
        else:
            rx, rz = addPoints(qx, qz, rx, rz, px, pz, n)
            qx, qz = duplicatePoint(qx, qz, a, n)

    return qx, qz

# Algorithme de Lenstra (version Montgomery)
# Limit specifie le maximum des premiers à tester
# Primes est une liste de nombres premiers < limit
def montgomery(n, limit=1000, primes=None):
  
    if primes is None:
        primes = sieve(limit)
   
    g = n
    while g == n:
        # x et z sont aléatoires
        x,z = q = randint(0, n - 1), randint(0, n - 1)
        # a est aléatoire, b est calculé pour que le point soit sur la courbe
        a = randint(0, n - 1)
        b = (x*x*x + a * x*x*z + x*z*z) % n
        g = gcd(( a * a - 4)* b, n)  # vérification de singularité
    # Par chance, g peut être un facteur
    if g > 1:
        return g
    # Multiplier le point q par le PPCM des nombres (1, ..., limit)
    i=0
    for p in primes:
        pp = p
        while pp < limit:
            x,z = q            
            q = x,z = montgomeryLadder(p, x, z, a, n)                                    
            i += 1    
            if i%1000==0 :                
                g = gcd(z,n)
                if g>1:
                    #print('g=',g, 'n=',n)
                    return g                   
            pp = p * pp
    
    g = gcd(z, n)
    #print('montgomery : g=',g)    
    return g

def factorECM(n, k=25, limit=10000, div=montgomery, primes=None):
    # Décompose n en facteurs probablement premiers
    if primes is None:
         primes = sieve(limit)    
    if n in primes:
       return [n]
    for p in primes:
        if n%p==0:
            return [p]+factorECM(n//p, k=k, limit=limit, div=div, primes=primes)
    if isPrime(n,k): # on pense que n est premier
        return [n]
    g = 1
    i = 0  
    while g==1 or g==n: # tant qu'on n'a pas de facteur
        # on essaye avec d'autres valeurs
        i +=1
        if i%10==0:
            print("retry",i,div.__name__)
        g=div(n, limit=limit, primes=primes)  
                       
    return factorECM(g, k=k, limit=limit, div=div, primes=primes)+factorECM(n//g, k=k, limit=limit, div=div, primes=primes)
 
for i in range(5):
    a=nextPrime(randint(1e8,1e12))    
    b=nextPrime(randint(1e8,1e10))
    c=randint(1e8,1e10) 
    n=a*b*c
    print('n=',n)
    t=time()
    f=factorECM(n)  
    t=time()-t  
    print('les facteurs de n sont',f,'t=',t,'s')
    print()


```
