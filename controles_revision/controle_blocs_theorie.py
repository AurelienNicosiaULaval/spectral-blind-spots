"""Independent finite-dimensional verification of the stated Fourier ranks.
Numerical ranks illustrate, but do not prove, the infinite-dimensional theorem.
"""
import numpy as np
from math import pi

def row(k,origin):
    return np.array([np.cos(k*origin),np.sin(k*origin)])

def matrix(k,l,angles):
    return np.array([np.kron(row(k,a),row(l,b)) for a,b in angles])

def rank(k,l,angles):
    return np.linalg.matrix_rank(matrix(k,l,angles),tol=1e-10)

# Four sine/cosine views at each common shift, with an ample, non-resonant
# grid of shifts, for all ordered positive Fourier blocks up to order 12.
origins=[(a+j*pi/2,a+m*pi/2) for a in np.arange(53)*2*pi/53 for j in (0,1) for m in (0,1)]
for k in range(1,13):
    for l in range(1,13):
        expected=3 if k==l and k%2==0 else 4
        assert rank(k,l,origins)==expected,(k,l,rank(k,l,origins),expected)
        if expected==3:
            # Matrix of sin(k(x-theta)): cos(x)sin(theta) - sin(x)cos(theta).
            blind=np.array([0.,-1.,1.,0.])
            assert np.max(np.abs(matrix(k,l,origins)@blind))<1e-11

# Published finite configuration: count null dimensions for all ordered blocks.
finite=[(a+j*pi/2,a+m*pi/2) for a in (0,pi/8,pi/3,pi/4) for j in (0,1) for m in (0,1)]
nullity=sum(4-rank(k,l,finite) for k in range(1,9) for l in range(1,9))
assert nullity==11,nullity
for delta in (pi/7,pi/np.sqrt(2)):
    crossed=[(a,b) for a in (0,delta) for b in (0,delta)]
    for k in range(1,9):
        for l in range(1,9):
            expected=(1 if abs(np.sin(k*delta))<1e-10 else 2)*(1 if abs(np.sin(l*delta))<1e-10 else 2)
            assert rank(k,l,crossed)==expected,(delta,k,l)
print('144 blocs testés : dimensions 3 sur la diagonale paire, 4 ailleurs.')
print('Quatre origines finies, ordre 8 : dimension aveugle totale =',nullity)
print('Origines croisées : rangs concordants pour pi/7 et pi/sqrt(2), ordres 1 à 8.')
