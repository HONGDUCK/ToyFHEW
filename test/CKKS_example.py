from include.CKKSEncoder import CKKSEncoder
from include.RLWE import RLWE
from include.core.Ring import Rq
import numpy as np
from numpy.polynomial import Polynomial

if __name__ == '__main__':
    n = 16  # power of 2
    scale = 64
    q = 67108289  # prime number, q = 1 (mod 2n)
    t = 4095  # prime number, t < q
    std = 3  # standard deviation of Gaussian distribution

    encoder = CKKSEncoder(n, scale)
    
    z1 = np.random.randint(t, size=n // 4)
    z2 = np.random.randint(t, size=n // 4)
    p1 = encoder.encode(z1)
    p2 = encoder.encode(z2)
    m1 = Rq(p1.coef, t)
    m2 = Rq(p2.coef, t)
    p = m1 + m2

    print(z1)
    print(z2)

    print(p1.coef)
    print(p2.coef)

    print(p1)
    print(p2)    

    print(m1)
    print(m2)

    print(p)

    decoded_z = encoder.decode(Polynomial(p.poly))
    print(decoded_z.real)