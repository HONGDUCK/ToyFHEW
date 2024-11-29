import numpy as np
from numpy.polynomial import Polynomial
from include.core import CRTRq
from include.core.ModuliChainGen import generate_pairwise_coprime
from math import prod

if __name__ == '__main__':
    n = 8  # power of 2
    q = 67108289  # prime number, q = 1 (mod 2n)
    t = 2063  # prime number, t < q
    std = 3  # standard deviation of Gaussian distribution
    d = 5 # Number of coprime pairs
    
    moduli_chain = generate_pairwise_coprime(d)
    Q = prod(moduli_chain)

    f = np.poly1d([1] + [0] * (n - 1) + [1])  # x^N + 1

    coeff0 = np.random.randint(t, size=n)
    coeff1 = np.random.randint(t, size=n)

    # print("Poly 0 : ", coeff0)
    # print("Poly 1 : ", coeff1)

    m0 = CRTRq(np.poly1d(coeff0), moduli_chain)
    m1 = CRTRq(np.poly1d(coeff1), moduli_chain)

    print(m0)
    print(m1)

    m_add = m0 + m1
    print("Origin     Add Poly : ", coeff0 + coeff1)
    print("Recomposed Add Poly : ", CRTRq.crt_recompose(m_add))

    m_mult  = m0 * m1
    m_mult2 = CRTRq.fftmul(m0, m1)
    _, r = np.polydiv(np.polymul(coeff0, coeff1), f)
    print("Origin     Mult Poly : ", np.mod(r, Q))
    print("Recomposed Mult Poly : ", CRTRq.crt_recompose(m_mult))
    print("Recomposed Mult with FFT Poly : ", CRTRq.crt_recompose(m_mult2))
