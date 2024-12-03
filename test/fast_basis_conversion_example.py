import numpy as np
from include.core import CRTRq
from include.core.ModuliChainGen import generate_pairwise_coprime

if __name__ == '__main__':
    n = 8  # power of 2
    q = 67108289  # prime number, q = 1 (mod 2n)
    t = 2063  # prime number, t < q
    std = 3  # standard deviation of Gaussian distribution
    d = 5 # Number of coprime pairs
    
    moduli_chain = generate_pairwise_coprime(d * 2)
    moduli_chain1 = moduli_chain[:5]
    moduli_chain2 = moduli_chain[5:]

    print("Moduli Chain1 : ", moduli_chain1)
    print("Moduli Chain2 : ", moduli_chain2)


    f = np.poly1d([1] + [0] * (n - 1) + [1])  # x^N + 1
    coeff0 = np.random.randint(t, size=n)
    print("Poly 0 : ", coeff0)

    m0 = CRTRq(np.poly1d(coeff0), moduli_chain1)
    m1 = CRTRq(np.poly1d(coeff0), moduli_chain2)

    print(m0)
    print(m1)

    CRTRq.fast_basis_conversion(m0, moduli_chain2)

    print(m0)