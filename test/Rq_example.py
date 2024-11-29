import numpy as np
from numpy.polynomial import Polynomial
from include import RLWE, Rq

if __name__ == '__main__':
    n = 8  # power of 2
    q = 67108289  # prime number, q = 1 (mod 2n)
    t = 37  # prime number, t < q
    std = 3  # standard deviation of Gaussian distribution

    coeff1 = np.random.randint(t, size=n)
    coeff2 = np.random.randint(t, size=n)

    m0 = Rq(coeff1, t)
    m1 = Rq(coeff2, t)

    res = m0 * m1

    print(m0)
    print(m1)
    print(res)

    p1 = Polynomial(coeff1)  
    p2 = Polynomial(coeff2)

    f = np.zeros((n+1), dtype=np.int64)  # x^n + 1
    f[0] = f[-1] = 1
    f = np.poly1d(f)

    p = p1 * p2

    p_coeff = p.coef  # Reverse coefficients for descending order
    f_coeff = f.coeffs  # Extract coefficients from poly1d and reverse
    q_coeff, r_coeff = np.polydiv(p_coeff, f_coeff)

    r_mod = np.mod(r_coeff, 37)  # Reverse again for ascending order
    print(r_mod)
