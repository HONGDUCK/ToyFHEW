## Not yet
import numpy as np

def mod_inverse(a, m):
    """Compute the modular inverse of a modulo m using the Extended Euclidean Algorithm."""
    m0 = m
    x0, x1 = 0, 1
    if m == 1:
        return 0
    while a > 1:
        q = a // m
        a, m = m, a % m
        x0, x1 = x1 - q * x0, x0
    return x1 + m0 if x1 < 0 else x1

def prime_factors(n):
    """Return the set of prime factors of n."""
    factors = set()
    # Factor out 2
    while n % 2 == 0:
        factors.add(2)
        n //= 2
    # Factor out odd primes
    i = 3
    while i * i <= n:
        while n % i == 0:
            factors.add(i)
            n //= i
        i += 2
    if n > 2:
        factors.add(n)
    return factors

def primitive_root(q):
    """
    Find a primitive root modulo q.
    Assumes q is prime.
    """
    if q == 2:
        return 1
    phi = q - 1
    factors = prime_factors(phi)
    for g in range(2, q):
        flag = True
        for factor in factors:
            if pow(g, phi // factor, q) == 1:
                flag = False
                break
        if flag:
            return g
    raise ValueError("No primitive root found.")

def nth_roots_of_unity(n, q):
    """
        For a given dimension n and prime modulus q:
          - Compute the inverse of n modulo q.
          - Compute the primitive n-th root of unity, omega.
          - Compute its inverse, omega_inv.
          - List all n-th roots of unity in Z_q.

        Note: n-th roots of unity exist only if n divides q-1.
    """
    # Check necessary condition
    if (q - 1) % n != 0:
        raise ValueError("n does not divide q-1. n-th roots of unity do not exist in Z_q.")
    
    inv_n = mod_inverse(n, q)
    
    # Find a primitive root modulo q
    g = primitive_root(q)
    
    # Compute omega: a primitive n-th root of unity
    exponent = (q - 1) // n
    omega = pow(g, exponent, q)
    
    # Compute the inverse of omega
    omega_inv = mod_inverse(omega, q)
    
    # Compute all n-th roots of unity: {omega^0, omega^1, ..., omega^(n-1)}
    roots = [pow(omega, k, q) for k in range(n)]
    
    return inv_n, roots, omega, omega_inv


def positive_wrapped_ntt(a, omega, q, n):
    # Initialize the result vector with zeros
    hat_a = np.zeros(n, dtype=int)

    # Loop over each j to compute hat_a[j]
    for j in range(n):
        # Summation: hat_a[j] = sum(omega^(ij) * a_i) mod q
        summation = 0
        for i in range(n):
            exponent = (i * j) % (n) # exponent ij modulo n
            summation += pow(omega, exponent, q) * a[i]

        # Take modulo q for the result
        hat_a[j] = summation % q
    return hat_a

def positive_wrapped_intt(hat_a, omega_inv, q, n, n_inv):
    # Initialize the result vector with zeros
    a = np.zeros(n, dtype=int)

    # Loop over each j to compute hat_a[j]
    for i in range(n):
        # Summation: a[i] = sum(omega_inv^{-(ij)} * hat_a[j]) mod q
        summation = 0
        for j in range(n):
            exponent = (i * j) % (n) # exponent 2ij + i modulo 2n
            summation += pow(omega_inv, exponent, q) * hat_a[j]
            
        # Take modulo q for the result
        a[i] = (n_inv * summation) % q
            
    return a