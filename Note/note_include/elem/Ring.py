import numpy as np

class Ring:
    '''
        ring : Zq[x]/(x^n+1)
    '''
    def __init__(self, dimension, modulus, coeffs):
        self.n = dimension
        self.q = modulus
        self.coeffs = [coeff % modulus for coeff in coeffs]
        divisor = np.zeros((dimension+1), dtype=np.int64)  # x^n + 1
        divisor[0] = divisor[-1] = 1
        self.divisor = divisor
        # self.inv_n, self.roots, self.omega, self.omega_inv = nth_roots_of_unity(dimension, modulus)

    def pad_coeffs(self, coeffs):
        """Pads coeffs with zeros to length n if necessary."""
        if len(coeffs) < self.n:
            return np.pad(coeffs, (0, self.n - len(coeffs)), 'constant')
        return coeffs

    def NTT_like_PWC(self, other):
        assert self.omega == other.omega
        assert self.omega_inv == other.omega_inv
        assert self.q == other.q
        assert self.n == other.n

        ntt_poly1 = self.positive_wrapped_ntt(self.coeffs, self.omega, self.q, self.n)
        ntt_poly2 = self.positive_wrapped_ntt(other.coeffs, other.omega, other.q, other.n)

        ntt_res = np.zeros(self.n)
        for i in range(self.n):
            ntt_res[i] = (ntt_poly1[i] * ntt_poly2[i]) % self.q

        poly_nega_conv = self.positive_wrapped_intt(ntt_res, self.omega_inv, self.q, self.n, self.inv_n)
        return Ring(self.n, self.q, poly_nega_conv)

    # def negative_ring_mult_q(self, other):
    #     assert self.q == other.q
    #     assert self.n == other.n # !! watchout
        
    #     prod = np.polymul(self.coeffs[::-1], other.coeffs[::-1])
    #     _, r = np.polydiv(prod, self.divisor)
    #     r = np.mod(r, self.q)

    #     return Ring(self.n, self.q, r[::-1])
    
    def negative_ring_mult_q(self, other):
        result = []
        poly1 = self.pad_coeffs(self.coeffs)
        poly2 = self.pad_coeffs(other.coeffs)
        # poly1 = self.coeffs[::-1]
        # poly2 = self.coeffs[::-1]
        for k in range(self.n):
            v = 0
            for i in range(k+1):
                v += poly1[i] * poly2[k-i]
            for i in range(k+1, self.n):
                v -= poly1[i] * poly2[(k + self.n - i)]

            result.append(v)
        return Ring(self.n, self.q, result)
    
    def ring_add_q(self, other):
        assert self.q == other.q
        assert self.n == other.n

        add_coeffs = [int((a + b) % self.q) for a,b in zip(self.coeffs, other.coeffs)]
        return Ring(self.n, self.q, add_coeffs)
    
    def ring_sub_q(self, other):
        assert self.q == other.q
        assert self.n == other.n

        sub_coeffs = [int((a - b) % self.q) for a,b in zip(self.coeffs, other.coeffs)]
        return Ring(self.n, self.q, sub_coeffs)

    def __mul__(self, other):
        return self.negative_ring_mult_q(other)
    
    def __rmul__(self, integer : int):
        coeffs = [(c * integer) % self.q for c in self.coeffs]
        return Ring(self.n, self.q, coeffs)
    
    def __mod__(self, integer):
        coeffs = [co % integer for co in self.coeffs]
        return Ring(self.n, self.q, coeffs)
    
    def __add__(self, other):
        return self.ring_add_q(other)
    
    def __sub__(self, other):
        return self.ring_sub_q(other)

    def __repr__(self, spilt=False):
        terms = []
        for i, coef in enumerate(self.coeffs):
            # if coef != 0:
                if i == 0:
                    terms.append(f"{coef}")
                elif i == 1:
                    terms.append(f"{coef}x")
                else:
                    terms.append(f"{coef}x^{i}")
    
        poly_str = " + ".join(terms) if terms else "0"
    
        # Split into two lines
        if spilt == True:
            mid = len(terms) // 2  # Find the middle index
            first_half = " + ".join(terms[:mid])
            second_half = " + ".join(reversed(terms[mid:]))
            second_half = " + ".join(terms[mid:])
            return f"R(n={self.n}, q={self.q}, coeffs=\n  {first_half}\n  {second_half})"
    
        else:
            return f"R(n={self.n}, q={self.q}, coeffs= {poly_str})"

    def __getitem__(self, index):
        return self.coeffs[index]  # Allow indexing into the coefficient array
    



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