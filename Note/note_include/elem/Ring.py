import numpy as np

def pad_coeffs(coeffs, n):
    """Pads coeffs with zeros to length n if necessary."""
    if len(coeffs) < n:
        return np.pad(coeffs, (0, n - len(coeffs)), 'constant')
    return coeffs

class Ring:
    def __init__(self, dimension, modulus, coeffs):
        self.n = dimension
        self.q = modulus
        self.coeffs = [int(coeff % modulus) for coeff in coeffs]

    def nega_conv(self, other):
        result = []
        poly1 = pad_coeffs(self.coeffs, self.n)
        poly2 = pad_coeffs(other.coeffs, self.n)
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

        add_coeffs = [a + b for a,b in zip(self.coeffs, other.coeffs)]
        return Ring(self.n, self.q, add_coeffs)
    
    def ring_sub_q(self, other):
        assert self.q == other.q
        assert self.n == other.n

        sub_coeffs = [a - b for a,b in zip(self.coeffs, other.coeffs)]
        return Ring(self.n, self.q, sub_coeffs)

    def __mul__(self, other):
        return self.nega_conv(other)
    
    def __rmul__(self, integer : int):
        coeffs = [c * integer for c in self.coeffs]
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
            return f"({poly_str} | n={self.n}, q={self.q})"

    def __getitem__(self, index):
        return self.coeffs[index]  # Allow indexing into the coefficient array
    

