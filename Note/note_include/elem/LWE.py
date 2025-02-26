import numpy as np
from note_include.utils.types import LWEctxt
from note_include.utils.noise_generator import discrete_uniform, discrete_gaussian

class LWE:
    def __init__(self, dimension, modulus, s_std, e_std):
        self.n     = dimension
        self.q     = modulus
        self.s_std = s_std
        self.e_std = e_std

    def keygen(self):
        if self.s_std == "Gaussian":
            s = [int(_s) for _s in discrete_gaussian(self.n, self.q, std=3.2).coeffs] # default
        elif self.s_std == "Binary":
            s = [int(_s) for _s in discrete_uniform(self.n, 2).coeffs]

        return s

    def encrypt(self, msg : int, sk : list[int]) -> LWEctxt:
        a   = np.random.randint(0, self.q-1, self.n)
        e   = np.round(self.e_std * np.random.rand()) % self.q
        as_ = sum(x * y for x, y in zip(a, sk))


        b   = (as_ + msg + e) % self.q
        return (a, b)

    def decrypt(self, ctxt : LWEctxt, sk : list[int]) -> int:
        a, b = ctxt
        as_  = sum(x * y for x, y in zip(a, sk))

        result = b - as_
        return result % self.q
    
    def add(self, c1:LWEctxt, c2:LWEctxt) -> LWEctxt:
        a1, b1 = c1
        a2, b2 = c2

        return ((a1 + a2) % self.q, (b1 + b2) % self.q)
    
    def sub(self, c1:LWEctxt, c2:LWEctxt) -> LWEctxt:
        a1, b1 = c1
        a2, b2 = c2

        return ((a1 - a2) % self.q, (b1 - b2) % self.q)
    
    def mult_cnst(self, ctxt:LWEctxt, d:int) -> LWEctxt:
        a, b = ctxt
        return ((a*d) % self.q, (b*d) % self.q)
