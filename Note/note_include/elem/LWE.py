import numpy as np
from note_include.utils.types import LWEctxt


class LWE:
    def __init__(self, dimension, modulus, std):
        self.n   = dimension
        self.q   = modulus
        self.std = std

    def keygen(self):
        s = np.round(3.2 * np.random.randn(self.n)) % self.q
        # s = np.round(self.std * np.random.randn(self.n)) % self.q
        return s

    def encrypt(self, msg : int, sk : list[int]) -> LWEctxt:
        a   = np.random.randint(0, self.q-1, self.n)
        e   = np.round(self.std * np.random.rand()) % self.q
        as_ = np.sum(a * sk)

        b   = (as_ + msg + e) % self.q
        return (a, b)

    def decrypt(self, ctxt : LWEctxt, sk : list[int]) -> int:
        a, b = ctxt
        as_  = np.sum(a * sk)

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
