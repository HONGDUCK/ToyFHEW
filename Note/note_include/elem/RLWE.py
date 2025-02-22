from note_include.elem.Ring import Ring
from note_include.utils.noise_generator import discrete_gaussian, discrete_uniform
from note_include.utils.types import RLWEctxt

class RLWE:
    def __init__(self, dimension, modulus, std):
        # assert np.log2(dimension) == int(np.log2(dimension)) # power-of-two dimension
        self.n   = dimension
        self.q   = modulus
        self.std = std

    def keygen(self):
        s = discrete_gaussian(self.n, self.q, std=3.2)     # default
        e = discrete_gaussian(self.n, self.q, std=self.std)

        a0 = discrete_uniform(self.n, self.q)
        a1 = (a0 * s + e)

        return (s, (a0, a1)) # (secret key, public key)
    
    # Only for symmetric encryption
    def encrypt(self, msg : Ring, sk : Ring) -> RLWEctxt:
        # e = discrete_gaussian(self.n, self.q, std=self.std) # Noise
        e = discrete_gaussian(self.n, self.q, mean = 0, std = 0) # Noise
        a = discrete_uniform(self.n, self.q, 0, self.q)                # Random Num


        b = a * sk + msg + e
        return (a, b)
    
    def pk_encrypt(self, msg:Ring, pk:RLWEctxt) -> RLWEctxt:
        a0, a1 = pk
        return (a0, a1 + msg)
    
    def decrypt(self, ctxt : RLWEctxt, sk : Ring) -> Ring:
        a, b = ctxt
        return b - (a * sk)
    
    def add_ctxt_ctxt(self, ct1 : RLWEctxt, ct2 : RLWEctxt) -> RLWEctxt:
        a1, b1 = ct1
        a2, b2 = ct2

        return [a1 + a2, b1 + b2]
    
    def add_ctxt_ptxt(self, ct : RLWEctxt, Z_ring : Ring) -> RLWEctxt:
        a, b = ct
        return [a, b + Z_ring]

    def mult_ring_ptxt(self, ct : RLWEctxt, Z_ring : Ring) -> RLWEctxt:
        '''
            Note that multiplier should be small to keep the noise be small.
        '''
        a, b = ct
        return [a * Z_ring, b * Z_ring]
    
    def mult_ring_cnst(self, ct : RLWEctxt, const : int) -> RLWEctxt:
        a, b = ct
        return [const * a, const * b]
    