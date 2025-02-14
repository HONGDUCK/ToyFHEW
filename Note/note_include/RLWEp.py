import numpy as np
from note_include.RLWE import RLWE
from note_include.Ring import Ring
from note_include.utils.gadget_decomposition import gadget_decomposition
from note_include.utils.types import RLWEpctxt, RLWEctxt

class RLWEp:
    def __init__(self, dimension, modulus, std, base):
        self.n   = dimension
        self.q   = modulus
        self.std = std
        self.B   = base
        self.d   = int(np.ceil(np.log(modulus) / np.log(base)))
        self.CCrlwe = RLWE(dimension, modulus, std)

    def encrypt(self, msg : Ring, sk : Ring) -> RLWEpctxt:
        ctxts = []
        for i in range(self.d):
            ctxt = self.CCrlwe.encrypt((self.B ** i) * msg, sk)
            ctxts.append(ctxt)

        return ctxts
    
    def mult_constant(self, ctxts : RLWEpctxt, num : int) -> RLWEctxt:
        '''
            Note that RLWEp multiplication uses only for multiply (huge) constant.
            -> Sure?
        '''
        num_coeff = np.zeros(self.n)
        num_coeff[0] = num
        num_ring = Ring(self.n, self.q, num_coeff)

        decomposed_polys = gadget_decomposition(num_ring, self.B, self.d)
        zero_coeffs = np.zeros(self.n)
        zero_poly   = Ring(self.n, self.q, zero_coeffs)
        result      = [zero_poly, zero_poly]
        for ctxt, poly in zip(ctxts, decomposed_polys):
            result = self.CCrlwe.add_ctxt_ctxt(result, self.CCrlwe.mult_ring_ptxt(ctxt, poly))
        
        return result
    
    def mult_poly(self, ctxts : RLWEpctxt, poly : Ring) -> RLWEctxt:
        decomposed_polys = gadget_decomposition(poly, self.B, self.d)
        zero_coeffs = np.zeros(self.n)
        zero_poly   = Ring(self.n, self.q, zero_coeffs)
        result      = [zero_poly, zero_poly]
        for ctxt, poly in zip(ctxts, decomposed_polys):
            result = self.CCrlwe.add_ctxt_ctxt(result, self.CCrlwe.mult_ring_ptxt(ctxt, poly))
        
        return result