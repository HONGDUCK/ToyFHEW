import numpy as np
from note_include.elem.RLWE import RLWE
from note_include.elem.Ring import Ring
from note_include.utils.gadget_decomposition import gadget_decomposition, format_ring_list, gadget_decomposition_int
from note_include.utils.types import RLWEpctxt, RLWEctxt

class RLWEp:
    def __init__(self, dimension, modulus, std, base, d):
        self.n   = dimension
        self.q   = modulus
        self.std = std
        self.B   = base
        # self.d   = int(np.ceil(np.log(modulus) / np.log(base)))
        self.d   = d
        self.CCrlwe = RLWE(dimension, modulus, std)

    def encrypt(self, msg : Ring, sk : Ring) -> RLWEpctxt:
        ctxts = []
        for i in range(self.d):
            ctxt = self.CCrlwe.encrypt((self.B ** i) * msg, sk)
            ctxts.append(ctxt)

        return ctxts
    
    def decrypt(self, ctxt:RLWEpctxt, sk:Ring) -> list[Ring]:
        res = []
        for rlwe_ctxt in ctxt:
            msg = self.CCrlwe.decrypt(rlwe_ctxt, sk)
            res.append(msg)
        return res
    

    ## There are some problems in this function, do not use this
    def mult_constant(self, ctxts : RLWEpctxt, num : int) -> RLWEctxt:
        '''
            Note that RLWEp multiplication uses only for multiply (huge) constant.
            -> Sure?
        '''
        decomposed_vals = gadget_decomposition_int(num, self.B, self.d)
        zero_coeffs = np.zeros(self.n)
        zero_poly   = Ring(self.n, self.q, zero_coeffs)
        result      = [zero_poly, zero_poly]
        for ctxt, v in zip(ctxts, decomposed_vals):
            if v == 0 : continue
            result = self.CCrlwe.add_ctxt_ctxt(result, self.CCrlwe.mult_ring_cnst(ctxt, v))
        
        return result
    
    def mult_poly(self, ctxts : RLWEpctxt, poly : Ring) -> RLWEctxt:
        decomposed_polys = gadget_decomposition(poly, self.B, self.d)
        zero_coeffs = np.zeros(self.n)
        zero_poly   = Ring(self.n, self.q, zero_coeffs)
        result      = [zero_poly, zero_poly]
        
        for ctxt, d_poly in zip(ctxts, decomposed_polys):
            if all(c == 0 for c in d_poly.coeffs):continue
            tmp = result
            result = self.CCrlwe.add_ctxt_ctxt(tmp, self.CCrlwe.mult_ring_ptxt(ctxt, d_poly))
        
        return result