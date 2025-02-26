import numpy as np
from note_include.elem.Ring  import Ring
from note_include.elem.RLWE  import RLWE
from note_include.elem.RLWEp import RLWEp
from note_include.utils.types import RGSWctxt, RLWEctxt

class RGSW:
    def __init__(self, 
                 dimension, 
                 modulus, 
                 s_std,
                 e_std, 
                 base, 
                 d):
        self.n       = dimension
        self.q       = modulus
        self.s_std   = s_std
        self.e_std   = e_std
        self.B       = base
        self.d       = d
        self.CCrlwe  = RLWE(dimension, modulus, s_std, e_std)
        self.CCrlwep = RLWEp(dimension, modulus, s_std, e_std, base, d)

    def encrypt(self, msg : Ring, sk : Ring) -> RGSWctxt:
        e0 = self.CCrlwep.encrypt(-1 * sk * msg, sk)
        e1 = self.CCrlwep.encrypt(msg, sk)

        return [e0, e1]
    
    def mult_rlwe(self, rlwe_ctxt : RLWEctxt, rgsw_ctxt : RGSWctxt) -> RLWEctxt:
        a, b     = rlwe_ctxt
        ct0, ct1 = rgsw_ctxt

        tmp1 = self.CCrlwep.mult_poly(ct0, a) # -> Return RLWEctxt
        tmp2 = self.CCrlwep.mult_poly(ct1, b) # -> Return RLWEctxt
        res  = self.CCrlwe.add_ctxt_ctxt(tmp1, tmp2)

        return res
