import numpy as np
from note_include.elem.Ring   import Ring
from note_include.elem.LWE    import LWE
from note_include.elem.RLWE   import RLWE
from note_include.elem.RLWEp  import RLWEp
from note_include.elem.RGSW   import RGSW
from note_include.utils.types import RGSWctxt, RLWEctxt, RLWEpctxt, LWEctxt
from joblib import Parallel, delayed

class CGGI:
    def __init__(self, lwe_dimension, lwe_modulus, dimension, modulus, s_std, e_std, base_gd):
        assert lwe_modulus == dimension * 2
        self.n        = lwe_dimension
        self.q        = lwe_modulus
        self.N        = dimension
        self.Q        = modulus
        self.s_std    = s_std
        self.e_std    = e_std
        self.B        = base_gd
        self.d_g      = int(np.ceil(np.log(lwe_modulus) / np.log(base_gd)))

        self.RLWE_CC  = RLWE (dimension,     modulus,     s_std, e_std)
        self.RGSW_CC  = RGSW (dimension, modulus, s_std, e_std, base_gd, self.d_g)

    # -------------------------------- KeyGen ---------------------------- #

    def keygen(self, s:list[int], s_ring:Ring):
        brk = self.BRKgen(s, s_ring) # Blind rotation key generation.
        return brk
    
    def BRKgen(self, s:list[int], s_ring:Ring) -> list[RGSWctxt]: # RGSW(X^{-s_i}) for i \in [0, n)
        brk = []
        for _s in s:
            monomial = np.zeros(self.N)
            monomial[0] = _s
            monomial = Ring(self.N, self.Q, monomial)
            brk.append(self.RGSW_CC.encrypt(monomial, s_ring))
        return brk
    
    def Blindrotation(self, ctxt_operand:LWEctxt, acc:RLWEctxt, brk:list[RGSWctxt]) -> RLWEctxt:
        a,_ = ctxt_operand
        N   = int(self.N)        

        # Blind rotation
        for i, _a in enumerate(a):
            if(_a == 0):continue
            monomial = np.zeros(self.N)

            reduced = False
            if _a >= self.N:
                _a -= self.N
                reduced = True

            if _a == 0 and reduced == True:
                monomial[0]    = -2
            elif _a == 0 and reduced == False:
                continue
            elif reduced == True:           
                monomial[N-_a] =  1
                monomial[0]    = -1            
            elif reduced == False:          
                monomial[N-_a] = -1
                monomial[0]    = -1

            monomial     = Ring(self.N, self.Q, monomial)
            tmp = acc
            acc = self.RGSW_CC.mult_rlwe(acc, brk[i])
            acc = self.RLWE_CC.mult_ring_ptxt(acc, monomial)
            acc = self.RLWE_CC.add_ctxt_ctxt(tmp, acc)

        return acc
    