import numpy as np
from note_include.elem.Ring   import Ring
from note_include.elem.LWE    import LWE
from note_include.elem.RLWE   import RLWE
from note_include.elem.RLWEp  import RLWEp
from note_include.elem.RGSW   import RGSW
from note_include.utils.types import RGSWctxt, RLWEctxt, RLWEpctxt, LWEctxt
# from note_include.FHEW        import FHEW
from joblib import Parallel, delayed

class DM:
    # Method -> DM, CGGI, LMKCDEY
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

        self.RGSW_CC  = RGSW (dimension, modulus, s_std, e_std, base_gd, self.d_g)

    # -------------------------------- KeyGen ---------------------------- #

    def keygen(self, s:list[int], s_ring:Ring):
        brk = self.BRKgen_parallel(s, s_ring) # Blind rotation key generation.
        return brk

    # # Its too slow, not used. # It has a little bit of problems.
    # def BRKgen(self, s:list[int], s_ring:Ring) -> list[list[RGSWctxt]]:
    #     brk = []
    #     for _s in s:
    #         matrix = []
    #         for v in range(self.B):
    #             row = []
    #             for j in range(self.d_g):
    #                 coef = np.zeros(self.N)
    #                 tmp = (v * (self.B**j) * _s) % self.q
    #                 coef[(self.N//2) - int(tmp)] = 1
    #                 poly = Ring(self.N, self.Q, coef)
    #                 row.append(self.RGSW_CC.encrypt(poly, s_ring))
    #             matrix.append(row)
    #         brk.append(matrix)
    #     return brk

    def BRKgen_parallel(self, s:list[int], s_ring:Ring):
        return Parallel(n_jobs=-1)(
            delayed(generate_rgsw_matrix)(_s, self.B, self.d_g, self.q, self.N, self.Q, self.RGSW_CC, s_ring)
            for _s in s
        )
    
    def Blindrotation(self, ctxt_operand:LWEctxt, acc:RLWEctxt, brk:list[list[RGSWctxt]]) -> RLWEctxt:
        a,_ = ctxt_operand
        
        # Blind rotation
        for i, _a in enumerate(a):
            if(_a == 0):continue
            for j in range(self.d_g):
                v = _a % self.B
                _a = _a // self.B

                acc = self.RGSW_CC.mult_rlwe(acc, brk[i][v][j])

        return acc

# -------------------------------- Blind rotation key generation : parallel ---------------------------- #
def generate_rgsw_row(v, B, d_g, _s, q, N, Q, RGSW_CC, s_ring):
    row = []
    for j in range(d_g):
        _as = int((v * (B**j) * _s) % q)

        reduced = False
        if (_as >= N) :
            _as -= N
            reduced = True

        if _as == 0 and reduced == True:
            monomial = np.zeros(N)
            monomial[0] = -1
            monomial = Ring(N, Q, monomial)
            row.append(RGSW_CC.encrypt(monomial, s_ring))
            continue

        elif _as == 0 and reduced == False:
            monomial = np.zeros(N)
            monomial[0] = 1
            monomial = Ring(N, Q, monomial)
            row.append(RGSW_CC.encrypt(monomial, s_ring))
            continue

        if reduced == False:
            monomial = np.zeros(N)
            monomial[(N)-(_as)] = -1
            monomial = Ring(N, Q, monomial)
        else:
            monomial = np.zeros(N)
            monomial[(N)-(_as)] = 1
            monomial = Ring(N, Q, monomial)

        row.append(RGSW_CC.encrypt(monomial, s_ring))
    return row

def generate_rgsw_matrix(_s, B, d_g, q, N, Q, RGSW_CC, s_ring):
    return Parallel(n_jobs=-1)(
        delayed(generate_rgsw_row)(v, B, d_g, _s, q, N, Q, RGSW_CC, s_ring)
        for v in range(B)
    )