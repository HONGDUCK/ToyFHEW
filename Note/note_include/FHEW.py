import numpy as np
from note_include.elem.Ring  import Ring
from note_include.elem.RLWE  import RLWE
from note_include.elem.RLWEp import RLWEp
from note_include.utils.types import RGSWctxt, RLWEctxt, RLWEpctxt

class FHEW:
    # Method -> DM, CGGI, LMKCDEY
    def __init__(self, dimension, modulus, std, base, method):
        self.n       = dimension
        self.q       = modulus
        self.std     = std
        self.B       = base
        self.d       = int(np.ceil(np.log(modulus) / np.log(base)))

        self.RLWE_CC  = RLWE(dimension, modulus, std)
        self.RLWEp_CC = RLWEp(dimension, modulus, std, base)
        self.RGSW_CC  = RLWE(dimension, modulus, std)

    