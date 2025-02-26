import numpy as np
from note_include.elem.Ring           import Ring
from note_include.elem.LWE            import LWE
from note_include.elem.RLWE           import RLWE
from note_include.elem.RLWEp          import RLWEp
from note_include.elem.RGSW           import RGSW
from note_include.utils.types         import RGSWctxt, RLWEctxt, RLWEpctxt, LWEctxt
from note_include.blindrotations.DM   import DM
from note_include.blindrotations.CGGI import CGGI
from joblib import Parallel, delayed

class FHEW:
    # Method -> DM, CGGI, LMKCDEY
    def __init__(self, lwe_dimension, 
                       lwe_modulus, 
                       dimension, 
                       modulus, 
                       s_std, 
                       e_std, 
                       base_gd, 
                       base_ks, 
                       method = 'DM'):
        
        # it should be satisfy q = 2N
        assert lwe_modulus == dimension * 2
        assert method == "DM" or method == "CGGI" or method == "LMKCDEY", "Method should be one of [DM, CGGI, LMKCDEY]."
        assert s_std == "Gaussian" or s_std == "Binary", "Secret key distribution should be one of [Gaussian, Binary]."
        
        if method == "CGGI":
            assert method == "CGGI" and s_std == "Binary", "We impl CGGI blind rotation only for binary key distribution."

        self.n        = lwe_dimension
        self.q        = lwe_modulus
        self.N        = dimension
        self.Q        = modulus
        self.s_std    = s_std
        self.e_std    = e_std
        self.B        = base_gd
        self.d_g      = int(np.ceil(np.log(lwe_modulus) / np.log(base_gd)))
        self.B_ks     = base_ks
        self.d_ks     = int(np.ceil(np.log(modulus) / np.log(base_ks)))

        self.LWE_CC   = LWE  (lwe_dimension, lwe_modulus, s_std, e_std)
        self.LWEQ_CC  = LWE  (lwe_dimension, modulus,     s_std, e_std)
        self.RLWE_CC  = RLWE (dimension,     modulus,     s_std, e_std)
        # self.RLWEp_CC = RLWEp(dimension,     modulus,     std, base_gd, self.d_g)
        self.RGSW_CC  = RGSW (dimension,     modulus,     s_std, e_std, base_gd, self.d_g)

        if method == 'DM':
            self.method =   DM(lwe_dimension, lwe_modulus, dimension, modulus, s_std, e_std, base_gd) 
        elif method== 'CGGI':
            self.method = CGGI(lwe_dimension, lwe_modulus, dimension, modulus, s_std, e_std, base_gd) 

    # -------------------------------- KeyGen ---------------------------- #
    def keygen(self): # Further impl
        s              = self.LWE_CC.keygen()
        s_ring,pk_ring = self.RLWE_CC.keygen()

        self.brk       = self.method.keygen(s, s_ring)
        self.ksk       = self.KSKgen(s, s_ring)
        self.pk        = pk_ring

        return s, s_ring
    
    def KSKgen(self, lwe_sk:list[int], rlwe_sk:Ring):
        ksk = []
        for _s in rlwe_sk.coeffs:
            matrix = []
            for v in range(self.B_ks):
                row = []
                for j in range(self.d_ks):
                    row.append(self.LWEQ_CC.encrypt((v * (self.B_ks**j) * _s) % self.Q, lwe_sk))
                matrix.append(row)
            ksk.append(matrix)    
        return ksk
    
    # ------------------------------------------------------------------- #

    # --------------------Encryption and decryption---------------------- #

    def encrypt(self, m:int, s:list[int]) -> LWEctxt:
        assert m == 0 or m == 1, "message must be 0 or 1"
        return self.LWE_CC.encrypt(m * (self.q // 4), s)
    
    def decrypt(self, ctxt : LWEctxt, s : list[int]) -> int:
        m = self.LWE_CC.decrypt(ctxt, s)
        if m <= self.q // 8 or m > self.q * 5 // 8 :
            m = 0
        else:
            m = 1
        return m

    def acc_init(self, ctxt:LWEctxt, gate = "AND") -> RLWEctxt:
        _, b = ctxt
        q0, q1, q2, q3 = mapping_function(self.q, gate)

        swap1 = False
        swap2 = False
        if q0 > q1 : swap1 = True
        if q2 > q3 : swap2 = True    

        acc = np.zeros(self.N)

        for i in range(self.q//2):
            b = (b - 1) % self.q

            if swap1 == False:
                if b >= q0 and b < q1:
                    acc[i] = self.Q // 8
            
            elif swap1 == True:
                if b > q0 or b <= q1:
                    acc[i] = self.Q // 8

            if swap2 == False:
                if b >= q2 and b < q3:
                    acc[i] = -self.Q // 8
            
            elif swap2 == True:
                if b > q2 or b <= q3:
                    acc[i] = -self.Q // 8
            
        poly_acc = Ring(self.N, self.Q, acc)
        acc = self.RLWE_CC.pk_encrypt(poly_acc, self.pk)

        return acc
    
    # ------------------------------------------------------------------- #

    # --------------------Ket Switch and Mod Switch---------------------- #

    def KeySwitch(self, ctxt:LWEctxt) -> LWEctxt:
        a, b = ctxt

        tmp_ct = (np.zeros(self.n), b)
        for i, _a in enumerate(a):
            tmpa = int(_a)
            for j in range(self.d_ks):
                v  = tmpa % self.B_ks
                tmpa = tmpa // self.B_ks

                tmp_ct = self.LWE_CC.sub(tmp_ct, self.ksk[i][v][j])

        return tmp_ct
    
    def ModSwitch(self, ctxt:RLWEctxt) -> LWEctxt:
        a, b = ctxt
        a_   = [(np.round(coef*self.q/self.Q)) % self.q for coef in a]
        b_   = (np.round(b*self.q/self.Q)) % self.q

        return (a_, b_)    

    # ---------------------------------------------------------------------- #

    # --------------------Bootstrapping(BlindRotation)---------------------- #
    
    def LWEextract(self, ctxt:RLWEctxt) -> LWEctxt:
        a, b = ctxt
        b_0  = int(b.coeffs[0] + (self.Q // 8))                                                 # extract constant term
        a_coeffs = [a.coeffs[0]] + [(-x) % a.q for x in reversed(a.coeffs[1:])] # suit for negacyclic works

        return (a_coeffs, b_0)

    def evalBin(self, ctxt1:LWEctxt, ctxt2:LWEctxt, gate="AND"):
        assert gate == "AND" or gate == "OR" or gate == "XOR" or gate == "NAND" or gate == "NOR" or gate == "XNOR", "Gate should be one of [AND, OR, XOR, NAND, NOR, XNOR]."

        # Operate homomorphic addition and saclar mmultiplication
        if gate == "XOR" or gate == "XNOR":
            ctxt_operand =  self.LWE_CC.mult_cnst(self.LWE_CC.add(ctxt1, ctxt2), 2)
        else:
            ctxt_operand =  self.LWE_CC.add(ctxt1, ctxt2)

        acc = self.acc_init(ctxt_operand, gate)

        rotated_acc = self.method.Blindrotation(ctxt_operand, acc, self.brk)

        # LWE extraction
        extracted_lwe = self.LWEextract(rotated_acc)

        # KeySwitch
        key_switched_lwe = self.KeySwitch(extracted_lwe)

        # ModSwitch # Too many noise so we do not use in toy example
        # mod_switched_lwe = self.ModSwitch(key_switched_lwe)

        return key_switched_lwe
    
    def evalNOT(self, ct:LWEctxt) -> LWEctxt:
        a, b = ct
        b_   = b - self.q//4
        a_   = [int(-aa % self.q )for aa in a]

        return (a_, -b_)
    

# ------------------ Mapping Function ------------------ #

def mapping_function(q:int, gate = "AND"):
    assert gate == "AND" or gate == "OR" or gate == "XOR" or gate == "NAND" or gate == "NOR" or gate == "XNOR", "Gate should be one of [AND, OR, XOR, NAND, NOR, XNOR]."

    if   gate == "AND":
        q0, q1, q2, q3 = (3 * q / 8), (7 * q / 8), (-q / 8) % q, (3 * q / 8)
    elif gate == "NAND":
        q0, q1, q2, q3 = (-q / 8) % q, (3 * q / 8), (3 * q / 8), (7 * q / 8)
    elif gate == "OR":
        q0, q1, q2, q3 = (q / 8) % q, (5 * q / 8), (-3 * q / 8) % q, (q / 8)
    elif gate == "NOR":
        q0, q1, q2, q3 = (-3 * q / 8) % q, (q / 8), (q / 8) % q, (5 * q / 8)
    elif gate == "XOR":
        q0, q1, q2, q3 = (q / 4), (3 * q / 4), (-q / 4) % q, (q / 4)
    elif gate == "XNOR":
        q0, q1, q2, q3 = (-q / 4) % q, (q / 4), (q / 4), (3 * q / 4)
    
    return q0, q1, q2, q3
