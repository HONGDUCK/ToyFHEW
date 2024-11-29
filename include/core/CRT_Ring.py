import numpy as np
from math import prod

class CRTRq:
    def __init__(self, input_poly, M=None):
        """
        Initialize a CRTRq instance.

        Args:
            input_poly: Either a single np.poly1d object or a list of np.poly1d objects
                        corresponding to CRT components.
        Params:
            moduli
            Q
            N
            f
            crt_poly
        """

        if M is not None:
            self.moduli = M
            self.Q = prod(M)
        elif self.moduli is None and M is None:
            raise ValueError("Moduli must be initialized with CRTRq.__moduli_init__(d) before instantiation.")
        
        if isinstance(input_poly, np.poly1d):
            # Single polynomial: perform CRT decomposition
            self.N = input_poly.order + 1
            self.f = np.poly1d([1] + [0] * (self.N - 1) + [1])  # x^N + 1
            self.crt_poly = self.crt_decompose(input_poly)
        elif isinstance(input_poly, list) and all(isinstance(p, np.poly1d) for p in input_poly):
            # List of CRT components
            if len(input_poly) != len(self.moduli):
                raise ValueError("Input list must match the number of moduli.")
            self.N = max(poly.order + 1 for poly in input_poly)
            self.f = np.poly1d([1] + [0] * (self.N - 1) + [1])  # x^N + 1
            self.crt_poly = [np.mod(p, m) for p, m in zip(input_poly, self.moduli)]
        else:
            raise ValueError("Input must be either a single np.poly1d or a list of np.poly1d objects.")

    def __repr__(self):
        """
        String representation of the CRT polynomials and the recomposed polynomial.
        """
        result = "------------------------------------------------------------------------\n"
        template = "CRT Poly [{}] : {} (mod {})"
        template2 = "Recomposed Poly : {} (mod {})"
        for idx, (poly, m) in enumerate(zip(self.crt_poly, self.moduli)):
            result += template.format(idx, poly, m) + "\n"
        result += template2.format(self.crt_recompose(), self.Q) + "\n"
        return result

    def crt_decompose(self, poly):
        """
        Decompose a polynomial into its CRT components.

        Args:
            poly: np.poly1d to be decomposed.

        Returns:
            List of np.poly1d objects modulo each modulus.
        """
        return [np.mod(poly, m) for m in self.moduli]

    def crt_recompose(self):
        """
        Recompose the polynomial from its CRT components.

        Returns:
            np.poly1d object modulo Q.
        """
        poly = np.poly1d([0])
        for cp, modulus in zip(self.crt_poly, self.moduli):
            Ni = self.Q // modulus
            Mi = pow(Ni, -1, modulus)
            cp = np.polymul(cp, Ni)
            cp = np.polymul(cp, Mi)
            poly = np.polyadd(poly, cp)
        return np.mod(poly, self.Q)

    def __add__(self, other):
        """
        Add two CRTRq objects.

        Args:
            other: Another CRTRq object.

        Returns:
            A new CRTRq object representing the sum.
        """
        if not isinstance(other, CRTRq):
            raise TypeError("Addition is supported only between two CRTRq objects.")
        elif other.moduli != self.moduli:
            raise TypeError("Moduli chain should be same")
        
        result = [
            np.poly1d(np.mod(np.polyadd(c1, c2), m))
            for c1, c2, m in zip(self.crt_poly, other.crt_poly, self.moduli)
        ]
        return CRTRq(result, self.moduli)

    def __mul__(self, other):
        """
        Multiply two CRTRq objects.

        Args:
            other: Another CRTRq object.

        Returns:
            A new CRTRq object representing the product.
        """


        if isinstance(other, CRTRq):
            if other.moduli != self.moduli:
                raise TypeError("Moduli chain should be same")
            
            result = []
            for c1, c2, m in zip(self.crt_poly, other.crt_poly, self.moduli):
                product = np.polymul(c1, c2)
                _, remainder = np.polydiv(product, self.f)
                result.append(np.poly1d(np.mod(remainder, m)))
            return CRTRq(result, self.moduli)
        elif isinstance(other, int):
            result = []
            for p1, m in zip(self.crt_poly, self.moduli):
                product = np.poly1d(p1 * other)
                _, remainder = np.polydiv(product, self.f)
                result.append(np.poly1d(np.mod(remainder, m)))
            return CRTRq(result, self.moduli)

    def fast_basis_conversion(self, Moduli):
        """
        it just recompose and decompose again
        """
        result = []
        for Next_Mod in Moduli:
            acc_poly = np.poly1d([0])
            for poly, Cur_Mod in zip(self.crt_poly, self.moduli):
                acc_poly = poly * (Cur_Mod/self.Q)
                acc_poly = np.polymul(acc_poly, self.Q/Cur_Mod)
                acc_poly = np.poly1d(np.mod(acc_poly, Next_Mod))
            result.append(acc_poly)
        return CRTRq(result, Moduli)
                
    def fftmul(self, other):
        """
        Multiply two CRTRq objects using fft.

        Args:
            other: Another CRTRq object.

        Returns:
            A new CRTRq object representing the product.       
        """

        if not isinstance(other, CRTRq):
            raise TypeError("Multiplication is supported only between two CRTRq objects.")
        elif other.moduli != self.moduli:
            raise TypeError("Moduli chain should be same")
        
        result = []
        for p1, p2, m in zip(self.crt_poly, other.crt_poly, self.moduli):
            product = fft_multiply(np.poly1d(p1), np.poly1d(p2))
            _, remainder = poly_divide(product, self.f)
            result.append(np.poly1d(np.mod(remainder, m)))
        return CRTRq(result, self.moduli)
    


def poly_divide(dividend, divisor):
    """
    Divide two polynomials using iterative synthetic division.

    Args:
        dividend (np.poly1d): The dividend polynomial (numerator).
        divisor  (np.poly1d): The divisor polynomial (denominator).

    Returns:
        tuple: (quotient, remainder) as np.poly1d objects.
    """
    # Ensure divisor is not zero
    if np.allclose(divisor.coef, 0):
        raise ValueError("Divisor polynomial cannot be zero.")
    
    # Convert coefficients to float to handle division properly
    dividend_coefs = dividend.coef.astype(float)
    divisor_coefs = divisor.coef.astype(float)

    # Degree of the polynomials
    deg_dividend = len(dividend_coefs) - 1
    deg_divisor = len(divisor_coefs) - 1

    # Check if division is even possible
    if deg_dividend < deg_divisor:
        return np.poly1d([0]), dividend  # Quotient is zero, remainder is the dividend

    # Initialize the quotient coefficients
    quotient_coefs = np.zeros(deg_dividend - deg_divisor + 1)

    # Perform synthetic division
    for i in range(len(quotient_coefs)):
        # Compute the leading term of the quotient
        quotient_coefs[i] = dividend_coefs[i] / divisor_coefs[0]

        # Subtract the scaled divisor from the dividend
        for j in range(len(divisor_coefs)):
            dividend_coefs[i + j] -= quotient_coefs[i] * divisor_coefs[j]

    # The remainder is what's left in the dividend
    remainder_coefs = dividend_coefs[len(quotient_coefs):]

    # Convert quotient and remainder back to np.poly1d
    quotient_poly = np.poly1d(quotient_coefs)
    remainder_poly = np.poly1d(remainder_coefs)

    return quotient_poly, remainder_poly

def fft_multiply(poly1, poly2):
    """
    Multiply two polynomials using FFT.

    Args:
        poly1 (np.poly1d): First polynomial.
        poly2 (np.poly1d): Second polynomial.

    Returns:
        np.poly1d: Resultant polynomial after multiplication.
    """
    # Length of the resulting polynomial's coefficients
    result_length = len(poly1.coef) + len(poly2.coef) - 1

    # Find the next power of 2 for FFT (padding length)
    fft_length = 2**int(np.ceil(np.log2(result_length)))

    # Pad the coefficients of both polynomials to fft_length
    padded_poly1 = np.pad(poly1.coef, (0, fft_length - len(poly1.coef)))
    padded_poly2 = np.pad(poly2.coef, (0, fft_length - len(poly2.coef)))

    # Perform FFT on both polynomials
    fft_poly1 = np.fft.fft(padded_poly1)
    fft_poly2 = np.fft.fft(padded_poly2)

    # Multiply in the frequency domain
    fft_result = fft_poly1 * fft_poly2

    # Inverse FFT to get the result in the time domain
    result_coefficients = np.fft.ifft(fft_result).real.round().astype(int)[:result_length]

    # Create the resulting polynomial
    return np.poly1d(result_coefficients)