import numpy as np


def crange(coeffs, q):
    """
        For Modulus range [-q/2, q/2]
    """
    coeffs = np.where((coeffs >= 0) & (coeffs <= q//2),
                      coeffs,
                      coeffs - q)

    return coeffs