import numpy as np
from note_include.Ring import Ring

def discrete_gaussian(n, q, mean=0., std=1.):
    coeffs = np.round(std * np.random.randn(n)) % q
    return Ring(n, q, coeffs)

def discrete_uniform(n, q, min=0., max=None):
    if max is None:
        max = q
    coeffs = np.random.randint(min, max, size=n)
    return Ring(n, q, coeffs)