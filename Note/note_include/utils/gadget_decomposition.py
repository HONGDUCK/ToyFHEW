import numpy as np
from note_include.elem.Ring import Ring
# from note_include.utils.types import RGSWctxt, RLWEctxt, RLWEpctxt

def gadget_decomposition(poly : Ring, B : int, d : int) -> list[Ring]:
    coeffs     = np.array(poly.coeffs)
    decomposed = []

    for _ in range(d):  # Decompose into d components
        remainder = coeffs % B  # Get the lowest base-B component
        coeffs = coeffs   // B  # Reduce the coefficients
        decomposed.append(Ring(poly.n, poly.q, remainder))  # Store decomposition step

    return decomposed  # Return vector of polynomials

def gadget_composition(decomposed: list[Ring], B: int, modulus: int) -> Ring:
    coeffs = np.zeros_like(decomposed[0].coeffs)  # Initialize with zeros

    for i, d in enumerate(decomposed):
        coeffs += (B**i) * np.array(d.coeffs)  # Sum up each term
    
    coeffs = [coef % modulus for coef in coeffs]
    return Ring(decomposed[0].n, decomposed[0].q, coeffs)

def format_ring_list(ring_list):
    return "\n".join(repr(ring) for ring in ring_list)

def gadget_decomposition_int(v:int, B : int, d : int) -> list[int]:
    decomposed = []
    for _ in range(d):  # Decompose into d components
        remainder = v % B  # Get the lowest base-B component
        v = v // B  # Reduce the coefficients
        decomposed.append(remainder)  # Store decomposition step

    return decomposed  # Return vector of polynomials