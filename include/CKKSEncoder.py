import numpy as np
from numpy.polynomial import Polynomial

class CKKSEncoder:
    def __init__(self, M: int, scale: float) -> None:
        """
        Initialize CKKS encoder
        """
        self.xi = np.exp(2 * np.pi * 1j / M)
        self.M  = M
        self.create_sigma_R_basis()
        self.scale = scale

    @staticmethod
    def vandermonde(xi: np.complex128, M: int) -> np.array:
        N = M // 2
        matrix = []

        for i in range(N):
            root = xi ** (2 * i + 1)
            row = []

            for j in range(N):
                row.append(root ** j)
            matrix.append(row)
        return matrix
    
    def sigma_inverse(self, b: np.array) -> Polynomial:
        A = CKKSEncoder.vandermonde(self.xi, self.M)
        coeffs = np.linalg.solve(A, b)
        p = Polynomial(coeffs)
        return p
    
    def sigma(self, p: Polynomial) -> np.array:
        outputs = []
        N = self.M // 2

        for i in range(N):
            root = self.xi ** (2 * i + 1)
            output = p(root)
            outputs.append(output)
        return np.array(outputs)
    
    def pi(self, z: np.array) -> np.array:
        """ Project a vector of H into C^{N/2}. """
        N = self.M // 4
        return z[:N]
    
    def pi_inverse(self, z: np.array) -> np.array:
        """ Expands a vector of C^{N/2} by expanding it with its complex conjugate. """
        z_conjugate = z[::-1]
        z_conjugate = [np.conjugate(x) for x in z_conjugate]
        return np.concatenate([z, z_conjugate])
    
    def create_sigma_R_basis(self):
        """ The basis (sigma(1), sigma(X), ..., sigma(X^{N-1})) """
        self.sigma_R_basis = np.array(self.vandermonde(self.xi, self.M)).T
    
    def compute_basis_coordinates(self, z):
        """ Computes the coordinates of a vector with respect to the orthogonal lattice basis. """
        output = np.array([np.real(np.vdot(z,b) / np.vdot(b,b)) for b in self.sigma_R_basis])
        return output
    
    def sigma_R_discretization(self, z):
        """Projects a vector on the lattice using coordinate wise random rounding."""
        coordinates = self.compute_basis_coordinates(z)

        rounded_coordinates = coordinate_wise_random_rounding(coordinates)
        y = np.matmul(self.sigma_R_basis.T, rounded_coordinates)
        return y
    
    def encode(self, z: np.array) -> Polynomial:
        """
            Encodes a vector by expanding it first to H,
            scale it, project it on the lattice of sigma(R), and performs
            sigma inverse.
        """
        pi_z = self.pi_inverse(z)
        scaled_pi_z = self.scale * pi_z
        rounded_scale_pi_zi = self.sigma_R_discretization(scaled_pi_z)
        p = self.sigma_inverse(rounded_scale_pi_zi)

        coef = np.round(np.real(p.coef)).astype(int)
        p = Polynomial(coef)
        return p

    def decode(self, p: Polynomial) -> np.array:
        """
            Decodes a polynomial by removing the scale, 
            evaluating on the roots, and project it on C^(N/2)
        """
        rescaled_p = p / self.scale
        z = self.sigma(rescaled_p)
        pi_z = self.pi(z)
        return pi_z
        

def coordinate_wise_random_rounding(coordinates):
    """Rounds coordinates randonmly."""
    r = coordinates - np.floor(coordinates) # round coordinates
    f = np.array([np.random.choice([c, c-1], 1, p=[1-c, c]) for c in r]).reshape(-1)

    rounded_coordinates = coordinates - r
    rounded_coordinates = [int(coeff) for coeff in rounded_coordinates]
    return rounded_coordinates