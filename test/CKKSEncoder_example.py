from include.CKKSEncoder import CKKSEncoder
import numpy as np

if __name__ == '__main__':
    n = 16  # power of 2
    scale = 64

    encoder = CKKSEncoder(n, scale)
    
    # z = np.array([3 + 4j, 2 - 1j, 1 + 1j, 2 + 2j, 4 - 3j, 6 + 9j, 3 + 2j, 1 - 8j])
    # z1 = np.array([3 + 4j, 2 - 1j, 1 + 1j, 2 + 2j])
    # z2 = np.array([4 - 3j, 6 + 9j, 3 + 2j, 1 - 8j])

    z1 = np.array([3.3, 2.3, 1.3, 2.3])
    z2 = np.array([4, 6, 3, 1])

    p1 = encoder.encode(z1)
    p2 = encoder.encode(z2)

    p = p1 + p2

    decoded_z = encoder.decode(p)
    print(decoded_z.real)