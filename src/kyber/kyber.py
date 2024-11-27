import math
from hashlib import shake_128, shake_256, sha3_512, sha3_256
import os
import time
import tracemalloc

"""
The parameters (ML-KEM-512, ML-KEM-768, and ML-KEM-1024) defined in the FIPS 203 document
https://nvlpubs.nist.gov/nistpubs/FIPS/NIST.FIPS.203.pdf
"""
DEFAULT_PARAMETERS = {
    "ML512": {"k": 2, "eta_1": 3, "eta_2": 2, "du": 10, "dv": 4},
    "ML768": {"k": 3, "eta_1": 2, "eta_2": 2, "du": 10, "dv": 4},
    "ML1024": {"k": 4, "eta_1": 2, "eta_2": 2, "du": 11, "dv": 5},
}
ML_KEM_512 = DEFAULT_PARAMETERS["ML512"]
ML_KEM_768 = DEFAULT_PARAMETERS["ML768"]
ML_KEM_1024 = DEFAULT_PARAMETERS["ML1024"]


class ML_KEM:
    def __init__(self, parameters):
        self.q = 3329
        self.n = 256
        self.zeta = 17  # Denotes the integer 17, which is a primitive n-th root of unity modulo q
        self.k = parameters["k"] # dimensions of the matrix A and the dimensions of vectors s, e, e1 and y
        self.eta_1 = parameters["eta_1"] # specify the distribution for generating the vectors s, e and y
        self.eta_2 = parameters["eta_2"] # specify the distribution for generating the vectors e1 and e2
        self.du = parameters["du"] # serve as parameters and inputs for the functions compress, decompress, byte_encode, and byte_decode in k_pke_encrypt and k_pke_decrypt
        self.dv = parameters["dv"] # serve as parameters and inputs for the functions compress, decompress, byte_encode, and byte_decode in k_pke_encrypt and k_pke_decrypt
        self.random_bytes = os.urandom
        self.zetas = [pow(self.zeta, self._bit_rev_7(i), self.q) for i in range(128)]
        self.zetas_ntt = [pow(self.zeta, 2 * self._bit_rev_7(i) + 1, self.q) for i in range(128)]

    def _poly_add(self, lst, other):
        """
        Addition of two polynomials in (Z_q)^256
        """
        assert len(lst) == len(other)
        return [(x + y) % self.q for x, y in zip(lst, other)]

    def _poly_sub(self, lst, other):
        """
        Subtraction of two polynomials in (Z_q)^256
        """
        assert len(lst) == len(other)
        return [(x - y) % self.q for x, y in zip(lst, other)]

    def _matrix_mul(self, matrix, other, error_matrix=None, transpose=False):
        """
        Implementation of 2.12 and 2.13 in FIPS 203 (page 10)
        With addition of error vector if needed

        :param matrix: matrix of size k x k
        :param other: vector of size k
        :param error_matrix: error vector of size k
        :param transpose: transpose value for the matrix
        :return: dot product of matrix with vector
        """
        assert len(matrix[0]) == len(other)
        if error_matrix is None:
            error_matrix = [[0 for _ in range(self.n)] for _ in range(self.k)]
        tmp = []
        res = [[0 for _ in range(self.n)] for _ in range(self.k)]
        return_matrix = []
        for i in range(self.k):
            for j in range(self.k):
                if transpose:
                    tmp.append(self._multiply_ntts(matrix[j][i], other[j]))
                else:
                    tmp.append(self._multiply_ntts(matrix[i][j], other[j]))
                tmp[-1] = self._poly_add(tmp[-1], error_matrix[j])
                res[i] = self._poly_add(res[i], tmp[-1])
            return_matrix.append(res[i])
        return return_matrix

    def _bits_to_bytes(self, b):
        """
        Algorithm 3 in FIPS 203 (page 20)
        Converts a bit array (of a length that is a multiple of eight) into an array of bytes
        """
        assert len(b) % 8 == 0
        B = [0 for _ in range(len(b) // 8)]
        for i in range(len(b)):
            B[i // 8] += b[i] * 2 ** (i % 8)
        return bytes(B)

    def _bytes_to_bits(self, B):
        """
        Algorithm 4 in FIPS 203 (page 20)
        Performs the inverse of BitsToBytes, converting a byte array into a bit array
        """
        b = []
        C = [i for i in B]
        for i in range(len(B)):
            for j in range(8):
                b.append(C[i] % 2)
                C[i] //= 2
        return b

    def _byte_encode(self, F, d=12):
        """
        Encodes an array of d-bit integers into a byte array for 1 <= d <= 12

        :param F: integer array F from (Z_m)^256, where m = 2^d if d < 12 else q
        :return: byte array B from B^(32 * d)
        """
        assert len(F) == 256
        b = []
        for i in range(self.n):
            a = F[i]
            for j in range(d):
                b.append(a % 2)
                a = (a - b[-1]) // 2
        return self._bits_to_bytes(b)

    def _byte_decode(self, B, d=12):
        """
        Decodes a byte array into an array of d-bit integers for 1 <= d <= 12

        :param B: byte array B from B^(32 * d)
        :return: F: integer array F from (Z_m)^256, where m = 2^d if d < 12 else q
        """
        assert len(B) == (32 * d)
        m = 1 << d if d < 12 else self.q
        b = self._bytes_to_bits(B)
        F = []
        for i in range(self.n):
            res = 0
            for j in range(d):
                res = (res % m) + (b[i * d + j] * 2 ** j) % m
            F.append(res)
        return F

    def _compress_elem(self, x, d):
        """
        Implementation of 4.7 in FIPS 203 (page 21)
        Compress_d(Decompress_d(y)) = y for all y from Z_(2^d) and all d < 12
        """
        m = 2 ** d
        if d < 12:
            y = (m / self.q) * x % m
            return self._round_int(y)
        else:
            return self.q

    def _decompress_elem(self, x, d):
        """
        Implementation of 4.8 in FIPS 203 (page 21)
        Compress_d(Decompress_d(y)) = y for all y from Z_(2^d) and all d < 12
        """
        m = 2 ** d
        if d < 12:
            y = (self.q / m) * x
            return self._round_int(y)
        else:
            return self.q

    def _compress(self, F, d=12):
        """
        Helper function to compress a polynomial F from (Z_q)^256
        """
        return [self._compress_elem(c, d) for c in F]

    def _decompress(self, F, d=12):
        """
        Helper function to decompress a polynomial F from (Z_q)^256
        """
        return [self._decompress_elem(c, d) for c in F]

    def _round_int(self, x):
        """
        Helper function to round a float to the nearest integer
        """
        assert x >= 0
        if ((x * 10) % 10) < 5:
            return math.floor(x)
        return math.ceil(x)

    def _bit_rev_7(self, r):
        """
        Bit reversal of an unsigned 7-bit integer
        """
        assert isinstance(r, int)
        assert len(bin(r)[2:]) <= 7  # bin(r)[2:] because else it would add '0b' to the length
        bin_i = bin(r & (2 ** 7 - 1))[2:].zfill(7)
        return int(bin_i[::-1], 2)

    def _xof(self, bytes32, i, j):
        """
        eXtendable-Output Function (XOF) 4.9 in FIPS 203 (page 19)
        """
        input_bytes = bytes32 + i + j
        if len(input_bytes) != 34:
            raise ValueError("Input bytes should be one 32 byte array and 2 single bytes.")
        return shake_128(input_bytes).digest(840)

    def _prf(self, eta, s, b):
        """
        Pseudorandom function 4.3 in FIPS 203 (page 18)
        """
        input_bytes = s + b
        if len(input_bytes) != 33:
            raise ValueError("Input bytes should be one 32 byte array and one single byte.")
        return shake_256(input_bytes).digest(64 * eta)

    def _H(self, s):
        """
        𝔹^* -> 𝔹^32
        Hash function 4.4 in FIPS 203 (page 18)
        """
        return sha3_256(s).digest()

    def _J(self, s):
        """
        𝔹^* -> 𝔹^32
        Hash function 4.4 in FIPS 203 (page 18)
        """
        return shake_256(s).digest(32)

    def _G(self, c):
        """
        B^* -> B^32 x B^32
        Hash function 4.5 in FIPS 203 (page 19)
        """
        g = sha3_512(c).digest()
        return g[:32], g[32:]

    def _sample_ntt(self, B):
        """
        Takes a 32-byte seed and two indices as input and outputs a pseudorandom element of T_q
        Algorithm 7 in FIPS 203 (page 23)

        :return: the coefficients of the NTT of a polynomial from (Z_q)^256
        """
        j = 0
        i = 0
        a_hat = []
        while j < self.n:
            d1 = B[i] + self.n * (B[i + 1] % 16)
            d2 = (B[i + 1] // 16) + 16 * B[i + 2]
            if d1 < self.q:
                a_hat.append(d1)
                j += 1
            if d2 < self.q and j < self.n:
                a_hat.append(d2)
                j += 1
            i += 3
        return a_hat

    def _sample_poly_cbd(self, eta, B):
        """
        Takes a seed as input and outputs a pseudorandom sample from the distribution D_eta(R_q)
        Algorithm 8 in FIPS 203 (page 23)

        :param B: byte array from B^(64 * eta)
        :return: the coefficients of the sampled polynomial from (Z_q)^256
        """
        assert 64 * eta == len(B)
        b = self._bytes_to_bits(B)
        f = []
        for i in range(self.n):
            x = 0
            y = 0
            for j in range(eta):
                x += b[2 * i * eta + j]
                y += b[2 * i * eta + eta + j]
            f.append((x - y) % self.q)
        return f

    def _to_ntt(self, f):
        """
        Computes the NTT representation of the given polynomial from R_q
        Algorithm 9 in FIPS 203 (page 26)

        :param f: given polynomial from R_q
        :return: NTT representation of the given polynomial
        """
        f_hat = f.copy()
        i = 1
        length = 128
        while length >= 2:
            start = 0
            while start < self.n:
                zeta = self.zetas[i]
                i += 1
                for j in range(start, start + length):
                    t = (zeta * f_hat[j + length]) % self.q
                    f_hat[j + length] = (f_hat[j] - t) % self.q
                    f_hat[j] = (f_hat[j] + t) % self.q
                start += (2 * length)
            length >>= 1
        return f_hat

    def _from_ntt(self, f_hat):
        """
        Computes the polynomial from R_q that corresponds to the given NTT representation from T_q
        Algorithm 10 in FIPS 203 (page 26)

        :param f_hat: the coefficients of input NTT representation
        :return: the coefficients of the inverse NTT of the input
        """
        f = f_hat.copy()
        i = 127
        length = 2
        while length <= 128:
            start = 0
            while start < self.n:
                zeta = self.zetas[i]
                i -= 1
                for j in range(start, start + length):
                    t = f[j]
                    f[j] = (f[j + length] + t) % self.q
                    f[j + length] = (zeta * (f[j + length] - t)) % self.q
                start += (2 * length)
            length *= 2
        for i in range(len(f)):
            f[i] = (f[i] * 3303) % self.q
        return f

    def _multiply_ntts(self, f_hat, g_hat):
        """
        Computes the product (in the ring T_q) of two NTT representations
        Algorithm 11 in FIPS 203 (page 27)

        :param f_hat: coefficients of NTT representations from (Z_q)^256
        :param g_hat: coefficients of NTT representations from (Z_q)^256
        :return: coefficients of the product of the inputs
        """
        assert len(f_hat) == len(g_hat)
        h_hat = []
        for i in range(128):
            c0, c1 = self._base_case_multiply(f_hat[2 * i],
                                              f_hat[2 * i + 1],
                                              g_hat[2 * i],
                                              g_hat[2 * i + 1],
                                              self.zetas_ntt[i])
            h_hat.append(c0)
            h_hat.append(c1)
        return h_hat

    def _base_case_multiply(self, a0, a1, b0, b1, gamma):
        """
        Computes the product of two degree-one polynomials with respect to a quadratic modulus
        Algorithm 12 in FIPS 203 (page 27)

        :return: coefficients of the product of the two polynomials
        """
        c0 = (a0 * b0 + a1 * b1 * gamma) % self.q
        c1 = (a0 * b1 + a1 * b0) % self.q
        return c0, c1

    def _generate_matrix_from_seed(self, rho):
        """
        Helper function which generates a matrix of size k x k from a seed rho with all elements from (Z_q)^256
        """
        A_hat = [[0 for _ in range(self.k)] for _ in range(self.k)]
        for i in range(self.k):
            for j in range(self.k):
                xof_bytes = self._xof(rho, bytes([j]), bytes([i]))
                # noinspection PyTypeChecker
                A_hat[i][j] = self._sample_ntt(xof_bytes)
        return A_hat

    def _select_bytes(self, a, b, cond):
        """
        Helper function to select between the bytes a or b depending on whether cond is False or True
        """
        assert len(a) == len(b)
        out = [0 for _ in range(len(a))]
        cw = -cond % 256
        for i in range(len(a)):
            out[i] = a[i] ^ (cw & (a[i] ^ b[i]))
        return bytes(out)

    def _k_pke_keygen(self, d):
        """
        Uses randomness to generate an encryption key and a corresponding decryption key
        Algorithm 13 in FIPS 203 (page 29)

        :param d: randomness byte array from B^32
        :return: encryption and decryption key as byte arrays ek_pke from B^(384 * k + 32) and dk_pke from B^(384 * k)
        """
        rho, sigma = self._G(d + bytes([self.k]))
        N = 0
        self.A_hat = self._generate_matrix_from_seed(rho)
        A_hat = self.A_hat
        s = []
        for i in range(self.k):
            s.append(self._sample_poly_cbd(self.eta_1, self._prf(self.eta_1, sigma, bytes([N]))))
            N += 1
        e = []
        for i in range(self.k):
            e.append(self._sample_poly_cbd(self.eta_1, self._prf(self.eta_1, sigma, bytes([N]))))
            N += 1
        s_hat = []
        e_hat = []
        for i in range(self.k):
            s_hat.append(self._to_ntt(s[i]))
            e_hat.append(self._to_ntt(e[i]))
        t_hat = self._matrix_mul(A_hat, s_hat, e_hat)
        ek_pke = b''
        dk_pke = b''
        for i in range(self.k):
            ek_pke += self._byte_encode(t_hat[i])
            dk_pke += self._byte_encode(s_hat[i])

        return ek_pke + rho, dk_pke

    def _k_pke_encrypt(self, ek_pke, m, r):
        """
        Uses the encryption key to encrypt a plaintext message using the randomness r
        Algorithm 14 in FIPS 203 (page 30)

        :param ek_pke: encryption key as byte array from B^(384 * k + 32)
        :param m: message as byte array from B^32
        :param r: randomness as byte array from B^32
        :return: ciphertext as byte array from B^(32 * (d_u * k + d_v))
        """
        if len(ek_pke) != 384 * self.k + 32:
            raise ValueError(
                f"Type check failed, ek_pke has the wrong length, expected {384 * self.k + 32} bytes and received {len(ek_pke)}")
        N = 0
        t_hat = []
        for i in range(self.k):
            t_hat.append(self._byte_decode(ek_pke[384 * i: 384 * (i + 1)]))
        if b''.join(self._byte_encode(t_hat[i]) for i in range(self.k)) != ek_pke[:-32]:
            raise ValueError("Modulus check failed, t_hat does not encode correctly")
        # rho = ek_pke[-32:]
        A_hat = self.A_hat
        y = []
        for i in range(self.k):
            y.append(self._sample_poly_cbd(self.eta_1, self._prf(self.eta_1, r, N.to_bytes(1, 'big'))))
            N += 1
        e1 = []
        for i in range(self.k):
            e1.append(self._sample_poly_cbd(self.eta_2, self._prf(self.eta_2, r, N.to_bytes(1, 'big'))))
            N += 1
        e2 = self._sample_poly_cbd(self.eta_2, self._prf(self.eta_2, r, N.to_bytes(1, 'big')))
        y_hat = [self._to_ntt(poly) for poly in y]
        u_hat = self._matrix_mul(A_hat, y_hat, transpose=True)
        u = []
        for i in range(self.k):
            u.append(self._poly_add(self._from_ntt(u_hat[i]), e1[i]))
        mu = self._decompress(self._byte_decode(m, 1), 1)
        v = [0 for _ in range(self.n)]
        for i in range(self.k):
            tmp = self._multiply_ntts(t_hat[i], y_hat[i])
            v = self._poly_add(v, tmp)
        v = self._from_ntt(v)
        v = self._poly_add(v, e2)
        v = self._poly_add(v, mu)
        c1 = b''
        for i in range(self.k):
            c1 += self._byte_encode(self._compress(u[i], self.du), self.du)
        c2 = self._byte_encode(self._compress(v, self.dv), self.dv)
        return c1 + c2

    def _k_pke_decrypt(self, dk_pke, c):
        """
        Uses the decryption key to decrypt a ciphertext
        Algorithm 15 in FIPS 203 (page 31)

        :param dk_pke: decryption key as byte array from B^(384 * k)
        :param c: ciphertext as byte array from B^(32 * (d_u * k + d_v))
        :return: message as byte array from B^32
        """
        n = self.k * self.du * 32
        c1, c2 = c[:n], c[n:]
        u = []
        for i in range(self.k):
            u.append(
                self._decompress(self._byte_decode(c1[i * self.du * 32: (i + 1) * self.du * 32], self.du), self.du))
        v = self._decompress(self._byte_decode(c2, self.dv), self.dv)
        s_hat = []
        u_hat = []
        for i in range(self.k):
            s_hat.append(self._byte_decode(dk_pke[i * 384: (i + 1) * 384]))
            u_hat.append(self._to_ntt(u[i]))
        s_mul_u_ntt = [0 for _ in range(self.n)]
        for i in range(self.k):
            tmp = self._multiply_ntts(s_hat[i], u_hat[i])
            s_mul_u_ntt = self._poly_add(s_mul_u_ntt, tmp)
        w = self._poly_sub(v, self._from_ntt(s_mul_u_ntt))
        m = self._byte_encode(self._compress(w, 1), 1)
        return m

    def _keygen_internal(self, d, z):
        """
        Uses randomness to generate an encapsulation key and a corresponding decapsulation key
        Algorithm 16 in FIPS 203 (page 32)

        :param d: randomness byte array from B^32
        :param z: randomness byte array from B^32
        :return: encapsulation and decapsulation key as byte arrays ek from B^(384 * k + 32) and dk from B^(768 * k + 96)
        """
        ek_pke, dk_pke = self._k_pke_keygen(d)
        ek = ek_pke
        dk = dk_pke + ek + self._H(ek) + z
        return ek, dk

    def _encaps_internal(self, ek, m):
        """
        Uses the encapsulation key and randomness to generate a key and an associated ciphertext
        Algorithm 17 in FIPS 203 (page 33)

        :param ek: encapsulation key from B^(384 * k + 32)
        :param m: randomness from B^32
        :return: shared secret key K from B^32 and ciphertext c from B^(32 * (d_u * k + d_v))
        """
        K, r = self._G(m + self._H(ek))
        c = self._k_pke_encrypt(ek, m, r)
        return K, c

    def _decaps_internal(self, dk, c):
        """
        Uses the decapsulation key to produce a shared secret key from a ciphertext
        Algorithm 18 in FIPS 203 (page 34)

        :param dk: decapsulation key from B^(768 * k + 96)
        :param c: ciphertext from B^(32 * (d_u * k + d_v))
        :return: shared secret key K from B^32
        """
        if len(c) != 32 * (self.du * self.k + self.dv):
            raise ValueError(
                f"ciphertext type check failed. Expected {32 * (self.du * self.k + self.dv)} bytes and obtained {len(c)}")
        if len(dk) != 768 * self.k + 96:
            raise ValueError(
                f"decapsulation type check failed. Expected {768 * self.k + 96} bytes and obtained {len(dk)}")
        dk_pke = dk[: 384 * self.k]
        ek_pke = dk[384 * self.k: 768 * self.k + 32]
        h = dk[768 * self.k + 32: 768 * self.k + 64]
        z = dk[768 * self.k + 64:]
        if self._H(ek_pke) != h:
            raise ValueError("hash check failed")
        m_prime = self._k_pke_decrypt(dk_pke, c)
        K_prime, r_prime = self._G(m_prime + h)
        K_bar = self._J(z + c)
        c_prime = self._k_pke_encrypt(ek_pke, m_prime, r_prime)
        return self._select_bytes(K_bar, K_prime, c == c_prime)

    def keygen(self):
        """
        Generates an encapsulation key and a corresponding decapsulation key
        Algorithm 19 in FIPS 203 (page 35)

        :return: encapsulation and decapsulation key as byte arrays ek from B^(384 * k + 32) and dk from B^(768 * k + 96)
        """
        d = self.random_bytes(32)
        z = self.random_bytes(32)
        ek, dk = self._keygen_internal(d, z)
        return ek, dk

    def encaps(self, ek):
        """
        Uses the encapsulation key to generate a shared secret key and an associated ciphertext
        Algorithm 20 in FIPS 203 (page 37)

        :param ek: encapsulation key from B^(384 * k + 32)
        :return: shared secret key K from B^32 and ciphertext c from B^(32 * (d_u * k + d_v))
        """
        if len(ek) != self.k * 384 + 32:
            raise ValueError("")
        m = self.random_bytes(32)
        K, c = self._encaps_internal(ek, m)
        return K, c

    def decaps(self, dk, c):
        """
        Uses the decapsulation key to produce a shared secret key from a ciphertext
        Algorithm 21 in FIPS 203 (page 38)

        :param dk: decapsulation key from B^(768 * k + 96)
        :param c: ciphertext from B^(32 * (d_u * k + d_v))
        :return: shared secret key K from B^32
        """
        K_prime = self._decaps_internal(dk, c)
        return K_prime

def test_runtime(runs, parameters):
    ml = ML_KEM(parameters)
    total_time = []
    keygen_time = []
    encaps_time = []
    decaps_time = []
    for run in range(runs):
        t0 = time.time()
        encaps_key, decaps_key = ml.keygen()
        t1 = time.time()
        K, c = ml.encaps(encaps_key)
        t2 = time.time()
        K_prime = ml.decaps(decaps_key, c)
        t3 = time.time()
        assert K == K_prime
        total_time.append(t3 - t0)
        keygen_time.append(t1 - t0)
        encaps_time.append(t2 - t1)
        decaps_time.append(t3 - t2)
    print("NUMBER OF ITERATIONS: ", runs)
    print("\n######## TOTAL TIME ########")
    print("average time: ", sum(total_time) / len(total_time))
    print("total time: ", sum(total_time))
    print("number of iterations per second: ", 1 / (sum(total_time) / len(total_time)))
    #print(total_time)

    print("\n######## KEY GENERATION TIME ########")
    print("average time: ", sum(keygen_time) / len(keygen_time))
    print("total time: ", sum(keygen_time))
    print("number of iterations per second: ", 1 / (sum(keygen_time) / len(keygen_time)))
    #print(keygen_time)

    print("\n######## ENCAPSULATION TIME ########")
    print("average time: ", sum(encaps_time) / len(encaps_time))
    print("total time: ", sum(encaps_time))
    print("number of iterations per second: ", 1 / (sum(encaps_time) / len(encaps_time)))
    #print(encaps_time)

    print("\n######## DECAPSULATION TIME ########")
    print("average time: ", sum(decaps_time) / len(decaps_time))
    print("total time: ", sum(decaps_time))
    print("number of iterations per second: ", 1 / (sum(decaps_time) / len(decaps_time)))
    #print(decaps_time)

def test_memory_use(parameters):
    ml = ML_KEM(parameters)

    tracemalloc.start()
    encaps_key, decaps_key = ml.keygen()
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()

    tracemalloc.start()
    K, c = ml.encaps(encaps_key)
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()

    tracemalloc.start()
    K_prime = ml.decaps(decaps_key, c)
    assert K == K_prime
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()

if __name__ == '__main__':
    test_runtime(1000, ML_KEM_1024)
    #test_memory_use(ML_KEM_1024)
