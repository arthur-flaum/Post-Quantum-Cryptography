import math
import os
from hashlib import shake_128, shake_256
import time

DEFAULT_PARAMETERS = {
    "ML_DSA_44": {
        "d": 13,  # number of bits dropped from t
        "tau": 39,  # number of ±1 in c
        "gamma_1": 131072,  # coefficient range of y: 2^17
        "gamma_2": 95232,  # low order rounding range: (q-1)/88
        "k": 4,  # Dimensions of A = (k, l)
        "l": 4,  # Dimensions of A = (k, l)
        "eta": 2,  # Private key range
        "omega": 80,  # Max number of ones in hint
        "c_tilde_bytes": 32,
    },
    "ML_DSA_65": {
        "d": 13,  # number of bits dropped from t
        "tau": 49,  # number of ±1 in c
        "gamma_1": 524288,  # coefficient range of y: 2^19
        "gamma_2": 261888,  # low order rounding range: (q-1)/32
        "k": 6,  # Dimensions of A = (k, l)
        "l": 5,  # Dimensions of A = (k, l)
        "eta": 4,  # Private key range
        "omega": 55,  # Max number of ones in hint
        "c_tilde_bytes": 48,
    },
    "ML_DSA_87": {
        "d": 13,  # number of bits dropped from t
        "tau": 60,  # number of ±1 in c
        "gamma_1": 524288,  # coefficient range of y: 2^19
        "gamma_2": 261888,  # low order rounding range: (q-1)/32
        "k": 8,  # Dimensions of A = (k, l)
        "l": 7,  # Dimensions of A = (k, l)
        "eta": 2,  # Private key range
        "omega": 75,  # Max number of ones in hint
        "c_tilde_bytes": 64,
    },
}
ML_DSA_44 = DEFAULT_PARAMETERS["ML_DSA_44"]
ML_DSA_65 = DEFAULT_PARAMETERS["ML_DSA_65"]
ML_DSA_87 = DEFAULT_PARAMETERS["ML_DSA_87"]





class ML_DSA:
    def __init__(self, parameters):
        self.q = 8380417 # q = 2^23 - 2^13 + 1
        self.n = 256
        self.d = parameters["d"]
        self.k = parameters["k"]
        self.l = parameters["l"]
        self.eta = parameters["eta"]
        self.tau = parameters["tau"]
        self.omega = parameters["omega"]
        self.gamma_1 = parameters["gamma_1"]
        self.gamma_2 = parameters["gamma_2"]
        self.beta = self.tau * self.eta
        self.c_tilde_bytes = parameters["c_tilde_bytes"]
        self.random_bytes = os.urandom
        self.zetas = [pow(1753, self._bit_rev_8(i), self.q) for i in range(256)]

    def _bit_rev_8(self, m):
        b = self._int_to_bits(m, 8)
        b_rev = []
        for i in range(8):
            b_rev.append(b[7 - i])
        return self._bits_to_int(b_rev, 8)


    def _poly_add(self, lst, other):
        """
        Addition of two polynomials in (Z_q)^256
        """
        assert len(lst) == len(other)
        return [(x + y) for x, y in zip(lst, other)]

    def _poly_sub(self, lst, other):
        """
        Subtraction of two polynomials in (Z_q)^256
        """
        assert len(lst) == len(other)
        return [(x - y) for x, y in zip(lst, other)]

    def _int_to_bits(self, x, alpha):
        x_prime = x
        y = []
        for i in range(alpha):
            y.append(x_prime % 2)
            x_prime >>= 1
        return y

    def _bits_to_int(self, y, alpha):
        x = 0
        for i in range(1, alpha+1):
            x = (x << 1) + y[alpha - i]
        return x

    def _int_to_bytes(self, x, alpha):
        x_prime = x
        y = []
        for i in range(alpha):
            y.append(x_prime % self.n)
            x_prime >>= 8
        return bytes(y)

    def _bits_to_bytes(self, y):
        z = [0 for _ in range(math.ceil(len(y) / 8))]
        for i in range(len(y)):
            z[i // 8] += y[i] * 2**(i % 8)
        return bytes(z)

    def _bytes_to_bits(self, z):
        z_prime = [i for i in z]
        y = []
        for i in range(len(z)):
            for j in range(8):
                y.append(z_prime[i] % 2)
                z_prime[i] >>= 1
        return y

    def _coeff_from_three_bytes(self, b0, b1, b2):
        b2_prime = b2
        if b2_prime > 127:
            b2_prime -= 128
        z = (1<<16) * b2_prime + (1<<8) * b1 + b0
        if z < self.q:
            return z
        return False

    def _coeff_from_half_byte(self, b):
        if self.eta == 2 and b < 15:
            return 2 - (b % 5)
        if self.eta == 4 and b < 9:
            return 4 - b
        return False

    def _bit_len(self, b):
        return len(bin(b)[2:])

    def _simple_bit_pack(self, w, b):
        z = []
        bit_len = self._bit_len(b)
        for i in range(self.n):
            z.extend(self._int_to_bits(w[i], bit_len))
        return self._bits_to_bytes(z)

    def _bit_pack(self, w, a, b):
        z = []
        bit_len = self._bit_len(a + b)
        for i in range(self.n):
            z.extend(self._int_to_bits(b - w[i], bit_len))
        return self._bits_to_bytes(z)

    def _simple_bit_unpack(self, v, b):
        c = self._bit_len(b)
        assert len(v) == 32 * c
        z = self._bytes_to_bits(v)
        w = []
        for i in range(self.n):
            w.append(self._bits_to_int(z[i * c : i * c + c], c))
        return w

    def _bit_unpack(self, v, a, b):
        c = self._bit_len(a + b)
        z = self._bytes_to_bits(v)
        w = []
        for i in range(self.n):
            w.append(b - self._bits_to_int(z[i * c : i * c + c], c))
        return w

    def _hint_bit_pack(self, h):
        y = [0 for _ in range(self.omega + self.k)]
        index = 0
        for i in range(self.k):
            for j in range(self.n):
                if index >= len(y):
                    raise IndexError(f"There are more hint bits than allowed with omega = {self.omega}")
                if h[i][j] != 0:
                    y[index] = j
                    index += 1
            y[self.omega + i] = index
        return bytes(y)

    def _hint_bit_unpack(self, y):
        h = [[0 for _ in range(self.n)] for _ in range(self.k)]
        index = 0
        for i in range(self.k):
            tmp = y[self.omega + i]
            if tmp < index or tmp > self.omega:
                return False
            first = index
            while index < tmp:
                if index > first:
                    if y[index - 1] >= y[index]:
                        return False
                h[i][y[index]] = 1
                index += 1
        for i in range(index, self.omega):
            if y[i] != 0:
                return False
        return h

    def _pk_encode(self, rho, t1):
        max_num = (1 << (self._bit_len(self.q - 1) - self.d)) - 1
        pk = rho
        for i in range(self.k):
            tmp = self._simple_bit_pack(t1[i], max_num)
            pk += tmp
        return pk

    def _pk_decode(self, pk):
        z_length = (self._bit_len(self.q - 1) - self.d)
        max_num = (1 << z_length) - 1
        rho, z = pk[:32], pk[32:]
        t1 = []
        for i in range(self.k):
            z_i = z[32 * i * z_length : 32 * (i+1) * z_length]
            t1.append(self._simple_bit_unpack(z_i, max_num))
        return rho, t1

    def _sk_encode(self, rho, K, tr, s1, s2, t0):
        bit_pack_range = (1 << (self.d - 1))
        sk = rho + K + tr
        for i in range(self.l):
            sk += self._bit_pack(s1[i], self.eta, self.eta)
        for i in range(self.k):
            sk += self._bit_pack(s2[i], self.eta, self.eta)
        for i in range(self.k):
            sk += self._bit_pack(t0[i], bit_pack_range - 1, bit_pack_range)
        return sk

    def _sk_decode(self, sk):
        byte_length = 32 * self._bit_len(2 * self.eta)
        bit_pack_range = (1 << (self.d - 1))
        rho, K, tr = sk[:32], sk[32:64], sk[64:128]
        y = sk[128: 128 + self.l * byte_length]

        tmp = 128 + self.l * byte_length
        z = sk[tmp : tmp + self.k * byte_length]

        tmp += self.k * byte_length
        w = sk[tmp : tmp + 32 * self.d * self.k]

        s1 = []
        s2 = []
        t0 = []
        for i in range(self.l):
            y_i = y[i * byte_length : (i + 1) * byte_length]
            s1.append(self._bit_unpack(y_i, self.eta, self.eta))
        for i in range(self.k):
            z_i = z[i * byte_length : (i + 1) * byte_length]
            s2.append(self._bit_unpack(z_i, self.eta, self.eta))
        for i in range(self.k):
            w_i = w[i * 32 * self.d : (i + 1) * 32 * self.d]
            t0.append(self._bit_unpack(w_i, bit_pack_range - 1, bit_pack_range))
        return rho, K, tr, s1, s2, t0

    def _sig_encode(self, c_tilde, z, h):
        sigma = c_tilde
        for i in range(self.l):
            sigma += self._bit_pack(z[i], self.gamma_1 - 1, self.gamma_1)
        sigma += self._hint_bit_pack(h)
        return sigma

    def _sig_decode(self, sigma):
        byte_length = 32 * (1 + self._bit_len(self.gamma_1 - 1))
        c_tilde = sigma[:self.c_tilde_bytes]
        x = sigma[self.c_tilde_bytes : self.c_tilde_bytes + self.l * byte_length]

        tmp = self.c_tilde_bytes + self.l * byte_length
        y = sigma[tmp : tmp + self.omega + self.k]

        z = []
        for i in range(self.l):
            x_i = x[i * byte_length : (i + 1) * byte_length]
            z.append(self._bit_unpack(x_i, self.gamma_1 - 1, self.gamma_1))
        h = self._hint_bit_unpack(y)
        return c_tilde, z, h

    def _w1_encode(self, w1):
        w1_tilde = b''
        for i in range(self.k):
            w1_tilde += self._simple_bit_pack(w1[i], (self.q - 1) // (2 * self.gamma_2) - 1)
        return w1_tilde

    def _H(self, str, l):
        return shake_256(str).digest(l)

    def _G(self, str, l):
        return shake_128(str).digest(l)

    def _sample_in_ball(self, rho):
        c = [0 for _ in range(self.n)]
        tmp = 8
        s = self._H(rho, 229)
        h = self._bytes_to_bits(s[:tmp])
        for i in range(self.n - self.tau, self.n): # put n - 1 instead n to the end range
            j = s[tmp : tmp+1]
            tmp += 1
            while j[0] > i:
                j = s[tmp : tmp+1]
                tmp += 1
            c[i] = c[j[0]]
            c[j[0]] = (-1)**h[i + self.tau - self.n]
        return c

    def _rej_ntt_poly(self, rho):
        j = 0
        tmp = 0
        s = self._G(rho, 894)
        a_hat = []
        while j < self.n:
            s0, s1, s2 = s[tmp : tmp+1], s[tmp+1 : tmp+2], s[tmp+2 : tmp+3]
            tmp += 3
            a_hat.append(self._coeff_from_three_bytes(s0[0], s1[0], s2[0]))
            if a_hat[-1] is False:
                a_hat.pop()
                continue
            j += 1
        return a_hat

    def _rej_bounded_poly(self, rho):
        j = 0
        tmp = 0
        z = self._H(rho, 481)
        a = [0 for _ in range(self.n)]
        while j < self.n:
            z0 = self._coeff_from_half_byte(z[tmp : tmp+1][0] % 16)
            z1 = self._coeff_from_half_byte(z[tmp : tmp+1][0] // 16)
            tmp += 1
            if z0 is not False:
                a[j] = z0
                j += 1
            if z1 is not False and j < self.n:
                a[j] = z1
                j += 1
        return a

    def _expand_A(self, rho):
        A_hat = [[] for _ in range(self.k)]
        for r in range(self.k):
            for s in range(self.l):
                rho_prime = rho + self._int_to_bytes(s, 1) + self._int_to_bytes(r, 1)
                A_hat[r].append(self._rej_ntt_poly(rho_prime))
        return A_hat

    def _expand_S(self, rho):
        s1 = []
        s2 = []
        for r in range(self.l):
            s1.append(self._rej_bounded_poly(rho + self._int_to_bytes(r, 2)))
        for r in range(self.k):
            s2.append(self._rej_bounded_poly(rho + self._int_to_bytes(r + self.l, 2)))
        return s1, s2

    def _expand_mask(self, rho, mu):
        c = 1 + self._bit_len(self.gamma_1 - 1)
        y = []
        for r in range(self.l):
            rho_prime = rho + self._int_to_bytes(mu + r, 2)
            v = self._H(rho_prime, 32 * c)
            y.append(self._bit_unpack(v, self.gamma_1 - 1, self.gamma_1))
        return y

    def _mod_plus_minus(self, m, a):
        """
        If a is a positive integer and m from Z or m from Z_a, then m mod^(+-) a denotes the unique
        element m' from Z in the range −⌈a/2⌉ < m' <= ⌊a/2⌋ such that m and m' are
        congruent modulo a.
        """
        m_prime = m % a
        if m_prime > a // 2:
            m_prime -= a
        return m_prime

    def _pow_to_round(self, r):
        r_prime = r % self.q
        mod = (1 << self.d)
        r0 = self._mod_plus_minus(r_prime, mod)
        return (r_prime - r0) // mod, r0

    def _pow_to_round_vector(self, vector):
        t = [[self._pow_to_round(vector[i][j]) for j in range(self.n)] for i in range(len(vector))]
        t1 = [[t[i][j][0] for j in range(self.n)] for i in range(len(t))]
        t0 = [[t[i][j][1] for j in range(self.n)] for i in range(len(t))]
        return t1, t0

    def _decompose(self, r):
        r_prime = r % self.q
        mod = 2 * self.gamma_2
        r0 = self._mod_plus_minus(r_prime, mod)
        r1 = 0
        if r_prime - r0 == self.q - 1:
            r0 -= 1
        else:
            r1 = (r_prime - r0) // mod
        return r1, r0

    def _high_bits(self, r):
        r1, r0 = self._decompose(r)
        return r1

    def _high_bits_vector(self, vector):
        return [[self._high_bits(vector[i][j]) for j in range(self.n)] for i in range(len(vector))]

    def _low_bits(self, r):
        r1, r0 = self._decompose(r)
        return r0

    def _low_bits_vector(self, vector):
        return [[self._low_bits(vector[i][j]) for j in range(self.n)] for i in range(len(vector))]

    def _make_hint(self, z, r):
        r1 = self._high_bits(r)
        v1 = self._high_bits(r + z)
        return r1 != v1

    def _use_hint(self, h, r):
        m = (self.q - 1) //  (2 * self.gamma_2)
        r1, r0 = self._decompose(r)
        if h == 1 and r0 > 0:
            return (r1 + 1) % m
        if h == 1 and r0 <= 0:
            return (r1 - 1) % m
        return r1

    def _use_hint_vector(self, vector, other):
        return [[self._use_hint(vector[i][j], other[i][j]) for j in range(self.n)] for i in range(len(vector))]

    def _to_ntt(self, w):
        w_hat = w.copy()
        m = 0
        length = 128
        while length >= 1:
            start = 0
            while start < self.n:
                m += 1
                z = self.zetas[m]
                for j in range(start, start + length):
                    t = (z * w_hat[j + length]) % self.q
                    w_hat[j + length] = (w_hat[j] - t) % self.q
                    w_hat[j] = (w_hat[j] + t) % self.q
                start += (length << 1)
            length >>= 1
        return w_hat

    def _from_ntt(self, w_hat) -> list[list[int] | int]:
        w = w_hat.copy()
        m = 256
        length = 1
        while length < self.n:
            start = 0
            while start < self.n:
                m -= 1
                z = -self.zetas[m]
                for j in range(start, start + length):
                    t = w[j]
                    w[j] = (w[j + length] + t) % self.q
                    w[j + length] = (t - w[j + length]) % self.q
                    w[j + length] = (z * (w[j + length])) % self.q
                start += (length << 1)
            length <<= 1
        f = 8347681     # f = 256^(-1) mod q
        for j in range(self.n):
            w[j] = (f * w[j]) % self.q
        if len(w) == self.n:
            return [self._mod_plus_minus(i, self.q) for i in w]
        return [[self._mod_plus_minus(i, self.q) for i in row] for row in w]

    def _add_ntt(self, a_hat, b_hat):
        return [(x + y) % self.q for x, y in zip(a_hat, b_hat)]

    def _mult_ntt(self, a_hat, b_hat):
        return [(x * y) % self.q for x, y in zip(a_hat, b_hat)]

    def _add_vector_ntt(self, v_hat, w_hat):
        u_hat = []
        for i in range(self.l):
            u_hat.append(self._add_ntt(v_hat[i], w_hat[i]))
        return u_hat

    def _scalar_vector_ntt(self, c_hat, v_hat):
        w_hat = []
        for i in range(len(v_hat)):
            w_hat.append(self._mult_ntt(c_hat, v_hat[i]))
        return w_hat

    def _matrix_vector_ntt(self, M_hat, v_hat):
        w_hat = [[0 for _ in range(self.n)] for _ in range(self.k)]
        for i in range(self.k):
            for j in range(self.l):
                w_hat[i] = self._add_ntt(w_hat[i], self._mult_ntt(M_hat[i][j], v_hat[j]))
        return w_hat

    def _infinity_norm_polynomial(self, polynomial):
        return max(abs(x) for x in polynomial)

    def _infinity_norm_vector(self, vector):
        return max(self._infinity_norm_polynomial(x) for x in vector)

    def _count_1_in_hint(self, h):
        return sum(h[i][j]for j in range(self.n) for i in range(self.k))

    def _ml_dsa_keygen_internal(self, xi):
        bytes_output = self._H(xi + self._int_to_bytes(self.k, 1) + self._int_to_bytes(self.l, 1), 128)
        rho, rho_prime, K = bytes_output[:32], bytes_output[32:96], bytes_output[96: 128]
        A_hat = self._expand_A(rho)
        s1, s2 = self._expand_S(rho_prime)
        s1_hat = [self._to_ntt(s1[i]) for i in range(self.l)]
        As1 = [self._from_ntt(self._matrix_vector_ntt(A_hat, s1_hat)[i]) for i in range(self.k)]
        t = [self._poly_add(As1[i], s2[i]) for i in range(self.k)]
        t1, t0 = self._pow_to_round_vector(t)
        pk = self._pk_encode(rho, t1)
        tr = self._H(pk, 64)
        sk = self._sk_encode(rho, K, tr, s1, s2, t0)
        return pk, sk

    def _ml_dsa_sign_internal(self, sk, M_prime, rnd):
        rho, K, tr, s1, s2, t0 = self._sk_decode(sk)
        s1_hat = [self._to_ntt(s1[i]) for i in range(self.l)]
        s2_hat = [self._to_ntt(s2[i]) for i in range(self.k)]
        t0_hat = [self._to_ntt(t0[i]) for i in range(self.k)]
        A_hat = self._expand_A(rho)
        mu = self._H(tr + M_prime, 64)
        rho_prime = self._H(K + rnd + mu, 64)
        kappa = 0
        c_tilde = b''
        z, h = False, False
        while z is False or h is False:
            y = self._expand_mask(rho_prime, kappa)
            kappa += self.l
            y_hat = [self._to_ntt(y[i]) for i in range(self.l)]
            w = [self._from_ntt(self._matrix_vector_ntt(A_hat, y_hat)[i]) for i in range(self.k)]
            w1 = self._high_bits_vector(w)

            c_tilde = self._H(mu + self._w1_encode(w1), self.c_tilde_bytes)
            c = self._sample_in_ball(c_tilde)
            c_hat = self._to_ntt(c)
            cs1 = [self._from_ntt(self._mult_ntt(c_hat, s1_hat[i])) for i in range(self.l)]
            z = [self._poly_add(y[i], cs1[i]) for i in range(self.l)]
            if self._infinity_norm_vector(z) >= (self.gamma_1 - self.beta):
                z, h = False, False
                continue

            cs2 = [self._from_ntt(self._mult_ntt(c_hat, s2_hat[i])) for i in range(self.k)]
            low_bits_vector = [self._poly_sub(w[i], cs2[i]) for i in range(self.k)]
            r0 = self._low_bits_vector(low_bits_vector)
            if self._infinity_norm_vector(r0) >= (self.gamma_2 - self.beta):
                z, h = False, False
                continue

            ct0 = [self._from_ntt(self._mult_ntt(c_hat, t0_hat[i])) for i in range(self.k)]
            h = [[0 for _ in range(self.n)] for _ in range(self.k)]
            for i in range(self.k):
                for j in range(self.n):
                    h[i][j] = self._make_hint(-ct0[i][j], w[i][j] - cs2[i][j] + ct0[i][j])
            if self._infinity_norm_vector(ct0) >= self.gamma_2 or self._count_1_in_hint(h) > self.omega:
                z, h = False, False

        sigma = self._sig_encode(c_tilde, z, h)
        return sigma

    def _ml_dsa_verify_internal(self, pk, M_prime, sigma):
        rho, t1 = self._pk_decode(pk)
        c_tilde, z, h = self._sig_decode(sigma)
        if h is False:
            return False
        A_hat = self._expand_A(rho)
        tr = self._H(pk, 64)
        mu = self._H(tr + M_prime, 64)
        c = self._sample_in_ball(c_tilde)
        z_hat = [self._to_ntt(z[i]) for i in range(self.l)]
        c_hat = self._to_ntt(c)
        d = (1 << self.d)
        scalar_polynomial = [d for _ in range(self.n)]
        t1_hat = [self._to_ntt(self._scalar_vector_ntt(scalar_polynomial, t1)[i]) for i in range(self.k)]
        ct1_hat = self._scalar_vector_ntt(c_hat, t1_hat)
        Az_hat = self._matrix_vector_ntt(A_hat, z_hat)
        w_prime_approx = [self._from_ntt(self._poly_sub(Az_hat[i], ct1_hat[i])) for i in range(self.k)]
        w1_prime = self._use_hint_vector(h, w_prime_approx)
        c_tilde_prime = self._H(mu + self._w1_encode(w1_prime), self.c_tilde_bytes)
        return self._infinity_norm_vector(z) < (self.gamma_1 - self.beta) and c_tilde == c_tilde_prime

    def keygen(self):
        """
        Generates a public-private key pair

        :return:
        """
        e = self.random_bytes(32)
        return self._ml_dsa_keygen_internal(e)

    def sign(self, sk, M, ctx=b''):
        if len(ctx) > 255:
            raise ValueError("the context string is too long")
        rnd = self.random_bytes(32)
        M_prime = self._int_to_bytes(0, 1) + self._int_to_bytes(len(ctx), 1) + ctx + M
        sigma = self._ml_dsa_sign_internal(sk, M_prime, rnd)
        return sigma

    def verify(self, pk, M, sigma, ctx=b''):
        if len(ctx) > 255:
            raise ValueError("the context string is too long")
        M_prime = self._int_to_bytes(0, 1) + self._int_to_bytes(len(ctx), 1) + ctx + M
        return self._ml_dsa_verify_internal(pk, M_prime, sigma)



if __name__ == "__main__":
    ml = ML_DSA(ML_DSA_87)
    runs = 1
    t = []
    for run in range(runs):
        message = os.urandom(1568)
        t0 = time.time()
        pk, sk = ml.keygen()
        print(len(pk))
        print(len(sk))
        sigma = ml.sign(sk, message)
        print(len(sigma))
        verify = ml.verify(pk, message, sigma)
        t1 = time.time()
        assert verify is True
        t.append(t1 - t0)
    print(sum(t) / len(t))
    print(sum(t))
    print(t)

