import collections
import hashlib
import random
import time

EllipticCurve = collections.namedtuple('EllipticCurve', 'name p a b g n h')
# elliptic curve formula: y**2 = x**3 + a*x + b
curve = EllipticCurve(
    'secp256k1',
    # Field characteristic.
    p=0xfffffffffffffffffffffffffffffffffffffffffffffffffffffffefffffc2f,
    # Curve coefficients.
    a=0,
    b=7,
    # Base point.
    g=(0x79be667ef9dcbbac55a06295ce870b07029bfcdb2dce28d959f2815b16f81798,
       0x483ada7726a3c4655da4fbfc0e1108a8fd17b448a68554199c47d08ffb10d4b8),
    # Subgroup order.
    n=0xfffffffffffffffffffffffffffffffebaaedce6af48a03bbfd25e8cd0364141,
    # Subgroup cofactor.
    h=1,
)


# Modular arithmetic

def inverse_mod(k, p):
    """Returns the inverse of k modulo p.

    This function returns the only integer x such that (x * k) % p == 1.

    k must be non-zero and p must be a prime.
    """
    if k == 0:
        raise ZeroDivisionError('division by zero')

    if k < 0:
        # k ** -1 = p - (-k) ** -1  (mod p)
        return p - inverse_mod(-k, p)

    # Extended Euclidean algorithm.
    s, old_s = 0, 1
    t, old_t = 1, 0
    r, old_r = p, k

    while r != 0:
        quotient = old_r // r
        old_r, r = r, old_r - quotient * r
        old_s, s = s, old_s - quotient * s
        old_t, t = t, old_t - quotient * t

    gcd, x, y = old_r, old_s, old_t

    assert gcd == 1
    assert (k * x) % p == 1

    return x % p


# Functions that work on curve points

def is_on_curve(point):
    """Returns True if the given point lies on the elliptic curve."""
    if point is None:
        # None represents the point at infinity.
        return True

    x, y = point

    return (y * y - x * x * x - curve.a * x - curve.b) % curve.p == 0


def point_neg(point):
    """Returns -point."""
    assert is_on_curve(point)

    if point is None:
        # -0 = 0
        return None

    x, y = point
    result = (x, -y % curve.p)

    assert is_on_curve(result)

    return result


def point_add(point1, point2):
    """Returns the result of point1 + point2 according to the group law."""
    assert is_on_curve(point1)
    assert is_on_curve(point2)

    if point1 is None:
        # 0 + point2 = point2
        return point2
    if point2 is None:
        # point1 + 0 = point1
        return point1

    x1, y1 = point1
    x2, y2 = point2

    if x1 == x2 and y1 != y2:
        # point1 + (-point1) = 0
        return None

    if x1 == x2:
        # This is the case point1 == point2.
        m = (3 * x1 * x1 + curve.a) * inverse_mod(2 * y1, curve.p)
    else:
        # This is the case point1 != point2.
        m = (y1 - y2) * inverse_mod(x1 - x2, curve.p)

    x3 = m * m - x1 - x2
    y3 = y1 + m * (x3 - x1)
    result = (x3 % curve.p,
              -y3 % curve.p)

    assert is_on_curve(result)

    return result


def scalar_mult(k, point):
    """Returns k * point computed using the double and point_add algorithm."""
    assert is_on_curve(point)

    if k % curve.n == 0 or point is None:
        return None

    if k < 0:
        # k * point = -k * (-point)
        return scalar_mult(-k, point_neg(point))

    result = None
    addend = point

    while k:
        if k & 1:
            # Add.
            result = point_add(result, addend)

        # Double.
        addend = point_add(addend, addend)

        k >>= 1

    assert is_on_curve(result)

    return result


# Keypair generation and ECDSA

def make_keypair():
    """Generates a random private-public key pair."""
    private_key = random.randrange(1, curve.n)
    public_key = scalar_mult(private_key, curve.g)

    return private_key, public_key


def hash_message(message):
    """Returns the truncated SHA521 hash of the message."""
    message_hash = hashlib.sha512(message).digest()
    e = int.from_bytes(message_hash, 'big')

    # FIPS 180 says that when a hash needs to be truncated, the rightmost bits
    # should be discarded.
    z = e >> (e.bit_length() - curve.n.bit_length())

    assert z.bit_length() <= curve.n.bit_length()

    return z


def sign_message(private_key, message):
    z = hash_message(message)

    r = 0
    s = 0

    while not r or not s:
        k = random.randrange(1, curve.n)
        x, y = scalar_mult(k, curve.g)

        r = x % curve.n
        s = ((z + r * private_key) * inverse_mod(k, curve.n)) % curve.n

    return r, s


def verify_signature(public_key, message, signature):
    z = hash_message(message)

    r, s = signature

    w = inverse_mod(s, curve.n)
    u1 = (z * w) % curve.n
    u2 = (r * w) % curve.n

    x, y = point_add(scalar_mult(u1, curve.g),
                     scalar_mult(u2, public_key))

    if (r % curve.n) == (x % curve.n):
        #print('signature matches')
        return public_key
    else:
        #print('invalid signature')
        return None


def sign_ecdh_key_exchange(secret_key, public_key):
    """
    :param secret_key: private key from the signer
    :param public_key: base point g multiplied with the private key (only x value)
    :return: x value of public_key and the hash of it encrypted with the secret key
    """
    message = str(public_key).encode("utf-8")
    return sign_message(secret_key, message)


def verify_ecdh_key_exchange(signature, public_key):
    """
    :param signature: signed message
    :param public_key: base point g multiplied with the private key
    :return: x value of public key if signature matches
    """
    message = str(public_key[0]).encode("utf-8")
    return verify_signature(public_key, message, signature)


def signed_dh_key_exchange():
    # a -> alice, b -> bob
    a_priv, a_pub = make_keypair()
    b_priv, b_pub = make_keypair()
    a_pub_x, a_pub_y = a_pub
    b_pub_x, b_pub_y = b_pub

    a_signed_key = sign_ecdh_key_exchange(a_priv, a_pub_x)
    a_verified_key = verify_ecdh_key_exchange(a_signed_key, a_pub)
    assert a_pub == a_verified_key
    b_secret_exchanged_key = scalar_mult(b_priv, a_verified_key)
    """
    print(f"alice public key x value: {hex(a_pub_x)}")
    print(f"random number multiplied with the x value of the base point mod n from elliptic curve: {hex(a_signed_key[0])} \nhashed and signed key from alice: {hex(a_signed_key[1])}")
    print(f"verified key bob got from alice: ({hex(a_verified_key[0])}, {hex(a_verified_key[1])})")
    print(f"the shared x value secret key between alice and bob: {hex(b_secret_exchanged_key[0])}")
    """

    b_signed_key = sign_ecdh_key_exchange(b_priv, b_pub_x)
    b_verified_key = verify_ecdh_key_exchange(b_signed_key, b_pub)
    assert b_pub == b_verified_key
    a_secret_exchanged_key = scalar_mult(a_priv, b_verified_key)
    """
    print(f"\nbob public key x value: {hex(b_pub_x)}")
    print(f"random number multiplied with the x value of the base point mod n from elliptic curve: {hex(b_signed_key[0])} \nhashed and signed key from bob: {hex(b_signed_key[1])}")
    print(f"verified key alice got from bob: ({hex(b_verified_key[0])}, {hex(b_verified_key[1])})")
    print(f"the shared x value secret key between alice and bob: {hex(a_secret_exchanged_key[0])}")

    print(f"\nshared x value secret key between alice and bob are equal: {a_secret_exchanged_key == b_secret_exchanged_key}")
    """
    assert a_secret_exchanged_key == b_secret_exchanged_key
    return a_secret_exchanged_key


if __name__ == "__main__":
    runs = 1000
    t0 = time.time()
    for i in range(runs):
        signed_dh_key_exchange()
    t1 = time.time()
    print(f"time elapsed: {t1 - t0}")
    print(f"average time elapsed: {(t1 - t0) / runs}")