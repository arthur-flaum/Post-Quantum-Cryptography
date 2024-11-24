from src.kyber import kyber as k
from src.dilithium import dilithium as d
from hashlib import sha512

def signed_key_exchange(kyber, dilithium):
    # alice generates random keys
    # public key a_pk and encapsulation key a_enc are public
    # secret key a_sk and decapsulation key a_dec are private
    a_enc, a_dec = kyber.keygen()
    a_pk, a_sk = dilithium.keygen()

    #alice sends bob the signed encapsulation key a_sigma
    a_sigma = dilithium.sign(a_sk, a_enc)

    # bob verifies the send encapsulation key from alice
    print(f"a_sigma and a_enc in byte size that was send to bob: {len(a_sigma + a_enc)}")
    a_verify = dilithium.verify(a_pk, a_enc, a_sigma)
    assert a_verify is True

    # bob uses a_enc to get the shared key K and the cipher c
    K, c = kyber.encaps(a_enc)

    # bob generates random keys to sign his message and sends it to alice
    b_pk, b_sk = dilithium.keygen()
    b_sigma = dilithium.sign(b_sk, c)

    # alice verifies the send cipher from bob
    print(f"b_sigma and c in byte size that was send to alice: {len(b_sigma + c)}")
    b_verify = dilithium.verify(b_pk, c, b_sigma)
    assert b_verify is True

    # alice uses the cipher to get the shared key K
    K_prime = kyber.decaps(a_dec, c)
    assert K == K_prime

    # return the session key that both use
    print(f"byte size used to generate hash: {len(K + a_enc + a_sigma + c + b_sigma)}")
    return sha512(K + a_enc + a_sigma + c + b_sigma).digest()





if __name__ == "__main__":
    ky = k.ML_KEM(k.ML_KEM_1024)
    di = d.ML_DSA(d.ML_DSA_87)

    session_key = signed_key_exchange(ky, di)
    print(session_key)
    print(len(session_key))

