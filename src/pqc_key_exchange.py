from src.kyber import kyber as k
from src.dilithium import dilithium as d
import time
from hashlib import sha3_512
import tracemalloc
from kyber_py.ml_kem import ML_KEM_1024, ML_KEM_768, ML_KEM_512
from dilithium_py.ml_dsa import ML_DSA_44, ML_DSA_65, ML_DSA_87

def signed_key_exchange(kyber, dilithium):
    # alice generates random keys
    # public key a_pk and encapsulation key a_enc are public
    # secret key a_sk and decapsulation key a_dec are private
    a_enc, a_dec = kyber.keygen()
    a_pk, a_sk = dilithium.keygen()

    #alice sends bob the signed encapsulation key a_sigma
    a_sigma = dilithium.sign(a_sk, a_enc)

    # bob verifies the send encapsulation key from alice
    #print(f"a_sigma and a_enc in byte size that was send to bob: {len(a_sigma + a_enc)}")
    a_verify = dilithium.verify(a_pk, a_enc, a_sigma)
    assert a_verify is True

    # bob uses a_enc to get the shared key K and the cipher c
    K, c = kyber.encaps(a_enc)

    # bob generates random keys to sign his message and sends it to alice
    b_pk, b_sk = dilithium.keygen()
    b_sigma = dilithium.sign(b_sk, c)

    # alice verifies the send cipher from bob
    #print(f"b_sigma and c in byte size that was send to alice: {len(b_sigma + c)}")
    b_verify = dilithium.verify(b_pk, c, b_sigma)
    assert b_verify is True

    # alice uses the cipher to get the shared key K
    K_prime = kyber.decaps(a_dec, c)
    assert K == K_prime

    # return the session key that both use
    #print(f"byte size used to generate hash: {len(K + a_enc + a_sigma + c + b_sigma)}")
    return sha3_512(K + a_enc + a_sigma + c + b_sigma).digest()

def test_runtime(runs, kyber_parameters, dilithium_parameters):
    ky = k.ML_KEM(kyber_parameters)
    di = d.ML_DSA(dilithium_parameters)
    t = []

    for i in range(runs):
        t0 = time.time()
        signed_key_exchange(ky, di)
        t1 = time.time()
        t.append(t1 - t0)
    print(f"average key exchange time: {sum(t) / len(t)}")
    print(f"total time for {runs} key exchanges: {sum(t)}")
    print("number of iterations per second: ", 1 / (sum(t) / len(t)))
    #print(t)

def test_memory_use(kyber_parameters, dilithium_parameters):
    ky = k.ML_KEM(kyber_parameters)
    di = d.ML_DSA(dilithium_parameters)

    tracemalloc.start()
    signed_key_exchange(ky, di)
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()

def test_compare_runtime(runs, ky, di):
    t = []

    for i in range(runs):
        t0 = time.time()
        signed_key_exchange(ky, di)
        t1 = time.time()
        t.append(t1 - t0)
    print(f"average key exchange time: {sum(t) / len(t)}")
    print(f"total time for {runs} key exchanges: {sum(t)}")
    print("number of iterations per second: ", 1 / (sum(t) / len(t)))
    # print(t)

def test_compare_memory_use(ky, di):
    tracemalloc.start()
    signed_key_exchange(ky, di)
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()

if __name__ == "__main__":
    # own implementation
    #test_runtime(1000, k.ML_KEM_512, d.ML_DSA_44)
    #test_memory_use(k.ML_KEM_1024, d.ML_DSA_87)

    # compared implementation from GitHub kyber-py and dilithium-py
    test_compare_runtime(1000, ML_KEM_1024, ML_DSA_87)
    #test_compare_memory_use(ML_KEM_1024, ML_DSA_87)
