from dilithium_py.ml_dsa import ML_DSA_44, ML_DSA_65, ML_DSA_87
import time
import tracemalloc
import os

def test_runtime(runs, ml):
    total_time = []
    keygen_time = []
    sign_time = []
    verify_time = []
    for run in range(runs):
        message = os.urandom(1568)
        t0 = time.time()
        pk, sk = ml.keygen()
        t1 = time.time()
        sig = ml.sign(sk, message)
        t2 = time.time()
        verify = ml.verify(pk, message, sig)
        t3 = time.time()
        assert verify is True
        total_time.append(t3 - t0)
        keygen_time.append(t1 - t0)
        sign_time.append(t2 - t1)
        verify_time.append(t3 - t2)
    print("NUMBER OF ITERATIONS: ", runs)
    print("\n######## TOTAL TIME ########")
    print("average time: ", sum(total_time) / len(total_time))
    print("total time: ", sum(total_time))
    print("number of iterations per second: ", 1 / (sum(total_time) / len(total_time)))
    # print(total_time)

    print("\n######## KEY GENERATION TIME ########")
    print("average time: ", sum(keygen_time) / len(keygen_time))
    print("total time: ", sum(keygen_time))
    print("number of iterations per second: ", 1 / (sum(keygen_time) / len(keygen_time)))
    # print(keygen_time)

    print("\n######## SIGNING TIME ########")
    print("average time: ", sum(sign_time) / len(sign_time))
    print("total time: ", sum(sign_time))
    print("number of iterations per second: ", 1 / (sum(sign_time) / len(sign_time)))
    # print(sign_time)

    print("\n######## VERIFYING TIME ########")
    print("average time: ", sum(verify_time) / len(verify_time))
    print("total time: ", sum(verify_time))
    print("number of iterations per second: ", 1 / (sum(verify_time) / len(verify_time)))
    # print(verify_time)

def test_memory_use(ml):
    message = os.urandom(1568)
    tracemalloc.start()
    pk, sk = ml.keygen()
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()

    tracemalloc.start()
    sig = ml.sign(sk, message)
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()

    tracemalloc.start()
    verify = ml.verify(pk, message, sig)
    assert verify is True
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()

if __name__ == "__main__":
    test_runtime(1000, ML_DSA_87)
    #test_memory_use(ML_DSA_87)

