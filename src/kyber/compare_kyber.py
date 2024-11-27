from kyber_py.ml_kem import ML_KEM_1024, ML_KEM_768, ML_KEM_512
import time
import tracemalloc

def test_runtime(runs, ml):
    total_time = []
    keygen_time = []
    encaps_time = []
    decaps_time = []
    for run in range(runs):
        t0 = time.time()
        ek, dk = ml.keygen()
        t1 = time.time()
        k, ct = ml.encaps(ek)
        t2 = time.time()
        k_prime = ml.decaps(dk, ct)
        t3 = time.time()
        assert k == k_prime
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

def test_memory_use(ml):
    tracemalloc.start()
    ek, dk = ml.keygen()
    print(tracemalloc.get_traced_memory())
    #tracemalloc.stop()

    #tracemalloc.start()
    k, ct = ml.encaps(ek)
    print(tracemalloc.get_traced_memory())
    #tracemalloc.stop()

    #tracemalloc.start()
    k_prime = ml.decaps(dk, ct)
    assert k == k_prime
    print(tracemalloc.get_traced_memory())
    tracemalloc.stop()

if __name__ == "__main__":
    test_runtime(1000, ML_KEM_1024)
    #test_memory_use(ML_KEM_1024)


