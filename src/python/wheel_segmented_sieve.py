import numpy as np
import math

# Possible offsets for numbers coprime to 30
WHEEL = np.array([1, 7, 11, 13, 17, 19, 23, 29], dtype=int)

def simple_sieve(limit):
   
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[0:2] = False
    for i in range(2, int(math.isqrt(limit)) + 1):
        if sieve[i]:
            sieve[i*i:limit+1:i] = False
    return np.flatnonzero(sieve)

def segmented_sieve_wheel(n):
    
    limit = int(math.isqrt(n)) + 1
    base_primes = simple_sieve(limit)

    # Precompute all candidates in [0, 30)
    pattern = np.zeros(30, dtype=bool)
    pattern[WHEEL] = True

    segment_size = max(limit, 10_000_000)

    primes = []

    for low in range(0, n + 1, segment_size):
        high = min(low + segment_size, n + 1)

        size = high - low
        sieve = np.zeros(size, dtype=bool)

        # Mark numbers coprime to 30 first
        for i in range(size):
            if (low + i) < 7:
                sieve[i] = (low + i) in (2, 3, 5)
            else:
                sieve[i] = pattern[(low + i) % 30]

        for p in base_primes:
            if p * p >= high:
                break
            # Find the first multiple of p in [low, high)
            start = max(p * p, ((low + p - 1) // p) * p)
            for j in range(start, high, p):
                sieve[j - low] = False

        primes.extend(low + np.flatnonzero(sieve))

    return primes

# Example usage
if __name__ == "__main__":
    n = 1000
    primes = segmented_sieve_wheel(n)
    print(primes)
