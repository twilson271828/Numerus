import numpy as np
import math
from concurrent.futures import ThreadPoolExecutor, as_completed

def simple_sieve(limit):
    """Simple sieve of Eratosthenes to get primes up to sqrt(n)."""
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[0:2] = False
    for i in range(2, int(math.isqrt(limit)) + 1):
        if sieve[i]:
            sieve[i*i:limit+1:i] = False
    return np.flatnonzero(sieve)

def sieve_segment(low, high, base_primes):
    """Sieve a single segment from low to high."""
    size = high - low
    sieve = np.ones(size, dtype=bool)

    for prime in base_primes:
        if prime * prime >= high:
            break
        start = max(prime * prime, ((low + prime - 1) // prime) * prime)
        for multiple in range(start, high, prime):
            sieve[multiple - low] = False

    segment_primes = [i for i in range(low, high) if sieve[i - low] and i > 1]
    return segment_primes

def segmented_sieve_parallel(n, num_threads=4):
    """Parallel segmented sieve up to n using num_threads threads."""
    limit = int(math.isqrt(n)) + 1
    base_primes = simple_sieve(limit)

    segment_size = max(limit, 10_000_000)

    futures = []
    primes = []

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for low in range(2, n + 1, segment_size):
            high = min(low + segment_size, n + 1)
            futures.append(executor.submit(sieve_segment, low, high, base_primes))

        for future in as_completed(futures):
            primes.extend(future.result())

    # Sort the primes because segments may complete in any order
    primes.sort()
    return primes

# Example usage
if __name__ == "__main__":
    n = 100000
    primes = segmented_sieve_parallel(n, num_threads=4)
    print(primes)
