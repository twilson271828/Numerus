import numpy as np
import math

def simple_sieve(limit):
    """Simple sieve to get primes up to sqrt(n), optimized for odd numbers."""
    sieve = np.ones(limit // 2, dtype=bool)  # only odd numbers: index i represents (2*i + 1)
    for i in range(1, int(math.isqrt(limit)) // 2 + 1):
        if sieve[i]:
            p = 2 * i + 1
            sieve[p*p//2::p] = False
    primes = [2] + [2 * i + 1 for i, is_prime in enumerate(sieve) if is_prime and i > 0]
    return primes

def segmented_sieve(n):
    """Segmented sieve using NumPy and odd-only optimization."""
    limit = int(math.isqrt(n)) + 1
    primes = simple_sieve(limit)

    segment_size = max(limit, 10_000_000)  # Make segments large (depends on your RAM)

    result = []

    for low in range(2, n + 1, segment_size):
        high = min(low + segment_size, n + 1)

        # Only odds: map [low, high) to odds
        offset = 0
        if low % 2 == 0:
            offset = 1
        size = (high - low + offset) // 2

        sieve = np.ones(size, dtype=bool)

        for prime in primes:
            # Find the minimum multiple of prime >= low
            start = max(prime * prime, ((low + prime - 1) // prime) * prime)
            if start % 2 == 0:
                start += prime  # ensure start is odd

            # Mark multiples of prime
            for multiple in range(start, high, 2 * prime):
                sieve[(multiple - low) // 2] = False

        # Collect primes from this segment
        if low <= 2 < high:
            result.append(2)
        result.extend([low + 2 * i + offset for i, is_prime in enumerate(sieve) if is_prime])

    return result

# Example usage
if __name__ == "__main__":
    n = 10000
    primes = segmented_sieve(n)
    print(primes)
