import numpy as np
import math

def simple_sieve(limit):
    """Simple sieve to generate primes up to sqrt(n)."""
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[0:2] = False
    for i in range(2, int(math.isqrt(limit)) + 1):
        if sieve[i]:
            sieve[i*i:limit+1:i] = False
    return np.flatnonzero(sieve)

def save_primes_to_file(primes, file_handle):
    """Save a list of primes to a file efficiently."""
    for prime in primes:
        file_handle.write(f"{prime}\n")  # write each prime on a new line

def segmented_sieve_large(n, output_filename, buffer_size=100000):
    """Segmented sieve saving primes immediately to disk, handling very large n."""
    limit = int(math.isqrt(n)) + 1
    base_primes = simple_sieve(limit)

    segment_size = max(limit, 10_000_000)

    with open(output_filename, "w") as f:
        buffer = []

        for low in range(2, n + 1, segment_size):
            high = min(low + segment_size, n + 1)

            size = high - low
            sieve = np.ones(size, dtype=bool)

            for prime in base_primes:
                if prime * prime >= high:
                    break
                start = max(prime * prime, ((low + prime - 1) // prime) * prime)
                for multiple in range(start, high, prime):
                    sieve[multiple - low] = False

            # Collect primes from this segment into the buffer
            for i in range(size):
                if sieve[i]:
                    prime = low + i
                    if prime > 1:
                        buffer.append(prime)

                # Flush buffer if it's big enough
                if len(buffer) >= buffer_size:
                    save_primes_to_file(buffer, f)
                    buffer.clear()

        # Save any remaining primes
        if buffer:
            save_primes_to_file(buffer, f)

if __name__ == "__main__":
    n = 10**7   # You can change this to 10**9, 10**12 etc.
    output_file = "primes_up_to_10million.txt"

    print(f"Generating primes up to {n} and saving to {output_file}...")
    segmented_sieve_large(n, output_file)
    print("Done.")
