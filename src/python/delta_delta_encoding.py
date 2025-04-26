import numpy as np
import math
import struct
import gzip

def simple_sieve(limit):
    """Simple sieve to generate primes up to sqrt(n)."""
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[0:2] = False
    for i in range(2, int(math.isqrt(limit)) + 1):
        if sieve[i]:
            sieve[i*i:limit+1:i] = False
    return np.flatnonzero(sieve)

def save_deltas_binary(deltas, file_handle):
    """Save a list of signed deltas to a binary file."""
    data = struct.pack(f"{len(deltas)}h", *deltas)  # 'h' = signed short (2 bytes)
    file_handle.write(data)

def segmented_sieve_delta_of_delta_gzip(n, output_filename, buffer_size=100000):
    """Segmented sieve saving delta-of-delta encoded primes compressed to disk."""
    limit = int(math.isqrt(n)) + 1
    base_primes = simple_sieve(limit)

    segment_size = max(limit, 10_000_000)

    with gzip.open(output_filename, "wb") as f:
        buffer = []

        last_prime = None
        last_delta = None

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

            for i in range(size):
                if sieve[i]:
                    prime = low + i
                    if prime > 1:
                        if last_prime is None:
                            # First prime: store it directly
                            buffer.append(prime)
                        elif last_delta is None:
                            # First delta
                            delta = prime - last_prime
                            buffer.append(delta)
                            last_delta = delta
                        else:
                            delta = prime - last_prime
                            delta_of_delta = delta - last_delta
                            if not -32768 <= delta_of_delta <= 32767:
                                raise ValueError(f"Delta-of-delta {delta_of_delta} out of range for 2 bytes.")
                            buffer.append(delta_of_delta)
                            last_delta = delta
                        last_prime = prime

                if len(buffer) >= buffer_size:
                    save_deltas_binary(buffer, f)
                    buffer.clear()

        if buffer:
            save_deltas_binary(buffer, f)


def load_prime_delta_of_delta_gzip(filename):
    """Load delta-of-delta encoded primes from compressed binary."""
    with gzip.open(filename, "rb") as f:
        data = f.read()
    deltas = struct.unpack(f"{len(data) // 2}h", data)  # signed short
    primes = []
    running_sum = 0
    last_delta = None
    for i, delta in enumerate(deltas):
        if i == 0:
            running_sum = delta
            primes.append(running_sum)
        elif i == 1:
            running_sum += delta
            last_delta = delta
            primes.append(running_sum)
        else:
            last_delta = last_delta + delta
            running_sum += last_delta
            primes.append(running_sum)
    return primes






if __name__ == "__main__":
    n = 100
    output_file = "prime_delta_of_delta_up_to_"+str(n)+".bin.gz"

    print(f"Generating primes up to {n} with delta-of-delta encoding into compressed binary {output_file}...")
    segmented_sieve_delta_of_delta_gzip(n, output_file)

    # Example usage:
    primes = load_prime_delta_of_delta_gzip(output_file)
    print(primes)
    print("Done.")
