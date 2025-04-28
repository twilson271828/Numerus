import numpy as np
import math
import struct
import gzip

def simple_sieve(limit):
    
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[0:2] = False
    for i in range(2, int(math.isqrt(limit)) + 1):
        if sieve[i]:
            sieve[i*i:limit+1:i] = False
    return np.flatnonzero(sieve)

def save_deltas_binary(deltas, file_handle):
   
    data = struct.pack(f"{len(deltas)}H", *deltas)  # 'H' = unsigned short
    file_handle.write(data)

def segmented_sieve_delta_binary_gzip(n, output_filename, buffer_size=100000):
  
    limit = int(math.isqrt(n)) + 1
    base_primes = simple_sieve(limit)

    segment_size = max(limit, 10_000_000)

    with gzip.open(output_filename, "wb") as f:  # Binary gzip write
        buffer = []

        last_prime = None

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
                            buffer.append(prime)
                        else:
                            delta = prime - last_prime
                            if delta >= 65536:
                                raise ValueError(f"Delta {delta} too large for 2 bytes.")
                            buffer.append(delta)
                        last_prime = prime

                if len(buffer) >= buffer_size:
                    save_deltas_binary(buffer, f)
                    buffer.clear()

        if buffer:
            save_deltas_binary(buffer, f)



def load_binary_deltas_gzip(filename):
    
    with gzip.open(filename, "rb") as f:
        data = f.read()
    deltas = struct.unpack(f"{len(data) // 2}H", data)
    primes = []
    running_sum = 0
    for delta in deltas:
        running_sum += delta
        primes.append(running_sum)
    return primes



if __name__ == "__main__":
    n = 10000000
    
    output_file = "prime_deltas_up_to_"+str(n)+".bin.gz"

    print(f"Generating primes up to {n} with delta encoding into compressed binary {output_file}...")
    segmented_sieve_delta_binary_gzip(n, output_file)
    primes = load_binary_deltas_gzip(output_file)
    print(primes)
    print("Done.")


