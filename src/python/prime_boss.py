import numpy as np
import math
import struct
import gzip
from concurrent.futures import ThreadPoolExecutor
from queue import Queue
import threading

def simple_sieve(limit):
    """Simple sieve of Eratosthenes to get primes up to sqrt(n)."""
    sieve = np.ones(limit + 1, dtype=bool)
    sieve[0:2] = False
    for i in range(2, int(math.isqrt(limit)) + 1):
        if sieve[i]:
            sieve[i*i:limit+1:i] = False
    return np.flatnonzero(sieve)

def sieve_segment(low, high, base_primes, last_prime_info):
    """Sieve a single segment and return delta-of-delta encoded primes."""
    size = high - low
    sieve = np.ones(size, dtype=bool)

    for prime in base_primes:
        if prime * prime >= high:
            break
        start = max(prime * prime, ((low + prime - 1) // prime) * prime)
        for multiple in range(start, high, prime):
            sieve[multiple - low] = False

    buffer = []
    last_prime, last_delta = last_prime_info

    for i in range(size):
        if sieve[i]:
            prime = low + i
            if prime > 1:
                if last_prime is None:
                    buffer.append(prime)
                elif last_delta is None:
                    delta = prime - last_prime
                    buffer.append(delta)
                    last_delta = delta
                else:
                    delta = prime - last_prime
                    delta_of_delta = delta - last_delta
                    buffer.append(delta_of_delta)
                    last_delta = delta
                last_prime = prime

    return buffer, (last_prime, last_delta)

def writer_thread(filename, queue):
    """Thread that writes compressed delta-of-delta data to disk."""
    with gzip.open(filename, "wb") as f:
        while True:
            item = queue.get()
            if item is None:
                break  # Sentinel value means "done"
            data = struct.pack(f"{len(item)}h", *item)  # 'h' = signed short
            f.write(data)
            queue.task_done()

def multithreaded_sieve_and_compress(n, output_filename, num_threads=4):
    """Main function: Multithreaded segmented sieve and parallel compressed writing."""
    limit = int(math.isqrt(n)) + 1
    base_primes = simple_sieve(limit)

    segment_size = max(limit, 10_000_000)
    queue = Queue(maxsize=10)  # Buffer up to 10 segments

    writer = threading.Thread(target=writer_thread, args=(output_filename, queue))
    writer.start()

    last_prime_info = (None, None)

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = []
        for low in range(2, n + 1, segment_size):
            high = min(low + segment_size, n + 1)
            futures.append(executor.submit(sieve_segment, low, high, base_primes, last_prime_info))
            # Update last_prime_info sequentially
            if futures:
                buffer, last_prime_info = futures[-1].result()
                queue.put(buffer)

    queue.put(None)  # Signal the writer thread to finish
    writer.join()


def load_prime_delta_of_delta_gzip(filename):
    """Load delta-of-delta encoded primes from compressed binary (.gz) and reconstruct primes."""
    with gzip.open(filename, "rb") as f:
        data = f.read()
    # Each delta-of-delta was stored as a signed short ('h'), 2 bytes
    deltas = struct.unpack(f"{len(data) // 2}h", data)

    primes = []
    running_sum = 0
    last_delta = None

    for i, delta in enumerate(deltas):
        if i == 0:
            # First item is the first prime
            running_sum = delta
            primes.append(running_sum)
        elif i == 1:
            # Second item is the first delta
            running_sum += delta
            last_delta = delta
            primes.append(running_sum)
        else:
            # After that, reconstruct delta and prime
            last_delta += delta  # Update delta
            running_sum += last_delta
            primes.append(running_sum)

    return primes



if __name__ == "__main__":
    n = 10**6  # Or 10**8, 10**9, etc.
    output_file = "primes_multithreaded_delta_of_delta.bin.gz"
    num_threads = 4

    print(f"Generating primes up to {n} with {num_threads} threads into compressed binary {output_file}...")
    multithreaded_sieve_and_compress(n, output_file, num_threads=num_threads)
    primes = load_prime_delta_of_delta_gzip("primes_multithreaded_delta_of_delta.bin.gz")
    print(primes)
    print("Done.")
