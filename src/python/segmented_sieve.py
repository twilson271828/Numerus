import math

def simple_sieve(limit):
    """Simple sieve of Eratosthenes up to limit (√n)."""
    is_prime = [True] * (limit + 1)
    is_prime[0:2] = [False, False]
    for num in range(2, int(math.isqrt(limit)) + 1):
        if is_prime[num]:
            for multiple in range(num * num, limit + 1, num):
                is_prime[multiple] = False
    return [num for num, prime in enumerate(is_prime) if prime]

def segmented_sieve(n):
    """Segmented Sieve of Eratosthenes up to n."""
    limit = int(math.isqrt(n)) + 1
    primes = simple_sieve(limit)  # Small primes up to √n

    # Set segment size. You can experiment with it for faster performance.
    segment_size = max(limit, 10_000)

    # Resulting list of primes
    result = []

    # Process segments
    low = 2
    high = low + segment_size
    while low <= n:
        if high > n + 1:
            high = n + 1

        # Create a boolean array for this segment
        is_prime = [True] * (high - low)

        for prime in primes:
            # Find the minimum number in [low, high) that is a multiple of prime
            start = max(prime * prime, ((low + prime - 1) // prime) * prime)

            for multiple in range(start, high, prime):
                is_prime[multiple - low] = False

        # Collect primes in this segment
        for i in range(low, high):
            if is_prime[i - low]:
                result.append(i)

        # Move to the next segment
        low += segment_size
        high += segment_size

    return result

# Example usage
if __name__ == "__main__":
    n = 10000
    primes = segmented_sieve(n)
    print(primes)
