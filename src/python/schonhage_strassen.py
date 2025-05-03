import numpy as np

def split_number(n, base):
    """Split integer n into list of digits in given base (least significant digit first)."""
    digits = []
    while n:
        digits.append(n % base)
        n //= base
    return digits if digits else [0]

def combine_digits(digits, base):
    """Recombine digits to integer, handling carry."""
    carry = 0
    for i in range(len(digits)):
        digits[i] += carry
        carry = digits[i] // base
        digits[i] %= base
    while carry:
        digits.append(carry % base)
        carry //= base
    # remove leading zeros
    while len(digits) > 1 and digits[-1] == 0:
        digits.pop()
    return digits

def schonhage_strassen(x, y):
    """Multiply two integers using FFT (simulated Schönhage-Strassen style)."""
    base = 10**4  # use 4-digit chunks to reduce FFT size
    a = split_number(x, base)
    b = split_number(y, base)

    n = 1 << (len(a) + len(b) - 1).bit_length()
    A = np.fft.fft(a, n)
    B = np.fft.fft(b, n)
    C = A * B
    c = np.fft.ifft(C).real.round().astype(np.int64)

    # propagate carries
    c = combine_digits(list(c), base)

    # convert back to Python int
    result = 0
    for digit in reversed(c):
        result = result * base + digit
    return result


if __name__ == "__main__":
    a = 12345678
    b = 98765432
    product = schonhage_strassen(a, b)
    print(product)
    print("Correct:", product == a * b)