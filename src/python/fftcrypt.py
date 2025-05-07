import numpy as np
import matplotlib.pyplot as plt


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
    A = fft(a)
    B = fft(b)
    C = A * B
    c = ifft(C)

    c = np.real(c).round().astype(np.int64)

    # propagate carries
    c = combine_digits(list(c), base)

    # convert back to Python int
    result = 0
    for digit in reversed(c):
        result = result * base + digit
    return result

def ifft(signal):
    """
    Perform IFFT on the input signal.
    """
    N = len(signal)
    if N <= 1:
        return signal
    even = ifft(signal[0::2])
    odd = ifft(signal[1::2])
    T = [np.exp(2j * np.pi * k / N) * odd[k] for k in range(N // 2)]
    return [(even[k] + T[k]) / 2 for k in range(N // 2)] + \
           [(even[k] - T[k]) / 2 for k in range(N // 2)]
   

def fft(signal):
    """
    Perform FFT on the input signal.
    """
    N = len(signal)
    if N <= 1:
        return signal
    even = fft(signal[0::2])
    odd = fft(signal[1::2])
    T = [np.exp(-2j * np.pi * k / N) * odd[k] for k in range(N // 2)]
    return [even[k] + T[k] for k in range(N // 2)] + \
           [even[k] - T[k] for k in range(N // 2)]



if __name__=="__main__":

    a = 12345678
    b = 98765432
    product = schonhage_strassen(a, b)
    print(product)
    print("Correct:", product == a * b)



    
    # Plotting the original and transformed signals
    plt.subplot(3, 1, 1)
    plt.plot(x)
    plt.title('Open AI FFT')
    plt.subplot(3, 1, 2)
    plt.plot(A)
    plt.title('Numpy FFT')

    plt.subplot(3,1,3)
    plt.plot(np.abs(x-A))   
    plt.title('Difference between Open AI FFT and Numpy FFT')
    plt.tight_layout()
    
    plt.show()

    """
   
    """
    