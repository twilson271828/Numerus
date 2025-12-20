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
    C = np.array(A) * np.array(B)
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

    x = [even[k] + T[k] for k in range(N // 2)] + \
           [even[k] - T[k] for k in range(N // 2)]
    
    

    return x


import numpy as np

def fft1(x):
    """
    Compute the Fast Fourier Transform of a sequence x using the
    Cooley–Tukey recursive algorithm.
    
    Parameters:
        x (list or np.ndarray): A list of complex numbers, length must be a power of 2.

    Returns:
        np.ndarray: The FFT of the input sequence.
    """
    x = np.asarray(x, dtype=complex)
    n = x.shape[0]

    if n == 1:
        return x
    elif n % 2 != 0:
        raise ValueError("Size of x must be a power of 2")

    # Split into even and odd terms
    X_even = fft(x[::2])
    X_odd = fft(x[1::2])
    
    # Twiddle factors
    factor = np.exp(-2j * np.pi * np.arange(n) / n)
    return np.concatenate([
        X_even + factor[:n // 2] * X_odd,
        X_even - factor[:n // 2] * X_odd
    ])


def ifft1(X):
    """
    Compute the Inverse Fast Fourier Transform using the Cooley–Tukey algorithm.
    
    Parameters:
        X (list or np.ndarray): A list of complex numbers, length must be a power of 2.

    Returns:
        np.ndarray: The inverse FFT of the input sequence.
    """
    X = np.asarray(X, dtype=complex)
    n = X.shape[0]

    if n == 1:
        return X
    elif n % 2 != 0:
        raise ValueError("Size of X must be a power of 2")

    # Split into even and odd terms
    X_even = ifft(X[::2])
    X_odd = ifft(X[1::2])

    # Twiddle factors (conjugated)
    factor = np.exp(2j * np.pi * np.arange(n) / n)
    result = np.concatenate([
        X_even + factor[:n // 2] * X_odd,
        X_even - factor[:n // 2] * X_odd
    ])
    
    return result / 2  # Normalize (recursive step doubles size)


if __name__=="__main__":


    # Sample data (e.g., sine wave)
    N = 2048 # Length of the signal (must be a power of 2)
    t = np.linspace(0, 1, N, endpoint=False)
    signal = np.sin(2 * np.pi * 3 * t)  # 3 Hz sine wave

    fft_result = fft1(signal)
    fft_real_result = fft_result.real.round().astype(np.int64)
    print(fft_real_result)
    np_fft_result = np.fft.fft(signal,N).real.round().astype(np.int64)
    
    print("FFT result:", fft_real_result)

     # Plotting the origina and transformed signals
    plt.subplot(4, 1, 1)
    plt.plot(fft_real_result)
    plt.title('Open AI FFT')
    plt.subplot(4, 1, 2)
    plt.plot(np_fft_result)
    plt.title('Numpy FFT')
    plt.subplot(4, 1, 3)
    plt.plot(signal)
    plt.title('Numpy FFT')
    plt.subplot(4, 1, 4)
    plt.plot(ifft1(np.real(fft_result)))
    plt.title('Numpy FFT')

    plt.show()

    


    """
    a = 12345678
    b = 98765432
    product = schonhage_strassen(a, b)
    print(product)
    print(a*b)
    print("Correct:", product == a * b)

    

    a =9876543231243242342
    
    base = 10**2 # use 4-digit chunks to reduce FFT size
    asplit = split_number(a, base)

    print(asplit)
   
    n1 = 1 << (len(asplit)-1).bit_length()
    
    print("n1",n1)

    x1 = fft1(asplit)

    x2=np.fft.fft(asplit,n1)
    
    # Plotting the original and transformed signals
    plt.subplot(3, 1, 1)
    plt.plot(x1)
    plt.title('Open AI FFT')
    plt.subplot(3, 1, 2)
    plt.plot(x2)
    plt.title('Numpy FFT')

    #plt.subplot(3,1,3)
    #plt.plot(np.abs(x1-x3))   
    #plt.title('Difference between Open AI FFT and Numpy FFT')
    #plt.tight_layout()
    
    plt.show()
    
    """
   

    