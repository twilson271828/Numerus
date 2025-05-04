import numpy as np
import matplotlib.pyplot as plt

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
    # Example usage
    signal = np.random.randint(0,9,16)  # Example signal
    x = fft(signal)
    print(x)
    # Plotting the original and transformed signals
    plt.subplot(2, 1, 1)
    plt.plot(signal)
    plt.title('Original Signal')
    plt.subplot(2, 1, 2)
    plt.plot(np.abs(x))
    plt.title('FFT of Signal')
    plt.show()

    