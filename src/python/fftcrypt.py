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
    signal = np.random.randint(0,9,128)  # Example signal
    x = fft(signal)
    n = len(signal)
    A = np.fft.fft(signal, n)
    

    
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

    