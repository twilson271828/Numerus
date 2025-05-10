import numpy as np
import matplotlib.pyplot as plt

# Sampling settings
Fs = 500              # Sampling frequency (Hz)
T = 1 / Fs            # Sampling interval
t = np.arange(0, 1, T)  # 1 second of data

# Create a noisy signal composed of two sine waves
f1 = 50   # Frequency of first sine wave (Hz)
f2 = 120  # Frequency of second sine wave (Hz)
signal = 0.7 * np.sin(2 * np.pi * f1 * t) + 0.5 * np.sin(2 * np.pi * f2 * t)
signal += 0.5 * np.random.randn(len(t))  # Add Gaussian noise

# Compute FFT
n = len(signal)
fft_result = np.fft.fft(signal)
frequencies = np.fft.fftfreq(n, d=T)

# Take the magnitude and focus only on the positive half
magnitude = np.abs(fft_result)[:n // 2]
frequencies = frequencies[:n // 2]

# Plot
plt.figure(figsize=(12, 5))

plt.subplot(1, 2, 1)
plt.plot(t, signal)
plt.title("Time Domain Signal")
plt.xlabel("Time [s]")
plt.ylabel("Amplitude")

plt.subplot(1, 2, 2)
plt.plot(frequencies, magnitude)
plt.title("Frequency Domain (FFT)")
plt.xlabel("Frequency [Hz]")
plt.ylabel("Magnitude")
plt.grid()

plt.tight_layout()
plt.show()
