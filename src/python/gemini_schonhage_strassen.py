import math


class SchonhageStrassen:
    def __init__(self):
        pass

    def _fft(self, a, invert, n, w):
        """
        Performs the Fast Fourier Transform (or Inverse FFT) over a ring.
        """
        if n == 1:
            return a

        a_even = a[0::2]
        a_odd = a[1::2]
        
        # Recursive calls
        # Note: We square w for the sub-problems because the order is halved
        y_even = self._fft(a_even, invert, n // 2, (w * w)) 
        y_odd = self._fft(a_odd, invert, n // 2, (w * w))

        y = [0] * n
        w_curr = 1
        
        # Combine results (Butterfly operations)
        for i in range(n // 2):
            val = w_curr * y_odd[i]
            y[i] = y_even[i] + val
            y[i + (n // 2)] = y_even[i] - val
            w_curr *= w
        
        return y

    def multiply(self, x, y):
        """
        Main entry point for multiplication.
        """
        # 1. Handle base cases and signs
        sign = 1
        if x < 0: x = -x; sign = -sign
        if y < 0: y = -y; sign = -sign
        
        if x == 0 or y == 0: return 0
        
        # 2. Determine size for FFT
        # Length of result in bits
        n_bits = x.bit_length() + y.bit_length()
        
        # Find smallest power of 2 greater than n_bits
        # We need N to be a power of 2 for standard FFT
        N = 1
        while N < n_bits:
            N *= 2
            
        # To avoid floating point issues and simulate the ring structure strictly,
        # standard Schönhage-Strassen uses Fermat numbers. 
        # However, for a pure Python demonstration without a custom BigInt library
        # for modular arithmetic, we can use a simpler approach using 
        # complex number FFT (standard float-based) and rounding, 
        # OR we can implement the Number Theoretic Transform (NTT) with a prime modulus.
        
        # Below is the NTT approach using a prime modulus P = c * 2^k + 1
        # which supports the transform size we need.
        
        # 3. Setup NTT Parameters
        # We need a prime P such that P = c * 2^k + 1 where 2^k >= N
        # A common choice for large N is creating a large prime or using specific constants.
        # For simplicity in this demo, we use Python's large integer support 
        # and a standard FFT logic adapted for large integers.
        
        # Let's break the numbers into chunks to treat them as polynomials.
        # We process 'bits_per_chunk' bits at a time.
        bits_per_chunk = max(1, n_bits.bit_length() // 2) # Heuristic size
        
        # Convert integers to list of "coefficients" (digits in base 2^bits_per_chunk)
        A = self._int_to_chunks(x, bits_per_chunk)
        B = self._int_to_chunks(y, bits_per_chunk)
        
        # Resize A and B to size M (power of 2)
        M = 1
        while M < len(A) + len(B):
            M *= 2
            
        A += [0] * (M - len(A))
        B += [0] * (M - len(B))
        
        # 4. Perform FFT (using complex numbers for standard demonstration ease)
        # Note: A strict SSA implementation uses modulo (2^(2^n)+1). 
        # This implementation uses standard complex FFT for clarity on the 
        # "Convolution Theorem" aspect, which is the heart of the algorithm.
        
        fft_a = self._recursive_fft(A, False)
        fft_b = self._recursive_fft(B, False)
        
        # 5. Point-wise Multiplication
        fft_c = [fft_a[i] * fft_b[i] for i in range(M)]
        
        # 6. Inverse FFT
        C_raw = self._recursive_fft(fft_c, True)
        
        # 7. Process Carries
        # Because we used complex FFT, we round real parts to nearest integer.
        result = 0
        carry = 0
        
        for i in range(len(C_raw)):
            val = round(C_raw[i].real) + carry
            
            # The value at this position is val % (2^bits_per_chunk)
            digit = val % (1 << bits_per_chunk)
            carry = val // (1 << bits_per_chunk)
            
            result += digit * (1 << (i * bits_per_chunk))
            
        return result * sign

    def _int_to_chunks(self, num, chunk_size):
        """Splits an integer into chunks of bits."""
        chunks = []
        mask = (1 << chunk_size) - 1
        while num > 0:
            chunks.append(num & mask)
            num >>= chunk_size
        return chunks

    def _recursive_fft(self, a, invert):
        """
        Standard Recursive FFT implementation (Cooley-Tukey).
        """
        n = len(a)
        if n == 1:
            return a
        
        a0 = [a[i] for i in range(0, n, 2)]
        a1 = [a[i] for i in range(1, n, 2)]
        
        self._recursive_fft(a0, invert)
        self._recursive_fft(a1, invert)
        
        # Prepare roots of unity
        angle = 2 * math.pi / n * (-1 if invert else 1)
        w = complex(1, 0)
        wn = complex(math.cos(angle), math.sin(angle))
        
        for i in range(n // 2):
            a[i] = a0[i] + w * a1[i]
            a[i + n // 2] = a0[i] - w * a1[i]
            if invert:
                a[i] /= 2
                a[i + n // 2] /= 2
            w *= wn
            
        return a

# --- Usage ---
if __name__ == "__main__":
    ssa = SchonhageStrassen()
    
    # Test with two large numbers
    num1 = 123456789012345678901234567890
    num2 = 987654321098765432109876543210
    
    print(f"Number 1: {num1}")
    print(f"Number 2: {num2}")
    
    result = ssa.multiply(num1, num2)
    expected = num1 * num2
    
    print("-" * 20)
    print(f"SSA Result:     {result}")
    print(f"Python Native:  {expected}")
    print(f"Match:          {result == expected}")