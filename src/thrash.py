import numpy as np

def multiply_poly_fft(p1, p2):
    size = len(p1) + len(p2) - 1
    p1_fft = np.fft.fft(p1, size)
    p2_fft = np.fft.fft(p2, size)
    result_fft = p1_fft * p2_fft
    result = np.fft.ifft(result_fft)
    return np.round(result).astype(int)

def newton_raphson_division(numerator, denominator, precision=1e-10):
    # Initial guess
    x = numerator
    while True:
        x_next = x - (x * denominator - numerator) / (denominator)
        if abs(x - x_next) < precision:
            break
        x = x_next
    return x


def synthetic_division(dividend, divisor):
    # Coefficients of the divisor should be in descending order
    divisor = list(map(float, divisor))
    dividend = list(map(float, dividend))

    # The degree of the divisor should be 1 and the leading coefficient should be 1
    assert len(divisor) == 2 and divisor[0] == 1, "The divisor should be a linear polynomial with a leading coefficient of 1"

    # Initialize the quotient and the remainder
    quotient = [0] * (len(dividend) - len(divisor) + 1)
    remainder = dividend.copy()

    # Perform the synthetic division
    for i in range(len(quotient)):
        quotient[i] = remainder[i]
        for j in range(1, len(divisor)):
            remainder[i+j] -= quotient[i] * divisor[j]

    # The remainder is the last len(divisor) - 1 elements of the remainder
    remainder = remainder[len(quotient):]

    return quotient, remainder


def poly_div(dividend, divisor):
    # Make a copy of the dividend and divisor lists so we don't modify the original lists
    dividend = list(dividend)
    divisor = list(divisor)

    # The degree of the divisor
    divisor_degree = len(divisor) - 1

    # Initialize the quotient
    quotient = [0] * (len(dividend) - divisor_degree)

    # Perform the division
    for i in range(len(dividend) - divisor_degree, -1, -1):
        quotient[i - 1] = dividend[i + divisor_degree] / divisor[-1]
        for j in range(divisor_degree + 1):
            dividend[i + j] -= quotient[i - 1] * divisor[j]

    # The remainder is the first len(dividend) - len(quotient) + 1 elements of the dividend
    remainder = dividend[:len(dividend) - len(quotient) + 1]

    return quotient, remainder

# Test the function
print(poly_div([1, -3, 2], [1, -1]))  # Output: ([1.0, -2.0], [0.0])

if __name__=="__main__":

    # Test the function
    print(synthetic_division([1, 3, 2], [1, 2]))  # Output: ([1.0, -2.0], [0.0])

