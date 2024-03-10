
import numpy as np

def karatsuba(x, y):
    # Base case for the recursion
    if x < 10 or y < 10:
        return x * y

    # Calculate the size of the numbers
    n = max(len(str(x)), len(str(y)))
    m = n // 2

    # Split the input numbers
    a, b = divmod(x, 10**m)
    c, d = divmod(y, 10**m)

    # Recursive steps
    ac = karatsuba(a, c)
    bd = karatsuba(b, d)
    abcd = karatsuba(a+b, c+d)

    # Combine the results
    return ac * 10**(2*m) + (abcd - ac - bd) * 10**m + bd


if __name__=="__main__":
    # Test the function
    print("1234*5678 = ", 1234*5678)
    print(karatsuba(1234, 5678))  #7006652
    print("12345678*98765432=", 12345678*98765432)
    print(karatsuba(12345678, 98765432))  #1219326221002896
    print("123*5432=", 123*5432)
    print(karatsuba(123,5432)) #668136

