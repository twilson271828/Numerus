import numpy as np

import numpy as np

def split_number_np_decimal(x: np.ndarray, m: int):
    """
    Splits a large number `x` (stored as an array of digits between 0 and 9) into parts of size `m` digits.
    Returns the list of parts.
    """
    parts = []
    total_value = 0

    # Convert the numpy array of digits (base 10) into a large integer
    for digit in x:
        total_value = total_value * 10 + int(digit)  # Combine digits into a large number

    # Now split the total_value into parts of m digits
    while total_value > 0:
        part = total_value % (10**m)  # Extract the least significant m digits
        parts.insert(0, part)  # Insert the part at the beginning of the list
        total_value //= 10**m  # Remove the m digits we just extracted

    return np.array(parts)


def burnikel_ziegler_divide(u: np.ndarray, v: np.ndarray, m: int):
    """
    Burnikel-Ziegler division algorithm.
    Divides `u` by `v` where u and v are large numbers represented as arrays of digits (0-9).
    Returns the quotient and remainder.
    """
    if v.size == 0 or np.all(v == 0):
        raise ZeroDivisionError("Division by zero is undefined.")
    
    # Convert arrays of digits into parts
    u_parts = split_number_np_decimal(u, m)
    v_parts = split_number_np_decimal(v, m)

    # Initial remainder and quotient
    quotient = []
    remainder = 0

    # Burnikel-Ziegler division steps:
    # Step 1: Perform divide-and-conquer division on parts
    for part in u_parts:
        remainder = remainder * (10**m) + part  # Shift remainder and add next part
        quotient_part = remainder // int("".join(map(str, v)))  # Divide remainder by v
        remainder = remainder % int("".join(map(str, v)))  # Get the new remainder
        quotient.append(quotient_part)

    return np.array(quotient), remainder


if __name__=="__main__":
    u = np.array([9, 9, 8, 7, 6, 5, 4, 3, 2, 1], dtype=np.uint8)  # Dividend: 9987654321
    v = np.array([9, 8, 7,3,4], dtype=np.uint8)                      # Divisor: 987
    m = 3  # Number of digits in each part

    parts1=split_number_np_decimal(u,m)
    print("parts1 = ",parts1)
    parts2=split_number_np_decimal(v,m)
    # Perform the division using Burnikel-Ziegler method
    quotient, remainder = burnikel_ziegler_divide(u, v, m)
    print("Quotient:", quotient)
    print("Remainder:", remainder)
    
