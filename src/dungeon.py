def split_number(x, m):
    """
    Splits a large number `x` into k parts of size `m` where k = ceil(len(x)/m).
    Returns the list of parts.
    """
    parts = []
    while x > 0:
        parts.insert(0, x & ((1 << m) - 1))  # Extract the last m bits
        x >>= m
    return parts

def burnikel_ziegler_divide(u, v):
    """
    Burnikel-Ziegler division algorithm.
    Divides `u` by `v` where u and v are large integers.
    Returns the quotient and remainder.
    """
    if v == 0:
        raise ZeroDivisionError("Division by zero is undefined.")
    
    if u < v:
        return 0, u  # quotient is 0 and remainder is u if u < v

    # Length of u and v
    n = v.bit_length()
    m = n // 2

    # Splitting the numbers into parts of m bits
    u_parts = split_number(u, m)
    v_parts = split_number(v, m)

    # Burnikel-Ziegler division steps:
    # Step 1: Perform divide-and-conquer division on parts
    quotient = 0
    remainder = 0

    for part in u_parts:
        remainder = (remainder << m) + part  # Shift remainder and add next part
        quotient_part = remainder // v
        remainder = remainder % v
        quotient = (quotient << m) + quotient_part  # Build up the quotient

    return quotient, remainder

# Example usage:
#u = 2**512 + 2**256 + 1  # A large integer
#v = 2**256 + 1            # Another large integer

u = 2463464354353453455
v = 325234223
x = 358490345905438953424224266424434234432
m=10
parts=split_number(x,m)
print(parts)
print(len(parts))

#quotient, remainder = burnikel_ziegler_divide(u, v)
#print("Quotient:", quotient)
#print("Remainder:", remainder)
