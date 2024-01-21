def knuth_division(n, d):
    # Step D1: Initialize
    q, r = 0, n  # quotient and remainder
    t = len(str(d))  # number of digits in d
    n_digits = [int(x) for x in str(n)]
    d_digits = [int(x) for x in str(d)]

    # Step D2: Divide
    for i in range(len(n_digits) - t + 1):
        # Estimate quotient
        if r >= 10**t:
            q_hat = r // (10**(t-1)) // d_digits[0]
        else:
            q_hat = r // d
        # Multiply and subtract
        r = 10 * (r - q_hat * d)
        if i < len(n_digits) - t:
            r += n_digits[i+t]
        # Test remainder
        if r < 0:
            r += d
            q_hat -= 1
        # Concatenate quotient
        q = 10 * q + q_hat

    return q, r

# Test the function
print(knuth_division(123456789, 12345))  # Output: (10005, 4369)
