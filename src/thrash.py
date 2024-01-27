import numpy as np

def multiply_poly_fft(p1, p2):
    size = len(p1) + len(p2) - 1
    p1_fft = np.fft.fft(p1, size)
    p2_fft = np.fft.fft(p2, size)
    result_fft = p1_fft * p2_fft
    result = np.fft.ifft(result_fft)
    return np.round(result).astype(int)

def deg(a):
    return len(a)-1

def sum_poly(a,b):
    if deg(a)>deg(b):
        return [a[i]+b[i] for i in range(len(b))]+a[len(b):]
    else:
        return [a[i]+b[i] for i in range(len(a))]+b[len(a):]

def karatsuba_poly_mult(p1, p2):
    n = max(len(p1), len(p2))
    p1 = p1 + [0] * (n - len(p1))  # Pad with zeros to make the lengths equal
    p2 = p2 + [0] * (n - len(p2))  # Pad with zeros to make the lengths equal

    if n == 1:
        return [p1[0] * p2[0]]

    m = n // 2
    p1_low = p1[:m]
    p1_high = p1[m:]
    p2_low = p2[:m]
    p2_high = p2[m:]

    z0 = karatsuba_poly_mult(p1_low, p2_low)
    z1 = karatsuba_poly_mult([x + y for x, y in zip(p1_low, p1_high)], [x + y for x, y in zip(p2_low, p2_high)])
    z2 = karatsuba_poly_mult(p1_high, p2_high)

    return z0 + [z1[i] - z0[i] - z2[i] for i in range(len(z1))] + z2




def euclidean_division(a,b):
    q=0
    r=a
    d=len(b)-1
    c=b[0]
    while deg(r)>=d:
        s=(r[0]/c)*(x**(deg(r)-d))
        q=sum_poly(q,s)
        r=r - s*b




if __name__=="__main__":

    # Test the function
    #print(euclidean_division([1, 3, 2], [5, 3])) 
    print(karatsuba_poly_mult([1, 2, 3], [4, 5, 6]))  # Output: [4, 13, 28, 27, 18] 


   

