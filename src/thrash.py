import numpy as np

from math import *



def karatsuba_division(dividend, divisor):
    quotient = []
    while len(dividend) >= len(divisor):
        leading_term_multiplier = dividend[0] / divisor[0]
        quotient_term_degree = len(dividend) - len(divisor)
        quotient_term = [0] * quotient_term_degree + [leading_term_multiplier]

        quotient.append(leading_term_multiplier)
        product = karatsuba(divisor, quotient_term)
        dividend = [dividend[i] - (product[i] if i < len(product) else 0) for i in range(len(dividend))]

        # Remove leading zeros
        while len(dividend) > 0 and dividend[0] == 0:
            dividend.pop(0)

    remainder = dividend
    return quotient, remainder



def deg(a):
    return len(a)

def poly_divide(u,v):
    m = deg(u)
    n = deg(v)
    print("m = ",m,"n = ",n)
    if (m > n):
        p = m -n
    else:
        p = n - m
    q = np.zeros(p).tolist()
    print("q = ",q)
    for k in range(p,0,-1):
        q[k]=u[n+k]/v[n]
        for j in range(n+k-1,k,-1):
            u[j]=u[j]-q[k]*v[j-k]
    r = u[0:n]    

    return (q,r)




def sum_poly(a,b):
    if deg(a)>deg(b):
        return [a[i]+b[i] for i in range(len(b))]+a[len(b):]
    else:
        return [a[i]+b[i] for i in range(len(a))]+b[len(a):]




def m10(x,m,add_to_front=True):
    y=[]
    for c in x:
        y.append(c)
    for i in range(m):
        if (add_to_front):
            y.insert(0,0)
        else:
            y.append(0)
    return y
#x and y have to be same length
def add(x, y):
    s = len(x)
    t = len(y)
    if (s > t):
        d = s - t
        y = m10(y,d)
    else:
        d = t - s
        x = m10(x,d)
    
    z=[]
    carry = 0
    tot=0
    n = max(s,t)
    for i in range(n-1,-1,-1):
        tot=x[i]+y[i]+carry
        if (tot >=10):
            carry = 1
            z.insert(0,tot % 10 )
            if (n == 1):
                z.insert(0,1)
            if i == 0:
                z.insert(0,carry)    
        else:
            carry = 0
            z.insert(0,tot)
    return z
def swap(x,y):
    temp = x
    x = y
    y = temp
    return [x,y]

def smul(x,y):
    n = len(x)
    m = len(y)
    carry = 0
    shift = 0
    vecs =[]
    for i in range(m-1,-1,-1):
        z=[]
        carry = 0
        for j in range(n-1,-1,-1):
  
            t = x[j]*y[i]+carry
            carry = 0
            if i == 1:
                print("i = ",i,"j = ",j,"t = ",t,"carry=",carry)
                                       
            if (t>=10):
                carry,s=divmod(t,10)
                if (j == 0):
                    z.insert(0,s)
                    z.insert(0,carry)
                    
                else:    
                    z.insert(0,s)
            else:
                z.insert(0,t)        
        z = m10(z,shift,False)
        shift = shift+1
        vecs.append(z)
    a=[]
    for i in range(len(vecs)):    
        print(vecs[i])
        a = add(a,vecs[i])
        #print(a)
    
    return a



def karatsuba(x,y):

    n = len(x)
    if n <= 3:
        return smul(x, y)
    else:
        m = n // 2
        A = np.zeros(2*n + 1).tolist()
        x0, x1 = x[0 : m], x[m : n]
        y0, y1 = y[0 : m], y[m : n]
        u = add(x1, x0)
        v = add(y1, y0)
        p0 = karatsuba(x0, y0)
        p1 = karatsuba(x1, y1)
        p2 = karatsuba(u, v)
        A[0 : 2*m] = p0
        A[2*m : 2*n] = p1
        A[m : 2*n+1] = add(A[m : 2*n], p2)
        A[m : 2*n+1] = sub(A[m : 2*n+1], p0)
        A[m : 2*n+1] = sub(A[m : 2*n+1], p1)
        return A[0 : 2*n]


def barrett_reduction(a, n):
    """
    Perform Barrett Reduction on 'a' modulo 'n'.
    
    :param a: The integer to reduce.
    :param n: The modulus.
    :return: The result of a mod n using Barrett Reduction.
    """
    # k is the number of bits in n
    k = n.bit_length()
    
    # Precompute the Barrett constant mu = floor(2^(2k) / n)
    mu = (1 << (2 * k)) // n
    
    # Barrett Reduction
    q = ((a * mu) >> (2 * k)) * n
    r = a - q
    if r >= n:
        r -= n
    return r

def deg(a):
    return len(a)-1


def euclidean_division(a,b):
   anp = np.array(a)
   bnp = np.array(b)
   q = 0
   r = anp
   d = len(b)-1
   c=b[0]
   while deg(r) >= d:
       s=r[0]/c
       q = q+s
       r = r -np.multiply(s,b)
   return (q,r)



def newton_raphson_division(dividend, divisor, precision):
    # Initial guess for the reciprocal of the divisor
    #reciprocal_estimate = 1.0 / divisor

    reciprocal_estimate = 1e-10

    # Refine the estimate using Newton-Raphson iteration
    for _ in range(precision):
        reciprocal_estimate = reciprocal_estimate * (2 - divisor * reciprocal_estimate)
    
    # Calculate the quotient using the reciprocal of the divisor
    quotient = dividend * reciprocal_estimate
    
    return quotient

if __name__=="__main__":
    dividend = [1, 0, -12, 0]  # x^3 - 12
    divisor = [1, -3]  # x - 3
    quotient, remainder = karatsuba_division(dividend, divisor)

    print("Quotient:", quotient)
    print("Remainder:", remainder)
    
  


   

