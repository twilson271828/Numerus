def karatsubaMultiply(x,y):
	# Defining base case for recursive calls
	if len(str(x)) == 1 or len(str(y)) == 1:
		return x * y
	else:
		n = max(len(str(x)),len(str(y)))
		n_by_2 = n // 2

		# Step1: Rewriting x and y in terms of a,b,c,d
		a = x // 10**(n_by_2)  #right shift by n_by_2
		b = x % 10**(n_by_2) #remainder after right shift by n_by_2
		c = y // 10**(n_by_2) #right shift by n_by_2
		d = y % 10**(n_by_2) #remainder after right shift by n_by_2

		# Step2: Computing ac
		ac = karatsubaMultiply(a,c)

		# Step3: Computing bd
		bd = karatsubaMultiply(b,d)

		# Step4: Computing ad+bc
		ad_bc = karatsubaMultiply(a+b,c+d) - ac - bd

		# Step5: Calculating finalValue
		finalAns = ac * 10**(2*n_by_2) + ad_bc * 10**(n_by_2) + bd

		return finalAns

def karatsuba_copilot(x, y):
    
	if x < 10 or y < 10:
		#print("x = ",x)
		#print("y = ",y)
		return x * y


	n = max(len(str(x)), len(str(y)))
	m = n // 2
	

	x_high, x_low = divmod(x, 10**m)
	y_high, y_low = divmod(y, 10**m)
	
	#print("**********************************")
	#print("x: ",x)
	#print("y: ",y)
	#print("x_high: ",x_high)
	#print("x_low: ",x_low)
	#print("y_high: ",y_high)
	#print("y_low: ",y_low)
	#print("k = ",n)
	#print("k2 = ",m)

	#print("**********************************")
	#print("k: ",n)
	#print("k2:  ",m)
	#print("x: ",x)
	#print("y: ",y)
	#print("x_high: ",x_high)
	#print("x_low: ",x_low)
	#print("y_high: ",y_high)
	#print("y_low: ",y_low)
	#print("****************************************")
	
	z0 = karatsuba_copilot(x_low, y_low)
	c1 = x_low+x_high
	c2 = y_low+y_high
	z1 = karatsuba_copilot((x_low + x_high), (y_low + y_high))
	z2 = karatsuba_copilot(x_high, y_high)

	print("****************************************")
	#print("k2 =  ",m)
	#print("x = ",x)
	#print("y = ",y)
	print("x_high =",x_high)
	print("x_low =",x_low)
	print("y_high =",y_high)
	print("y_low =",y_low)
	print ("z0 =",z0)
	print("c1 =",c1)
	print("c2 =",c2)
	print ("z1 =",z1)
	print("z2 =",z2)
	print ("z2 * 10**(2*m) =",z2 * 10**(2*m))
	print ("(z1 - z2 - z0) * 10**m =", (z1 - z2 - z0) * 10**m)
	print("****************************************")

	return (z2 * 10**(2*m)) + ((z1 - z2 - z0) * 10**m) + z0
	




def main():
	#x=230
	#y=3204
	x = 90890809809
	y=48909889098

	z=karatsuba_copilot(x,y)
	print(z)
	print(x*y)
	
if __name__ == "__main__":
	main()
	
