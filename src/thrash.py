def karatsubaMultiply(x,y):
	# Defining base case for recursive calls
	if len(str(x)) == 1 or len(str(y)) == 1:
		return x * y
	else:
		n = max(len(str(x)),len(str(y)))
		n_by_2 = n // 2

		# Step1: Rewriting x and y in terms of a,b,c,d
		a = x // 10**(n_by_2)
		b = x % 10**(n_by_2)
		c = y // 10**(n_by_2)
		d = y % 10**(n_by_2)

		# Step2: Computing ac
		ac = karatsubaMultiply(a,c)

		# Step3: Computing bd
		bd = karatsubaMultiply(b,d)

		# Step4: Computing ad+bc
		ad_bc = karatsubaMultiply(a+b,c+d) - ac - bd

		# Step5: Calculating finalValue
		finalAns = ac * 10**(2*n_by_2) + ad_bc * 10**(n_by_2) + bd

		return finalAns
	
def main():
	x=2322533454334543534534534
	y=324423443423555534534534534553455345355453553454
	print(karatsubaMultiply(x,y))
	print(x*y)
	
if __name__ == "__main__":
	main()
	
