from objects import * 
from sympy import *
from random import randint  
import sys 


if len(sys.argv) < 4:
	raise Exception("Input format: <Protocol> <Object> <FiniteField> || <mode flag>")

__Protocol   	= sys.argv[1]
__Object 		= sys.argv[2]
__FiniteField 	= sys.argv[3] 

# Input option
if len(sys.argv) == 5:
	__Input_Flag = sys.argv[4]
else: 
	__Input_Flag = "NA"

x,y = symbols('x,y')
__Sample_EC_Curves 	= [y**2 - x**3 - x + 1]
__Sample_HEC_Curves = [y**2 - x**5 + 3*x**4 + 2*x**3 - 6*x**2 + 3*x - 1]

if __Object == "EC":
	__Group = EllipticCurve([__Sample_EC_Curves[randint(0,len(__Sample_EC_Curves)-1)]],FiniteField(int(__FiniteField)))
if __Object == "HEC":
	__Group = HyperEllipticCurve([__Sample_HEC_Curves[randint(0,len(__Sample_HEC_Curves)-1)]],FiniteField(int(__FiniteField)))

if __Protocol == "ECDH":

	if __Input_Flag == "-i":
		private_key 	= int(input("Enter private key: "))
		curve 			= input("Enter elliptic Curve: ")
		__Group.set_polynomials([curve])
		public_point 	= tuple(input("Enter public point: "))
		public_key 		= __Group.scalar_mult(private_key,public_point)

		print "Public key: P=%s " 	% str(public_key)
		print "Curve: 0 = %s " 			% curve
		print "Public Point: G=%s"	% str(public_point)


	elif __Input_Flag == "-t":
		curve 			= input("Enter elliptic Curve: ")
		__Group.set_polynomials([curve])

		private_key_A 	= int(input("Enter private key of A: "))
		public_key_B 	= tuple(input("Enter public key of B: "))

		key_1 = __Group.scalar_mult(private_key_A,public_key_B )

		print "Key is: %s" % str(key_1)

		private_key_B 	= int(input("Enter private key of B: "))
		public_key_A 	= tuple(input("Enter public key of A: "))

		key_2 = __Group.scalar_mult(private_key_B,public_key_A)

		print "Key is: %s" % str(key_2)

		if key_1 == key_2:
			print "Alice and Bob share a common key: %s " % str(key_1)
		else:
			print "Commmon key not shared."


	else:
		private_key = int(input("Enter private key: "))
		curve 			= __Sample_EC_Curves[0]
		public_point 	= __Group.random_point()
		public_key 		= __Group.scalar_mult(private_key,public_point)

		print "Public key: P=%s " 	% str(public_key)
		print "Curve: 0 = %s " 			% curve
		print "Public Point: G=%s"	% str(public_point)







 






