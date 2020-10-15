#I'm too lazy to make this work for both split and inert primes, so this is just for inert primes. I will give you slopes of U2 at level 8*p_11 where p_11|11. Make sure you have the SquareRoot5 file in the same place. If the code doesnt work blame magma and/or sage, if the code works blame me.

#This creates the matrices in magma, imports them into sage from which you can then compute the char polys and the slopes.


magma.eval('load "SquareRoot5"')

from sage.geometry.newton_polygon import NewtonPolygon



#Here the input 'pres' is the precisision at which you want your matrices to come out with. Use something between 1000 to 6000 if you want to get all the slopes

def OCPrecomp(pres): 
	L=magma('SageComputeOverconvergentUp(' + str(pres) + ')')
	def Up(k1,k2,n):
			LIT=L(k1,k2,n).sage()
		 	A=Zp(2,pres)
	        	S.<x>=A[]
	        	B.<t>=A.ext(x^2+x+1)
	        	P.<y>=B[]
	        	LI=[B(x) for x in LIT]
		    	M=Matrix(B,16*n,16*n,LI)  
		    	return M
 	return Up


#The output is a function Up(k_1,k_2,n) whose input is k_1,k_2 and n. It will then return a 16*n x 16*n matrix of Up acting in weight k_1,k_2 (make sure the parity of the weights is the same). Here 16 is what I call the class number of the quaternion algebra at this level, i.e. the size of the shim var Y_D(Level).


#Once you have the matrices, then you can feed them into the Slopes function to get its slopes. If start seeing slopes like 899/57 it means that the precision you are working with isnt hight enough.

def Slopes(M):
		 B=M.base_ring()
		 P.<y>=B[]	
		 ch=characteristic_polynomial(M)
		 O=[(i,valuation(P(ch)[i])) for i in range(M.ncols()+1)]
	 	 NP=NewtonPolygon(O)
		 E=NP.slopes()
		 E.reverse()
		 EE=[-1*x for x in E]		
		 return uniq([(x,EE.count(x)) for x in EE])
