#hashed cells corrected (values are active if and only if at least diff or value is active, when both are inactive, value is not active)
#corrected the Kin and Kout
#improved version of key recovery, guessing the internal state in effective cases
from gurobipy import *
DistRound=11
PreRound=5
PostRound=6
marker = 0
filename_model = "CRAFT_"+str(PreRound)+"_"+str(DistRound)+"_"+str(PostRound)+".lp"
fileobj = open(filename_model, "w")
fileobj.close()
class CRAFT:
	MC = [[1,0,1,1],[0,1,0,1],[0,0,1,0],[0,0,0,1]]
	Inv_MC = MC #inverse of MC
	def column(A, j):
	    return [A[j], A[j+4], A[j+8], A[j+12]]

	def PN(A):
	    return [A[15], A[12], A[13], A[14],\
		    A[10], A[9], A[8], A[11],\
		    A[6],A[5],A[4], A[7],\
		    A[1],A[2],A[3],A[0]]

		    
class Truncated:
	S_T = [(1, 2, 1, 3, -1, -2, 0, 0, -5, -2, 0),
 (-1, 0, -1, -3, 3, 2, 0, 0, 5, 2, 0),
 (0, 1, 0, 1, -2, -1, 0, 0, -2, -2, 2),
 (-1, -2, -1, 0, 2, 2, 0, 0, 3, 2, 0), 
(2, 1, 2, 3, -2, -1, 0, 0, -6, -2, 0),
 (0, 0, 0, -1, 0, 1, 0, 0, 1, 1, 0), 
(0, 0, 0, 0, -1, 0, 0, 0, -1, 0, 1),
 (-1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0),
(1, 0, 1, 1, -1, -2, 0, 0, -3, -2, 2)]


	def genEncSubjectionAtRound(r):

	    eqn = []

	    inX = BasicTools.inMixAtRound(r)
	    outY = BasicTools.inMixAtRound(r+1)

	    eqn = eqn + Truncated.mixColumnAndShuffle(inX,outY)

	    return eqn
	    
	def genEncSubjectionAtRoundPre(r):

	    eqn = []

	    inX = ['x_InPN_' + str(r) + '_' + str(i) for i in range(0,16)]
	    outY = ['x_InMC_' + str(r) + '_' + str(i) for i in range(0,16)]

	    eqn = eqn + Truncated.mixColumn(inX,outY)

	    return eqn


	def genEncSubjection(totalRound,PreRound,PostRound):
		eqn = []

		for i in range(PreRound):
			eqn = eqn + Truncated.genEncSubjectionAtRoundPre(i)
		for i in range(PreRound,totalRound + PreRound + PostRound):
			eqn = eqn + Truncated.genEncSubjectionAtRound(i)

		return eqn


	
	def mixColumnSubjection(inX,outY,ineq):
	    assert(len(inX) == len(outY) and len(inX) == 4)

	    global marker
	    eqn = []
	    T = ineq

	    for t in T:
	    	eqn.append((str(t[0]) + ' ' + inX[0] + ' + ' + str(t[1]) + ' ' + inX[1] + ' + ' + str(t[2]) + ' ' + inX[2] + ' + ' + str(t[3]) + ' ' + inX[3] + ' + ' + str(t[4]) + ' ' + outY[0] + ' + ' + str(t[5]) + ' ' + outY[1] + ' + ' + str(t[6]) + ' ' + outY[2] + ' + ' + str(t[7]) + ' ' + outY[3] + ' + ' +  str(t[8]) + ' ' + 'p1' + '_' + str(marker) + ' + ' + str(t[9]) + ' ' + 'p0' +  '_' + str(marker) + ' >= ' + str(-t[10])).replace('+ -','- '))
	    eqn.append(inX[2] + ' - ' + outY[2] + ' = 0')
	    eqn.append(inX[3] + ' - ' + outY[3] + ' = 0')

	    marker = marker +1

	    return eqn
	    
	def mixColumnAndShuffle(inX,outY):

	    eqn = []

	    #outY = [outY[15],outY[10],outY[9],outY[4],outY[3],outY[6],outY[5],outY[8],outY[7],outY[2],outY[1],outY[12],outY[11],outY[14],outY[13],outY[0]]
	    outY = [outY[15],outY[10],outY[6],outY[1],outY[12],outY[9],outY[5],outY[2],outY[13],outY[8],outY[4],outY[3],outY[14],outY[11],outY[7],outY[0]]
	    
	    inX1 = [inX[0],inX[4],inX[8],inX[12]]
	    inX2 = [inX[1],inX[5],inX[9],inX[13]]
	    inX3 = [inX[2],inX[6],inX[10],inX[14]]
	    inX4 = [inX[3],inX[7],inX[11],inX[15]]
	    outY1 = outY[0:4]
	    outY2 = outY[4:8]
	    outY3 = outY[8:12]
	    outY4 = outY[12:16]

	    eqn = eqn + Truncated.mixColumnSubjection(inX1,outY1,Truncated.S_T)
	    eqn = eqn + Truncated.mixColumnSubjection(inX2,outY2,Truncated.S_T)
	    eqn = eqn + Truncated.mixColumnSubjection(inX3,outY3,Truncated.S_T)
	    eqn = eqn + Truncated.mixColumnSubjection(inX4,outY4,Truncated.S_T)
	    
	    return eqn
	def mixColumn(inX,outY):

	    eqn = []

	    #outY = [outY[15],outY[10],outY[9],outY[4],outY[3],outY[6],outY[5],outY[8],outY[7],outY[2],outY[1],outY[12],outY[11],outY[14],outY[13],outY[0]]
	    #outY = [outY[15],outY[10],outY[6],outY[1],outY[12],outY[9],outY[5],outY[2],outY[13],outY[8],outY[4],outY[3],outY[14],outY[11],outY[7],outY[0]]
	    
	    inX1 = [inX[0],inX[4],inX[8],inX[12]]
	    inX2 = [inX[1],inX[5],inX[9],inX[13]]
	    inX3 = [inX[2],inX[6],inX[10],inX[14]]
	    inX4 = [inX[3],inX[7],inX[11],inX[15]]
	    outY1 = [outY[0],outY[4],outY[8],outY[12]]
	    outY2 = [outY[1],outY[5],outY[9],outY[13]]
	    outY3 = [outY[2],outY[6],outY[10],outY[14]]
	    outY4 = [outY[3],outY[7],outY[11],outY[15]]

	    eqn = eqn + Truncated.mixColumnSubjection(inX1,outY1,Truncated.S_T)
	    eqn = eqn + Truncated.mixColumnSubjection(inX2,outY2,Truncated.S_T)
	    eqn = eqn + Truncated.mixColumnSubjection(inX3,outY3,Truncated.S_T)
	    eqn = eqn + Truncated.mixColumnSubjection(inX4,outY4,Truncated.S_T)
	    
	    return eqn
	    
class BasicTools:
	def transpose(M):
		m = len(M)
		n = len(M[0])

		Mt = []
		for i in range(0, n):
			row = [M[k][i] for k in range(0, m)]
			Mt.append(row)

		return Mt
	def inMixAtRound(r):
    		#assert(r >= 1)
    		return ['x_InMC_' + str(r) + '_' + str(i) for i in range(0,16)]
	def VarGen(s,n):
	    return [str(s) + '_' + str(n) + '_' + str(i) for i in range(0,16)]
	def plusTerm(in_vars):
		"""
		>>> BasicTools.plusTerm(['x','y','z'])
		'x + y + z'
		>>> BasicTools.plusTerm(['x','y'])
		'x + y'
		>>> BasicTools.plusTerm(['x','y','z','a','b'])
		'x + y + z + a + b'
		>>>
		"""
		t = ''
		for v in in_vars:
		    t = t + v + ' + '

		return t[0:-3]

	def MinusTerm(in_vars):
		"""
		>>> BasicTools.plusTerm(['x','y','z'])
		'x + y + z'
		>>> BasicTools.plusTerm(['x','y'])
		'x + y'
		>>> BasicTools.plusTerm(['x','y','z','a','b'])
		'x + y + z + a + b'
		>>>
		"""
		t = ''
		for v in in_vars:
		    t = t + v + ' - '

		return t[0:-3]	    
	def equalConstraints(x, y):
		assert len(x) == len(y)
		c = []
		for i in range(0, len(x)):
	    		c = c + [x[i] + ' - ' + y[i] + ' = 0']
		return c
		
	def greaterConstraints(x, y):
		assert len(x) == len(y)
		c = []
		for i in range(0, len(x)):
	    		c = c + [x[i] + ' - ' + y[i] + ' >= 0']
		return c
		
	def greaterConstraints_xyz(x, y, z):
		assert len(x) == len(y)
		c = []
		for i in range(0, len(x)):
	    		c = c + [x[i] + ' + ' + y[i] + ' - ' + z[i] + ' >= 0']
		return c			    
		
	def transpose(M):
		"""
		Transpose the matrix M
		>>> M = [[1,0,1,1],[1,0,0,0],[0,1,1,0],[1,0,1,0]]
		>>>
		>>> BasicTools.transpose(M)
		[[1, 1, 0, 1], [0, 0, 1, 0], [1, 0, 1, 1], [1, 0, 0, 0]]
		>>>
		>>>
		"""
		m = len(M)
		n = len(M[0])

		Mt = []
		for i in range(0, n):
		    row = [M[k][i] for k in range(0, m)]
		    Mt.append(row)

		return Mt		
	def getVariables(C):
		V = set([])
		#for s in C :
		#print(s)
		temp = C.strip()
		#print(temp)
		        
		temp = temp.replace('+', ' ')
		temp = temp.replace('-', ' ')
		temp = temp.replace('>=', ' ')
		temp = temp.replace('<=', ' ')
		temp = temp.replace('=', ' ')
		        
		#print(temp)
		temp = temp.split()
		#print(temp)
		for v in temp :
		        if not v.isdecimal():
		                V.add(v)
		#print(V)
		return V

class Extension:	
	print("Extension")	
	def ForwardDiff_LinearLayer(M, V_in, V_out):
		"""
		>>> M = [[1,0,1,1],[1,0,0,0],[0,1,1,0],[1,0,1,0]]
		>>> a = ['a0', 'a1', 'a2', 'a3']
		>>> b = ['b0', 'b1', 'b2', 'b3']
		>>>
		>>> for c in MITMConstraints.ForwardDiff_LinearLayer(M, a, b): print(c)
		...
		3 b0 -  a0 - a2 - a3 >= 0
		a0 + a2 + a3 - b0 >= 0
		1 b1 -  a0 >= 0
		a0 - b1 >= 0
		2 b2 -  a1 - a2 >= 0
		a1 + a2 - b2 >= 0
		2 b3 -  a0 - a2 >= 0
		a0 + a2 - b3 >= 0
		>>>
		"""
		assert len(M[0]) == len(V_in), "The input is not compatible with the matrix"
		assert len(M) == len(V_out), "The output is not compatible with the matrix"


		m = len(M)
		n = len(M[0])

		constr = []
		for i in range(0, m):
		    s = sum(M[i]) # the number of 1s in row i
		    terms = [V_in[j] for j in range(0, n) if M[i][j] == 1]
		    constr = constr + [str(s) + ' ' + V_out[i] + ' - ' + ' ' + BasicTools.MinusTerm(terms) + ' >= 0']
		    constr = constr + [BasicTools.plusTerm(terms) + ' - ' + V_out[i] + ' >= 0']

		return constr
	def LinearLayer(M, V_in, V_out):
		"""
		>>> M = [[1,0,1,1],[1,0,0,0],[0,1,1,0],[1,0,1,0]]
		>>> a = ['a0', 'a1', 'a2', 'a3']
		>>> b = ['b0', 'b1', 'b2', 'b3']
		>>>
		>>> for c in MITMConstraints.ForwardDiff_LinearLayer(M, a, b): print(c)
		...
		3 b0 -  a0 - a2 - a3 >= 0
		a0 + a2 + a3 - b0 >= 0
		1 b1 -  a0 >= 0
		a0 - b1 >= 0
		2 b2 -  a1 - a2 >= 0
		a1 + a2 - b2 >= 0
		2 b3 -  a0 - a2 >= 0
		a0 + a2 - b3 >= 0
		>>>
		"""
		assert len(M[0]) == len(V_in), "The input is not compatible with the matrix"
		assert len(M) == len(V_out), "The output is not compatible with the matrix"


		m = len(M)
		n = len(M[0])

		constr = []
		for i in range(0, m):
		    s = sum(M[i]) # the number of 1s in row i
		    terms = [V_in[j] for j in range(0, n) if M[i][j] == 1]
		    constr = constr + [str(s) + ' ' + V_out[i] + ' - ' + ' ' + BasicTools.MinusTerm(terms) + ' >= 0']
		    constr = constr + [BasicTools.plusTerm(terms) + ' - ' + V_out[i] + ' >= 0']

		return constr	

	def McRelImprovementForward(inX,W):
		#T = [("x3+y2'"),("x1+y0'"),("x3+y0'"),("x0+y3'"),("x2+y1'"),("x2'+y1"),("x0'+y3"),("y2'+y3"),("x1'+x3'+y0+y3'"),("y1+y2'"),("x3'+y1'+y2+y3'"),("y0'+y3")]
		T=[("x1+y1'"),("x2'+y2'"),("x0+y0'"),("x3'+y3'"),("x2+x0'+y2"),("x3+x1'+y1"),("y2'+y0"),("y3+y1'"),("x3+x0'+y0"),("y3+y2+y0'"),("y3'+y1+y0"),("x3+y3+y0'")]
		#T=[("x2+y3'"),("x0+y1'"),("x2+y1'"),("x1+y2'"),("x3+y0'"),("x3'+y0"),("x1'+y2"),("y2+y3'"),("x0'+x2'+y1+y2'"),("y0+y3'"),("x2'+y0'+y2'+y3"),("y1'+y2")]
		eqn = []

		for t in T:
			cnt=0
			eq=""
			if ("x0'" in t):
				cnt+=1
				eq=eq+("- " + inX[3])
				#eq=eq+("1 - " + inX[0])
			elif ("x0" in t):
				eq=eq+(inX[3])
			if ("x1'" in t):
				cnt+=1
				eq=eq+(" - " + inX[2])
				#eq=eq+(" + 1 - " + inX[1])
			elif ("x1" in t):
				eq=eq+(" + " + inX[2])
			if ("x2'" in t):
				cnt+=1
				eq=eq+(" - " + inX[1])
				#eq=eq+(" + 1 - " + inX[2])
			elif ("x2" in t):
				eq=eq+(" + " + inX[1])
			if ("x3'" in t):
				cnt+=1
				eq=eq+(" - " + inX[0])
				#eq=eq+(" + 1 - " + inX[3])
			elif ("x3" in t):
				eq=eq+(" + " + inX[0])
			if ("y0'" in t):
				cnt+=1
				eq=eq+(" - " + W[3])
			elif ("y0" in t):
				eq=eq+(" + " + W[3])
			if ("y1'" in t):
				cnt+=1
				eq=eq+(" - " + W[2])
			elif ("y1" in t):
				eq=eq+(" + " + W[2])
			if ("y2'" in t):
				cnt+=1
				eq=eq+(" - " + W[1])
			elif ("y2" in t):
				eq=eq+(" + " + W[1])
			if ("y3'" in t):
				cnt+=1
				eq=eq+(" - " + W[0])
			elif ("y3" in t):
				eq=eq+(" + " + W[0])
			eq=eq + " >= " + str(-cnt+1)
			eqn.append(eq)
		return eqn
		
	def McRelImprovementBackward(inX,W):
		#T = [("x3+y2'"),("x1+y0'"),("x3+y0'"),("x0+y3'"),("x2+y1'"),("x2'+y1"),("x0'+y3"),("y2'+y3"),("x1'+x3'+y0+y3'"),("y1+y2'"),("x3'+y1'+y2+y3'"),("y0'+y3")]
		T=[("x1+y1'"),("x2'+y2'"),("x0+y0'"),("x3'+y3'"),("x2+x0'+y2"),("x3+x1'+y1"),("y2'+y0"),("y3+y1'"),("x3+x0'+y0"),("y3+y2+y0'"),("y3'+y1+y0"),("x3+y3+y0'")]
		#T=[("x2+y3'"),("x0+y1'"),("x2+y1'"),("x1+y2'"),("x3+y0'"),("x3'+y0"),("x1'+y2"),("y2+y3'"),("x0'+x2'+y1+y2'"),("y0+y3'"),("x2'+y0'+y2'+y3"),("y1'+y2")]
		eqn = []

		for t in T:
			cnt=0
			eq=""
			if ("x0'" in t):
				cnt+=1
				eq=eq+("- " + inX[3])
				#eq=eq+("1 - " + inX[0])
			elif ("x0" in t):
				eq=eq+(inX[3])
			if ("x1'" in t):
				cnt+=1
				eq=eq+(" - " + inX[2])
				#eq=eq+(" + 1 - " + inX[1])
			elif ("x1" in t):
				eq=eq+(" + " + inX[2])
			if ("x2'" in t):
				cnt+=1
				eq=eq+(" - " + inX[1])
				#eq=eq+(" + 1 - " + inX[2])
			elif ("x2" in t):
				eq=eq+(" + " + inX[1])
			if ("x3'" in t):
				cnt+=1
				eq=eq+(" - " + inX[0])
				#eq=eq+(" + 1 - " + inX[3])
			elif ("x3" in t):
				eq=eq+(" + " + inX[0])
			if ("y0'" in t):
				cnt+=1
				eq=eq+(" - " + W[3])
			elif ("y0" in t):
				eq=eq+(" + " + W[3])
			if ("y1'" in t):
				cnt+=1
				eq=eq+(" - " + W[2])
			elif ("y1" in t):
				eq=eq+(" + " + W[2])
			if ("y2'" in t):
				cnt+=1
				eq=eq+(" - " + W[1])
			elif ("y2" in t):
				eq=eq+(" + " + W[1])
			if ("y3'" in t):
				cnt+=1
				eq=eq+(" - " + W[0])
			elif ("y3" in t):
				eq=eq+(" + " + W[0])
			eq=eq + " >= " + str(-cnt+1)
			eqn.append(eq)
		return eqn
			
	def Determination_Decision(M, V_in, V_out , t_variable):
		assert len(M[0]) == len(V_in)
		assert len(M) == len(V_out)
		m = len(M)
		n = len(M[0])
		constr = []
		# constr2=[]
		# constr3=[]
		# constr4=[]
		for j in range(0, 4):#This Forloop Indicates That How Yi and Ti(variable) Related (Yi >= Ti).Yi is output of mixcolumn matrix, and for each active nibble we consider T variable
			constr=constr+[V_out[j] + ' - ' + t_variable[j] + ' >= 0' ]
		if True:
			for i in range(0, m):#This Forloop Shows That If Yi are active and Ti variables is non-active
				s = sum(M[i]) # the number of 1s in row i
				terms1=[V_in[j] for j in range(0, n) if M[i][j] == 1]
				constr = constr + [BasicTools.plusTerm(terms1) + ' - ' + str(s) + ' ' + V_out[i] + ' + ' + str(s) + ' ' + t_variable[i] + ' >= 0']
		if True:
			# for j in range(0, m):#This Forloop Is used to Prevent ai(input of mixcolumn matrix) from acting Unnecessarily
			#      terms2=[V_out[i] for i in range(0, n) if M[i][j] == 1]
			#      constr3 = constr3 + [Tools.plusTerm(terms2) + ' - ' + ' ' + V_in[j] + ' >= 0'] This Forloop Is not necessary because the Next One covers All Contraint That we Need In Our Model

			for j in range(0, m):#This Forloop shows that if Ti variable are now active and we want to describe the constraints
				s = sum(M[k][j] for k in range(0,len(M))) # the number of 1s in column j
				terms1=[t_variable[i] for i in range(0, n) if M[i][j] == 1]
				terms2 = [ V_out[i] for i in range(0, n) if M[i][j] == 1]
				constr = constr + [BasicTools.plusTerm(terms1) + ' + ' + ' ' + V_in[j] + ' ' + ' - ' +  (BasicTools.MinusTerm(terms2)) + ' <= 0']#original       
		return constr
        	
	def BackwardDet_LinearLayer(M, V_in, V_out):
		"""
		>>> M = [[1,0,1,1],[1,0,0,0],[0,1,1,0],[1,0,1,0]]
		>>> a = ['a0', 'a1', 'a2', 'a3']
		>>> b = ['b0', 'b1', 'b2', 'b3']
		>>> MITMConstraints.BackwardDet_LinearLayer(M, a, b)
		['3 a0 -  b0 - b1 - b3 >= 0',
		 'b0 + b1 + b3 - a0 >= 0',
		 '1 a1 -  b2 >= 0',
		 'b2 - a1 >= 0',
		 '3 a2 -  b0 - b2 - b3 >= 0',
		 'b0 + b2 + b3 - a2 >= 0',
		 '1 a3 -  b0 >= 0',
		 'b0 - a3 >= 0']
		>>>
		>>>
		"""
		return Extension.ForwardDiff_LinearLayer(BasicTools.transpose(M), V_out, V_in)	
	def genConstraints_backwardkeyrecovery(r): 
		Input_PN_diff = BasicTools.VarGen("x_InPN",r)
		Input_round_diff = BasicTools.VarGen("x_InMC",r)
		Output_round_diff = BasicTools.VarGen("x_InMC",r+1)
		
		Input_PN_val = BasicTools.VarGen("y_InPN",r)
		Input_round_val = BasicTools.VarGen("y_InMC",r)
		Output_round_val = BasicTools.VarGen("y_InMC",r+1)
		
		Input_PN_h = BasicTools.VarGen("h_InPN",r)
		Input_round_h = BasicTools.VarGen("h_InMC",r)
		Output_round_h = BasicTools.VarGen("h_InMC",r+1)
		
		DeterminationVar = BasicTools.VarGen("t",r)	     
		
		Constr = []
		#for j in range(4):
		#	Constr = Constr + Extension.LinearLayer(CRAFT.Inv_MC, CRAFT.column(Input_PN_diff,j), CRAFT.column(Input_round_diff,j))
		Constr = Constr + BasicTools.equalConstraints(CRAFT.PN(Input_PN_diff), Output_round_diff)
		for j in range(4):
			#Constr = Constr + Extension.LinearLayer(BasicTools.transpose(CRAFT.MC), CRAFT.column(Input_PN_val,j), CRAFT.column(Input_round_val,j))
			Constr = Constr + Extension.Determination_Decision(CRAFT.MC, CRAFT.column(Input_round_val,j), CRAFT.column(Input_PN_val,j), CRAFT.column(DeterminationVar,j))
		Constr = Constr + BasicTools.equalConstraints(CRAFT.PN(Input_PN_h), Output_round_val)
		#Constr = Constr + BasicTools.greaterConstraints(Output_round_val, Output_round_diff)
		Constr = Constr + BasicTools.greaterConstraints(Input_PN_val, Input_PN_diff)
		Constr = Constr + BasicTools.greaterConstraints(Input_PN_val, Input_PN_h)
		Constr = Constr + BasicTools.greaterConstraints_xyz(Input_PN_h,Input_PN_diff,Input_PN_val)		
		return Constr		
		
	def genConstraints_backwardkeyrecoveryLastR(r): 
		Input_round_diff = BasicTools.VarGen("x_InMC",r)
		Input_PN_diff = BasicTools.VarGen("x_InPN",r)
		Output_round_diff = BasicTools.VarGen("x_InMC",r+1)
		Input_round_val = BasicTools.VarGen("y_InMC",r)
		Input_PN_val = BasicTools.VarGen("y_InPN",r)
		Output_round_val = BasicTools.VarGen("y_InMC",r+1)     
		Input_MC_McRel = BasicTools.VarGen("w",r)
		Constr = []
		for j in range(4):
			#Constr = Constr + Extension.LinearLayer(CRAFT.Inv_MC, CRAFT.column(Input_PN_diff,j), CRAFT.column(Input_round_diff,j))		
			Constr = Constr + Extension.McRelImprovementBackward(CRAFT.column(Input_PN_diff,j), CRAFT.column(Input_MC_McRel,j))
		return Constr
				
	def genConstraints_forwardkeyrecovery(r): 
		Input_round_diff = BasicTools.VarGen("x_InMC",r)
		Input_PN_diff = BasicTools.VarGen("x_InPN",r)
		Output_round_diff = BasicTools.VarGen("x_InMC",r+1)
		Input_round_val = BasicTools.VarGen("y_InMC",r)
		Input_PN_val = BasicTools.VarGen("y_InPN",r)
		Output_round_val = BasicTools.VarGen("y_InMC",r+1)
		Input_round_h = BasicTools.VarGen("h_InMC",r)
		Input_PN_h = BasicTools.VarGen("h_InPN",r)
		Output_round_h = BasicTools.VarGen("h_InMC",r+1)
		DeterminationVar = BasicTools.VarGen("t",r)
	     
		Constr = []
		#for j in range(4):
		#	Constr = Constr + Extension.LinearLayer(CRAFT.MC, CRAFT.column(Input_round_diff,j), CRAFT.column(Input_PN_diff,j))
		Constr = Constr + BasicTools.equalConstraints(CRAFT.PN(Input_PN_diff), Output_round_diff)
		for j in range(4):
			#Constr = Constr + Extension.LinearLayer(BasicTools.transpose(CRAFT.Inv_MC), CRAFT.column(Input_round_h,j), CRAFT.column(Input_PN_h,j))
			Constr = Constr + Extension.Determination_Decision(CRAFT.Inv_MC, CRAFT.column(Input_PN_h,j), CRAFT.column(Input_round_h,j), CRAFT.column(DeterminationVar,j))
		Constr = Constr + BasicTools.equalConstraints(Output_round_h, CRAFT.PN(Input_PN_val))
		Constr = Constr + BasicTools.greaterConstraints(Input_PN_val, Input_PN_diff)
		Constr = Constr + BasicTools.greaterConstraints(Input_PN_val, Input_PN_h)
		Constr = Constr + BasicTools.greaterConstraints_xyz(Input_PN_h,Input_PN_diff,Input_PN_val)
		
		return Constr
	def genConstraints_forwardkeyrecoveryFirstR(r): 
		Input_round_diff = BasicTools.VarGen("x_InMC",r)
		Input_PN_diff = BasicTools.VarGen("x_InPN",r)
		Output_round_diff = BasicTools.VarGen("x_InMC",r+1)
		Input_round_val = BasicTools.VarGen("y_InMC",r)
		Input_PN_val = BasicTools.VarGen("y_InPN",r)
		Output_round_val = BasicTools.VarGen("y_InMC",r+1)
		Output_MC_McRel = BasicTools.VarGen("w",r+1)
	     
		Constr = []
		for j in range(4):
			#Constr = Constr + Extension.LinearLayer(CRAFT.MC, CRAFT.column(Input_round_diff,j), CRAFT.column(Input_PN_diff,j))	
			Constr = Constr + Extension.McRelImprovementForward(CRAFT.column(Input_round_diff,j), CRAFT.column(Output_MC_McRel,j))		
		return Constr		
        	    		    
if __name__ == '__main__':
	print("HI")
	const=[]
	fileobj = open(filename_model, "w")	
	#----------------------------------------2
	fileobj.write("Minimize\n")
	fileobj.write("Dummy\n")
	fileobj.write("Subject to\n")
	#fileobj.write(' + '.join( ['x_InMC_3_' + str(j) for j in range(0,16)]))
	#fileobj.write(" = 5\n x_InMC_3_0 + x_InMC_3_4 + x_InMC_3_8 + x_InMC_3_12 + x_InMC_3_3 = 5\n")
	fileobj.write("Dummy - ")
	fileobj.write(' - '.join( ['40 p0' + '_' + str(i) for i in range(PreRound*4,(PreRound+DistRound)*4)] + ['80 p1' + '_' + str(i) for i in range(PreRound*4,(PreRound+DistRound)*4)]))#p
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['40 p0' + '_' + str(i) for i in range((PreRound+DistRound)*4,(PreRound+DistRound+PostRound)*4)] + ['80 p1' + '_' + str(i) for i in range((PreRound+DistRound)*4,(PreRound+DistRound+PostRound)*4)]))#p_out
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['40 k_pre_0_' + str(j) for j in range(16)]))#k_in
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['40 k_pre_1_' + str(j) for j in range(16)]))#k_in
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['41 t' + '_' + str(i) + '_' + str(j) for i in range(0,PreRound-1) for j in range(16)]))#t in k_in
	
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['t' + '_' + str(i) + '_' + str(j) for i in range(PreRound+DistRound+1,PreRound+DistRound+PostRound) for j in range(16)]))#t in k_out
	
	fileobj.write(" >= 0 \n")
	
	
	fileobj.write(' + '.join( ['x_InMC_' + str(PreRound) + '_' + str(j) for j in range(0,16)]))#delta_in
	fileobj.write(" >= 1\n")
	"""for i in range(DistRound):
		fileobj.write(' + '.join( ['x_InMC_' + str(PreRound+i) + '_' + str(j) for j in range(0,16)]))
		fileobj.write(" <= 13\n")"""
		
	fileobj.write(' + '.join( ['40 p0' + '_' + str(i) for i in range(0,(PreRound+DistRound+PostRound)*4)] + ['80 p1' + '_' + str(i) for i in range(0,(PreRound+DistRound+PostRound)*4)]))#p
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['40 x_InMC_' + str(PreRound) + '_' + str(j) for j in range(0,16)]))#delta_in
	fileobj.write(" <= 630 \n") #data complexity p-delta_in<n=64
	####################################
	fileobj.write("Dummy - ")
	fileobj.write(' - '.join( ['40 p0' + '_' + str(i) for i in range(PreRound*4,(PreRound+DistRound)*4)] + ['80 p1' + '_' + str(i) for i in range(PreRound*4,(PreRound+DistRound)*4)]))#p
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['40 p0' + '_' + str(i) for i in range(0,PreRound*4)] + ['80 p1' + '_' + str(i) for i in range(0,PreRound*4)]))#p_in
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['40 k_post_0_' + str(j) for j in range(16)]))#k_out
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['40 k_post_1_' + str(j) for j in range(16)]))#k_out
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['41 t' + '_' + str(i) + '_' + str(j) for i in range(PreRound+DistRound+1,PreRound+DistRound+PostRound) for j in range(16)]))#t in k_out
	
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['t' + '_' + str(i) + '_' + str(j) for i in range(0,PreRound-1) for j in range(16)]))#t in k_in
	
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['40 x_InMC_' + str(PreRound+DistRound) + '_' + str(j) for j in range(0,16)]))#delta_out
	fileobj.write(" + ")
	fileobj.write(' + '.join( ['40 x_InMC_' + str(PreRound) + '_' + str(j) for j in range(0,16)]))#delta_in
	fileobj.write(" >= 0 \n")
	
	fileobj.write("Dummy2 - ")
	fileobj.write(' - '.join( ['40 p0' + '_' + str(i) for i in range(PreRound*4,(PreRound+DistRound)*4)] + ['80 p1' + '_' + str(i) for i in range(PreRound*4,(PreRound+DistRound)*4)]))#p
	fileobj.write(" = 0 \n")
	
	fileobj.write("Dummy3 - ")
	fileobj.write(' - '.join( ['40 x_InMC_' + str(PreRound+DistRound) + '_' + str(j) for j in range(0,16)]))#delta_out
	fileobj.write(" = 0 \n")
	
	#fileobj.write("Dummy3 <= 120 \n")
	
	fileobj.write("Dummy4 - ")
	fileobj.write(' - '.join( ['40 x_InMC_' + str(PreRound) + '_' + str(j) for j in range(0,16)]))#delta_in
	fileobj.write(" = 0 \n")
	
	fileobj.write("Dummy5 - ")
	fileobj.write(' - '.join( ['40 k_post_0_' + str(j) for j in range(16)]))#k_out
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['40 k_post_1_' + str(j) for j in range(16)]))#k_out
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['41 t' + '_' + str(i) + '_' + str(j) for i in range(PreRound+DistRound+1,PreRound+DistRound+PostRound) for j in range(16)]))#t in k_out
	fileobj.write(" = 0 \n")
	
	fileobj.write("Dummy6 - ")
	fileobj.write(' - '.join( ['40 k_pre_0_' + str(j) for j in range(16)]))#k_in
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['40 k_pre_1_' + str(j) for j in range(16)]))#k_in
	fileobj.write(" - ")
	fileobj.write(' - '.join( ['41 t' + '_' + str(i) + '_' + str(j) for i in range(0,PreRound-1) for j in range(16)]))#t in k_in
	fileobj.write(" = 0 \n")
	
	fileobj.write("Dummy7 - ")
	fileobj.write(' - '.join( ['40 p0' + '_' + str(i) for i in range(0,PreRound*4)] + ['80 p1' + '_' + str(i) for i in range(0,PreRound*4)]))#p_in
	fileobj.write(" = 0 \n")
	
	fileobj.write("Dummy8 - ")
	fileobj.write(' - '.join( ['40 p0' + '_' + str(i) for i in range((PreRound+DistRound)*4,(PreRound+DistRound+PostRound)*4)] + ['80 p1' + '_' + str(i) for i in range((PreRound+DistRound)*4,(PreRound+DistRound+PostRound)*4)]))#p_out
	fileobj.write(" = 0 \n")
	
	fileobj.write("Dummy9 - ")
	fileobj.write(' - '.join( ['40 p0' + '_' + str(i) for i in range(0,(PreRound+DistRound+PostRound)*4)] + ['80 p1' + '_' + str(i) for i in range(0,(PreRound+DistRound+PostRound)*4)]))#p
	fileobj.write(" + ")
	fileobj.write(' + '.join( ['40 x_InMC_' + str(PreRound) + '_' + str(j) for j in range(0,16)]))#delta_in
	fileobj.write(" = 0 \n") #data complexity p-delta_in<n=64
	####################################	
	fileobj.write(' + '.join( ['40 p0' + '_' + str(i) for i in range(PreRound*4,(PreRound+DistRound)*4)] + ['80 p1' + '_' + str(i) for i in range(PreRound*4,(PreRound+DistRound)*4)]))#p
	fileobj.write(" + ")
	fileobj.write(' + '.join( ['40 x_InMC_' + str(PreRound+DistRound) + '_' + str(j) for j in range(0,16)]))
	fileobj.write(" <= 630\n")
	const+=Truncated.genEncSubjection(DistRound,PreRound,PostRound)
	for i in range(0,PreRound+DistRound+PostRound):
		const+=BasicTools.equalConstraints(CRAFT.PN(BasicTools.VarGen("x_InMC",i+1)), BasicTools.VarGen("x_InPN",i))
	for i in range(0,PreRound-2):
		const+=BasicTools.greaterConstraints((BasicTools.VarGen("z_InPN",i)), BasicTools.VarGen("x_InPN",i))
		const+=BasicTools.greaterConstraints((BasicTools.VarGen("z_InPN",i)), BasicTools.VarGen("y_InPN",i))
	for i in range(PreRound-2,PreRound-1):
		#const+=BasicTools.greaterConstraints((BasicTools.VarGen("z_InPN",i)), BasicTools.VarGen("x_InPN",i))
		const+=BasicTools.greaterConstraints((BasicTools.VarGen("z_InPN",i)), BasicTools.VarGen("y_InPN",i))
	for i in range(PreRound+DistRound+2,PreRound+DistRound+PostRound):
		const+=BasicTools.greaterConstraints((BasicTools.VarGen("z_InMC",i)), BasicTools.VarGen("x_InMC",i))
		#const+=BasicTools.greaterConstraints((BasicTools.VarGen("z_InMC",i)), BasicTools.VarGen("y_InMC",i))
		const+=BasicTools.greaterConstraints((BasicTools.VarGen("z_InMC",i)), BasicTools.VarGen("h_InMC",i))
	for i in range(PreRound+DistRound+1,PreRound+DistRound+2):
		#const+=BasicTools.greaterConstraints((BasicTools.VarGen("z_InMC",i)), BasicTools.VarGen("x_InMC",i))
		#const+=BasicTools.greaterConstraints((BasicTools.VarGen("z_InMC",i)), BasicTools.VarGen("y_InMC",i))
		const+=BasicTools.greaterConstraints((BasicTools.VarGen("z_InMC",i)), BasicTools.VarGen("h_InMC",i))
					
	c = []
	x=BasicTools.VarGen("x_InPN",PreRound-1)
	y=BasicTools.VarGen("y_InMC",PreRound-1)
	w=BasicTools.VarGen("w",PreRound-1)
	for i in range(0, 16):#y=x&w
		c = c + [x[i] + ' + ' + w[i] + ' - ' + y[i] + ' >= 0']
		c = c + [x[i] + ' - ' + w[i] + ' - ' + y[i] + ' >= -1']
		c = c + [' - ' + x[i] + ' + ' + w[i] + ' - ' + y[i] + ' >= -1']
		c = c + [' - ' + x[i] + ' - ' + w[i] + ' + ' + y[i] + ' >= -1']
	const+=c	
	for i in range(PreRound-1):
		const+=Extension.genConstraints_backwardkeyrecovery(i)
	const+=Extension.genConstraints_backwardkeyrecoveryLastR(PreRound-1)
	
	
	for i in range(1,PreRound-1):
		zprime=['zp_InPN_' + str(i-1) + '_' + str(j) for j in range(0,16)]
		z=['z_InPN_' + str(i-1) + '_' + str(j) for j in range(0,16)]
		t=['t_' + str(i-1) + '_' + str(j) for j in range(0,16)]
		for j in range(16):
			eqn=[]
			eqn.append(z[j] + ' - ' + t[j] + ' - ' + zprime[j] + ' = 0 ')
			#eqn.append(z[j] + ' - ' + zprime[j] + ' = 0 ')
			const += eqn
			
	zprime=['zp_InPN_' + str(PreRound-2) + '_' + str(j) for j in range(0,16)]
	h=['h_InPN_' + str(PreRound-2) + '_' + str(j) for j in range(0,16)]
	t=['t_' + str(PreRound-2) + '_' + str(j) for j in range(0,16)]
	for j in range(16):
		eqn=[]
		eqn.append(h[j] + ' - ' + t[j] + ' - ' + zprime[j] + ' = 0 ')
		#eqn.append(h[j] + ' - ' + zprime[j] + ' = 0 ')
		const += eqn
		
	"""zprime=['zp_InMC_' + str(PreRound+DistRound) + '_' + str(j) for j in range(0,16)]
	z=['z_InMC_' + str(PreRound+DistRound) + '_' + str(j) for j in range(0,16)]
	for j in range(16):
		eqn=[]
		eqn.append(z[j] + ' - ' + zprime[j] + ' = 0 ')
		const += eqn	"""
	for i in range(PreRound+DistRound+2,PreRound+DistRound+PostRound):
		zprime=['zp_InMC_' + str(i) + '_' + str(j) for j in range(0,16)]
		z=['z_InMC_' + str(i) + '_' + str(j) for j in range(0,16)]
		t=['t_' + str(i) + '_' + str(j) for j in range(0,16)]
		for j in range(16):
			eqn=[]
			eqn.append(z[j] + ' - ' + t[j] + ' - ' + zprime[j] + ' = 0 ')
			#eqn.append(z[j] + ' - ' + zprime[j] + ' = 0 ')
			const += eqn
	
	zprime=['zp_InMC_' + str(PreRound+DistRound+1) + '_' + str(j) for j in range(0,16)]
	h=['h_InMC_' + str(PreRound+DistRound+1) + '_' + str(j) for j in range(0,16)]
	t=['t_' + str(PreRound+DistRound+1) + '_' + str(j) for j in range(0,16)]
	for j in range(16):
		eqn=[]
		eqn.append(h[j] + ' - ' + t[j] + ' - ' + zprime[j] + ' = 0 ')
		#eqn.append(h[j] + ' - ' + zprime[j] + ' = 0 ')
		const += eqn
			
	k0_pre=BasicTools.VarGen("k_pre",0)
	k1_pre=BasicTools.VarGen("k_pre",1)
	k0_post=BasicTools.VarGen("k_post",0)
	k1_post=BasicTools.VarGen("k_post",1)	
	if PreRound%2 == 0:	
		for i in range((PreRound//2)-1):
			for j in range(16):
				eqn=[]
				eqn.append(k0_pre[j] + ' - zp_InPN_' + str(2*i) + '_' + str(j) + ' >= 0')
				const += eqn
				eqn=[]
				eqn.append(k1_pre[j] + ' - zp_InPN_' + str(2*i+1) + '_' + str(j) + ' >= 0')
				const += eqn
		for j in range(16):
			eqn=[]
			eqn.append(k0_pre[j] + ' - zp_InPN_' + str((PreRound//2)) + '_' + str(j) + ' >= 0')
			const += eqn
	elif PreRound%2 == 1:	
		for i in range(((PreRound+1)//2)-1):
			for j in range(16):
				eqn=[]
				eqn.append(k0_pre[j] + ' - zp_InPN_' + str(2*i) + '_' + str(j) + ' >= 0')
				const += eqn
				eqn=[]
				eqn.append(k1_pre[j] + ' - zp_InPN_' + str(2*i+1) + '_' + str(j) + ' >= 0')
				const += eqn
				
	if PostRound%2 == 0:	
		for i in range((PostRound//2)-1):
			for j in range(16):
				eqn=[]
				eqn.append(k0_post[j] + ' - zp_InMC_' + str(2*i+PreRound+DistRound+1) + '_' + str(j) + ' >= 0')
				const += eqn
				eqn=[]
				eqn.append(k1_post[j] + ' - zp_InMC_' + str(2*i+1+PreRound+DistRound+1) + '_' + str(j) + ' >= 0')
				const += eqn
		for j in range(16):
			eqn=[]
			eqn.append(k0_post[j] + ' - zp_InMC_' + str(PostRound+PreRound+DistRound-1) + '_' + str(j) + ' >= 0')
			const += eqn
	elif PostRound%2 == 1:	
		for i in range(((PostRound+1)//2)-1):
			for j in range(16):
				eqn=[]
				eqn.append(k0_post[j] + ' - zp_InMC_' + str(2*i+PreRound+DistRound+1) + '_' + str(j) + ' >= 0')
				const += eqn
				eqn=[]
				eqn.append(k1_post[j] + ' - zp_InMC_' + str(2*i+1+PreRound+DistRound+1) + '_' + str(j) + ' >= 0')
				const += eqn

		
	
	####################################
	#const+=BasicTools.equalConstraints(BasicTools.VarGen("x_InMC",PreRound+DistRound+1), BasicTools.VarGen("h_InMC",PreRound+DistRound+1))
	c = []
	x=BasicTools.VarGen("x_InMC",PreRound+DistRound+1)
	h=BasicTools.VarGen("h_InMC",PreRound+DistRound+1)
	w=CRAFT.PN(BasicTools.VarGen("w",PreRound+DistRound+1))
	for i in range(0, 16):#h=x&w
		c = c + [x[i] + ' + ' + w[i] + ' - ' + h[i] + ' >= 0']
		c = c + [x[i] + ' - ' + w[i] + ' - ' + h[i] + ' >= -1']
		c = c + ['- ' + x[i] + ' + ' + w[i] + ' - ' + h[i] + ' >= -1']
		c = c + ['- ' + x[i] + ' - ' + w[i] + ' + ' + h[i] + ' >= -1']
	const+=c
	const+=Extension.genConstraints_forwardkeyrecoveryFirstR(PreRound+DistRound)
	for i in range(PreRound+DistRound+1,PreRound+DistRound+PostRound):
		const+=Extension.genConstraints_forwardkeyrecovery(i)
	####################################		
	for c in const:
		fileobj.write(str(c))
		fileobj.write("\n")
	fileobj.write("Binary\n")
	for c in const:
		Var = BasicTools.getVariables(c)
		#print(c)
		#print(Var)
		for v in Var:
			fileobj.write(v)
			fileobj.write("\n")
	fileobj.write("General\n")
	fileobj.write("Dummy\n")
	fileobj.write("Dummy2\n")
	fileobj.write("Dummy3\n")
	fileobj.write("Dummy4\n")
	fileobj.write("Dummy5\n")
	fileobj.write("Dummy6\n")
	fileobj.write("Dummy7\n")
	fileobj.write("Dummy8\n")
	fileobj.write("Dummy9\n")
	fileobj.write("End")
	fileobj.close()
	
	
	m = read(filename_model)
	m.Params.threads=192
	#m.Params.PoolSolutions=2000000000
	m.optimize()
	#print("m.solcount=*************")
	#print(m.solcount)	
	m.write("t.sol")
