execfile("Header.py")
execfile("Constantes.py")
import warnings

#to make a continuous function from a table of values
def interpolation(V,x,stepx):
    N = int(x/stepx)
    t = x/stepx - N
    if N >=len(V)-1:
        return V[-1]
    return (1-t)*V[N+1] + t*V[N]

#to limit spatially the research of the eigenmodes of the well (numerical stability), when the potential excedes the energy the solution can diverge ...
def lastindex(E,V):
    for i in range(len(V)):
        if E < V[i]:
            return i
    return len(V)

def schrodinger2(sol,z,stepz,E,V,mt):
	dpsi = sol[1]
	el = 1.6*10**(-19)
	g = mt/(h**2)*el*(- E + interpolation(V,z,stepz))*sol[0]
	return [dpsi,g]

#define the schrodinger equation in a potential V. The third component is to tell the program to stop computing in diverging cases
#When the solver cuts to avoid divergence, the wavefunctions are let to their final value for the end of the z tab
def schrodinger(sol,z,stepz,E,V,mt):
	dpsi = sol[1]
	el = 1.6*10**(-19)
	if - E + interpolation(V,z,stepz) > 0 and z > 4*stepz:
		if signum(sol[0]) == signum(sol[1]) or sol[2]!= 0:			
			return [0,0,1]
		else:
			g = mt/(h**2)*el*(- E + interpolation(V,z,stepz))*sol[0]
	else:
		g = mt/(h**2)*el*(- E + interpolation(V,z,stepz))*sol[0]
	return [dpsi,g,0]

#solver of the schrodinger equation, take into account the limitations of lastindex
def short_solver2(energy,V,z,m):
	CI_tableau = [0,0,0]
	CI_tableau[0] = 0.
	CI_tableau[1] = 0.1
	CI_tableau[2] = 0
	stepz = abs(z[1]-z[0])
    #numerical stability
	j = lastindex(energy,V)
	y_ode = [0 for i in range(len(z))]
    #spatial limitation of the solution
	ind = max(min(int(1.3*j),len(z)),200)
	#ind = len(z)
    #spatial limitation of the z coordinate to solve only where the solution is stable numerically

	bo = 0 
	warnings.filterwarnings('error')
	while bo == 0:
		bo = 1
		try:
			partialz = [z[k] for k in range(ind)]
			sol_ode = odeint(schrodinger,CI_tableau,partialz,(stepz,energy,V,m))
		except:
			ind -= 10
			bo = 0
			print(ind)
	t_ode = [0 for i in range(len(z))]
	#maybe all these steps aint necessary anymore
	for i in range(ind):
		y_ode[i] = sol_ode[i,0]
		t_ode[i] = abs(y_ode[i])
	for i in range(ind,len(y_ode)):
		y_ode[i] = y_ode[ind-1]
	y_ode /= sqrt(norm(y_ode))
	warnings.filterwarnings('default')
	return y_ode
	
def short_solver(energy,V,z,m):
	CI_tableau = [0,0,0]
	CI_tableau[0] = 0.
	CI_tableau[1] = 0.1
	stepz = abs(z[1]-z[0])
    
	y_ode = [0 for i in range(len(z))]
    
	sol_ode = odeint(schrodinger,CI_tableau,z,(stepz,energy,V,m))
	#maybe all these steps aint necessary anymore
	for i in range(len(z)):
		y_ode[i] = sol_ode[i,0]
	normy = norm(y_ode)
	if normy != 0:
		for i in range(len(z)):
			y_ode[i] /= normy
	return y_ode

def short_solverderiv(energy,V,z,m):
	CI_tableau = [0,0,0]
	CI_tableau[0] = 0.
	CI_tableau[1] = 0.1
	stepz = z[1]-z[0]
    
	y_ode = [0 for i in range(len(z))]
	y_ode2 = [0 for i in range(len(z))]
    
	sol_ode = odeint(schrodinger,CI_tableau,z,(stepz,energy,V,m))
	#maybe all these steps aint necessary anymore
	for i in range(len(z)):
		y_ode[i] = sol_ode[i,0]
		y_ode2[i] = sol_ode[i,1]
	normy = norm(y_ode)
	if normy != 0:
		for i in range(len(z)):
			y_ode2[i] /= normy
	return y_ode2
    
    
#to find the eigenvalues of the schrodinger equation which are vanishing at infinity, uses tricks to count sign changes of the wavefunction. dichotomy. Searching for the i-th evanescent state in the potential V, seeking from the energy Egs, with the accuracy acc and the mass m for the particle
def finding_eigenval_i(z,Egs,V,m,i,acc):
	#initializing the dichotomy loop
    Emin = Egs
    Emax = Egs+0.25
    psi = short_solver(Emin,V,z,m)
    nmin = counting_sign(psi)
    psi = short_solver(Emax,V,z,m)
    #how many sign change are there in the wf (careful it must be a real function !!!)
    nmax = counting_sign(psi)
    #end of the initialization of the dichotomy
    while nmax <= i:
        Emin = Emax
        Emax+=0.25
        psi = short_solver(Emax,V,z,m)
        nmax = counting_sign(psi)
    #dichotomy loop
    while Emax - Emin > acc:
        Emil = (Emax + Emin)/2.
        nmil = counting_sign(short_solver(Emil,V,z,m))
        if nmil <= i:
            Emin = Emil
        else:
            Emax = Emil
    return Emax
    
#compute the eigenstate associated to an energy level E, in the approximation of a localized state. fills the end of the values of the eigenfunctions if needed. Need to use renorm then to have proper "vanishing" eigenfunctions
def eigenstate(z,E,V,m):
    CI_tableau = np.zeros(3)
    CI_tableau[0] = 0.
    CI_tableau[1] = 0.1
    stepz = z[1]-z[0]
    y_ode = short_solver(E,V,z,m)
    length = len(y_ode)
    y_ode = [y_ode[i] for i in range(length)]
    nory = norm(y_ode)
    y_ode = [y/nory for y in y_ode]
    return y_ode
      

#see Schrodinger for more explanation
#set to zero all the value that are equal at the end, then normalize the wf
#use renorm if the function has the correct number of zeros
def renorm(psi):
	N = counting_sign(psi)
	ite = 0
	psirenorm = [psi[i] for i in range(len(psi))]
	ini_value = psi[0]
	trig = 0
	i = 0
	while psi[len(psi)-1-i] == psi[-1]:
		psirenorm[-i] = 0
		i+=1
	if norm(psirenorm) != 0:
		npsi = norm(psirenorm)
		psirenorm = [p/npsi for p in psirenorm]
	return psirenorm

def indexlastzero(psi):
	index0 = len(psi)-1
	inisign = signum(psi[-1])
	while signum(psi[index0]) == inisign:
		index0 -= 1

#see Schrodinger for more explanation
#set to zero all the value that are after the last zero of the wf
#use renorm if the function has a zero in excess
def renorm2(psi):
	ite = 0
	psirenorm = [psi[i] for i in range(len(psi))]
	ini_value = psi[0]
	trig = 0
	i = 1
	while signum(psi[-i]) == signum(psi[-1]):
		psirenorm[-i] = 0
		i+=1
	if norm(psirenorm) != 0:
		npsi = norm(psirenorm)
		psirenorm = [p/npsi for p in psirenorm]
	return psirenorm

def indexlastzero(psi):
	index0 = len(psi)-1
	inisign = signum(psi[-1])
	while signum(psi[index0]) == inisign:
		index0 -= 1

#to count the sign change of an wavefunction   
def counting_sign(psi):
    counter = 0
    sign_ini = signum(psi[1])
    for i in range(len(psi)):
        if signum(psi[i]) != sign_ini:
            sign_ini = signum(psi[i])
            counter += 1
    return counter

#sign of a real number, could be -1 or 1, must check that it does not change anything
def signum(x):
    if x >= 0:
        return 1
    return -1

#do the computation of the ith eigenfunction of the well by taking into account the number of zeros. That is if the function has more zeros than expected then cutting the end. For the first eigenvalue for instance, if it has a zero then cut after the zero. If it has no zero then cut when it becomes flat
def ith_eigenfunction(energy,i,mz,V,z):
	psi = eigenstate(z,energy,V,mz)
	if counting_sign(psi) <= i:
		return renorm(psi)
	else:
		return renorm2(psi)

#compute the index of the point where E = V
def VequalsEindex(V,E):
	index = 0
	while V[index] - E < 0:
		index += 1
		if index == len(V)-1:
			return len(V)
	return index
	
def VequalsEmaxindex(V,E,z,mz):
	index = 0
	index0 = VequalsEindex(V,E)
	while (2*mz*el*(V[index0+index]-E))/h**2 * (z[index0+index]-z[index0])**2 < 900:
		index+=1
		if index+index0 >=len(V)-1:
			return len(V)
	return index+index0

#to solve the part right to the point where V = E of the Schrodinger equation
def rightSolver(energy,V,z,m):
	index = VequalsEindex(V,energy)
	zreverse = [z[len(z)-1-i] for i in range(len(z)-index)]
	
	CI_tableau = [0.1,-0.1]
	
	sol_odetild = odeint(schrodinger2,CI_tableau,zreverse,(stepz,energy,V,m))
	y_odetild = sol_odetild[:,0]
	y_odetild2 = sol_odetild[:,1]
	
	normtild = norm(y_odetild)
	
	for i in range(index):
		y_odetild[i] /= normtild
		y_odetild2[i] /= normtild
	
	return [y_odetild,y_odetild2]
	
#to solve the part right to the point where V = E of the Schrodinger equation
def rightSolver2(energy,V,z,m):

	index = VequalsEindex(V,energy)
	stepz = abs(z[1]-z[0])
	index2 = VequalsEmaxindex(V,energy,z,m)
	zreverse = [z[index2-i] for i in range(index2-index)]
	
	Vreverse = [V[index2-i] for i in range(index2-index)]
	
	CI_tableau = [0.1,-0.1]
	
	sol_odetild = odeint(schrodinger2,CI_tableau,zreverse,(stepz,energy,V,m))
	y_odetild = sol_odetild[:,0]
	y_odetild2 = sol_odetild[:,1]
	
	
	return [y_odetild,y_odetild2]
	
	
#more accurate method to compute the eigenvalues of the problem. Compute the number of sign change, and then if the number of sign change is correct one computes the right part of the wavefunction. Then one looks at the branching of the wf at the contact point 
def finding_accurate_eigen(z,Egs,V,mr,nzeros,acc):
	#initializing the dichotomy loop
	Emin = Egs
	Emax = Egs+0.25
	psi = short_solver(Emin,V,z,mr)
	nmin = counting_sign(psi)
	psi = short_solver(Emax,V,z,mr)
	#how many sign change are there in the wf (careful it must be a real function !!!)
	nmax = counting_sign(psi)
	#end of the initialization of the dichotomy
	while nmax <= nzeros:
		Emax+=0.25
		psi = short_solver(Emax,V,z,mr)
		nmax = counting_sign(psi)
	#dichotomy loop
	while Emax - Emin > acc:
		Emil = (Emax + Emin)/2.
		psi = short_solver(Emil,V,z,mr)
		nmil = counting_sign(psi)
		if nmil == nzeros:
			Ty = rightSolver2(Emil,V,z,mr)
			yy = Ty[0]
			yyderiv = Ty[1]
			index = VequalsEindex(V,Emil)
			psideriv = short_solverderiv(Emil,V,z,mr)
			const = yy[-1]
			for k in range(len(yy)):
				yy[k] *= psi[index]/const
			#inutile, mais pourra permettre de reconstruire la fonction dans la prochaine fonction
			psicomplete = [0 for j in range(len(z))]
			for k in range(index):
				psicomplete[k] = psi[k]
			for k in range(len(yy)):
				psicomplete[index+k] = yy[len(yy)-1-k]
				
			ncomplete = norm(psicomplete)
			
			for j in range(len(z)):
				psicomplete[j] /= ncomplete
			
			#plot(z,psicomplete,label=Emil)
			
			if (-1)**nzeros * yyderiv[-1]*psi[index]/const < (-1)**nzeros * psideriv[index]:
				Emin = Emil
			else:
				Emax = Emil
		elif nmil < nzeros:
			Emin = Emil
		else:
			Emax = Emil
	return Emil
	
	
def ith_eigenfuncacc(energy,mz,V,z):
	psi = short_solver(energy,V,z,mz)
	Ty = rightSolver2(energy,V,z,mz)
	yy = Ty[0]
	index = VequalsEindex(V,energy)
	const = yy[-1]
	for k in range(len(yy)):
		yy[k] *= psi[index]/const
	psicomplete = [0 for j in range(len(z))]
	for k in range(index):
		psicomplete[k] = psi[k]
	for k in range(len(yy)):
		psicomplete[index+k] = yy[len(yy)-1-k]	
	ncomplete = norm(psicomplete)
			
	for j in range(len(z)):
		psicomplete[j] /= ncomplete
		
	return psicomplete