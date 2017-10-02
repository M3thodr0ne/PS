execfile("Header.py")
execfile("Constantes.py")

#the function giving epsilonr E, used to determine E in shooting poisson
def W(E):
	B = 25462.
	E0 = 82213.
	return E*(1+B/((1+(E/E0)**2)**(1./3)))

#solve the implicit equation epsilon(E)*E = e*n2d/epsilon0
def shooting_poisson(W_func,n2d,iniguess):
	el = 1.6*10**(-19)
	func = lambda x : el*n2d/epsilon0 - W_func(x)
	E_initial_guess = iniguess
	Esol = fsolve(func,E_initial_guess)
	return Esol[0]
    
#return the primitive of the function T, cancelling on 0 and ranging until the index index with the step of integration step
def integration(T,index,step):
    return sum(T[i] for i in range(index))*step

def double_integration(T,index,step):
    return sum(integration(T,i,step) for i in range(index))*step

B = 25462.
E0 = 82213.

#to compute the relative permittivity of the material
def epsilonr(E):
    B = 25462.
    E0 = 82213.
    return 1+B/((1+(E/E0)**2)**(1./3))

#solving the Poisson equation from the LDOS rho and the given 2D density n2d
def newpotential(rho,n2d,z):
	#compute the step of z to have the step of integration
	step = z[1]-z[0]
	#initial value of the electric field computed from the 2D density of the gas
	electric_field = shooting_poisson(W,n2d,E0/2)
	DI = [0 for i in range(len(z))]
	I = [0 for i in range(len(z))]
	s = 0
	origin_field = electric_field
	for i in range(len(z)):
		#integration of the LDOS to compute the displacement field
		s += rho[i]
		I[i] = s
	s = 0
	for i in range(len(z)):
		if i >=1:
        	#compute the electric field from the value of the displacement field
			electric_field = shooting_poisson(W,n2d-I[i],electric_field)
			
            #integrating the electric field to have the potential
			s += electric_field*step
			DI[i] = s
		else:
			s += electric_field*step
			DI[i] = s
    #satisfying the boundary condition
	V = [-DI[0]+DI[i] for i in range(len(z))]
	return V

#computing the corrected potential
def potential_corrected(Vanc,Vnew,ksi):
    V = [(1-ksi)*Vanc[i] + ksi*Vnew[i] for i in range(len(Vanc))]
    return V
    
#computing the error in the potential
def potential_error(Vanc,Vnew):
    I = sum(((Vanc[i+1]-Vnew[i+1])/Vanc[i+1])**2 for i in range(len(Vanc)-1))
    return I/(len(Vanc))