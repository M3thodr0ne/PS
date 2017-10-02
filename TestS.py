import os
from scipy.special import airy

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src/")

execfile("Header.py")
execfile("Constantes.py")
execfile("Schrodinger.py")
execfile("Poisson.py")

def schrodingertest(sol,z,stepz,E,V,mt):
	dpsi = sol[1]
	el = 1.6*10**(-19)
	if - E + interpolation(V,z,stepz) > 0 and z > 4*stepz:
		if signum(sol[0]) == signum(sol[1]) or sol[1] == 0:			
			return [0,0]
		else:
			g = mt/(h**2)*el*(- E + interpolation(V,z,stepz))*sol[0]
	else:
		g = mt/(h**2)*el*(- E + interpolation(V,z,stepz))*sol[0]
	return [dpsi,g]
	
	
def testSolve(E):
	n2d = 1.6*10**18
	electric_field = shooting_poisson(W,n2d,E0)
	stepz = 10**(-10)
	siz = 10**3
	z = [stepz * i for i in range(siz)]
	V = [electric_field*z[i] for i in range(len(z))]
	
	m = 9.1*10**(-31)
	
	CI_tableau = [0,0]
	CI_tableau[0] = 0.
	CI_tableau[1] = 0.1
	
	sol_ode = odeint(schrodingertest,CI_tableau,z,(stepz,E,V,m))
	psi = [sol_ode[i,0] for i in range(len(z))]
	plot(z,psi)
	plt.show()
	return psi

tabz = [0.01 * i for i in range(500)]
tabairy = [airy(zz-2.3381)[0]**2 for zz in tabz]

plot(tabz,tabairy)
plt.show()