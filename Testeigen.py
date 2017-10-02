import os
from scipy.special import airy
chemin = "/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src"
os.chdir(chemin)
execfile("Poisson.py")
execfile("Schrodinger.py")
execfile("DOS.py")
execfile("Files.py")


def listeigen(V,n,m,z):
	acc = 0.001
	tabE = [0 for i in range(n)]
	for i in range(n):
		tabE[i] = finding_eigenval_i(z,0,V,m,i,acc)
	return tabE

def schro(sol,z):
	dpsi = sol[1]
	g = (z-1)*sol[0]
	return [dpsi,g]

surf_density = 1.6*10**17

electric_field = shooting_poisson(W,surf_density,E0)
stepz = 10**(-10)
siz = 10**3
z = [stepz * i for i in range(siz)]
V = [electric_field*z[i] for i in range(len(z))]
acc = 0.00005

zplot = [zz*10**9 for zz in z]

K = 1

zplot2 = [2*4/26.*18.5/20.1*zplot[i]-2.3381 for i in range(len(zplot))]

airynorm = sum(airy(zi)[0]**2 for zi in zplot2)

tabairy = [airy(zz)[0]/sqrt(airynorm) for zz in zplot2]

energ = finding_accurate_eigen(z,0,V,7*m,0,acc)
energ = round(energ,6)


psi = ith_eigenfuncacc(energ,7*m,V,z)
plot(zplot,psi,label = "Numeric")

#energ = finding_accurate_eigen(z,0,V,7*m,1,acc)
#energ = round(energ,6)

#psi = ith_eigenfuncacc(energ,7*m,V,z)
#plot(zplot,psi,label = str(10**3*(energ))+" meV")

#energ = finding_accurate_eigen(z,0,V,7*m,2,acc)
#energ = round(energ,6)

#psi = ith_eigenfuncacc(energ,7*m,V,z)
plot(zplot,tabairy,label = "Theory")

plt.title("Shooting method")
plt.xlabel("z (nm)")
plt.ylabel("psi")

plt.legend()

plt.show()