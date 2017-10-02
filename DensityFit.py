import os
from scipy.optimize import curve_fit

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src/")

execfile("Schrodinger.py")
execfile("Poisson.py")
execfile("DOS.py")
execfile("Files.py")


def expFit(x, a, b, d, L):
	return d/(x+L)**(12/7.)

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/Results16/")

name = 19


T = readIteration(str(name))
well_energies = T[0]
Ef = T[1]
n2d = T[2]
z = T[3]
zrev = [10**9*z[i] for i in range(len(z))]
V = T[4]
	
P = rho3d_band(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d,1)
	
Vtild = newpotential(P,n2d,z)
	
stepzzz = z[1]-z[0]
	
Pscaled = [p/stepzzz * 10**(-27) for p in P]

popt, pcov = curve_fit(expFit,zrev,Pscaled)
	
plt.plot(zrev,Pscaled,color = "Black",label = 'data')
print(popt)
print(pcov)

Pfit = [expFit(zrev[i],*popt) for i in range(len(z))]

plt.plot(zrev,Pfit, 'r-', label='fit')

plt.xlabel("z (nm)")
plt.ylabel("LDOS ($nm^{-3}$)")
plt.title("Density near the interface")
#plt.savefig('Density.pdf')


plt.legend()
plt.show()