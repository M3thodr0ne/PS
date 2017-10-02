import os

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src/")

execfile("Poisson.py")


tabN = [10**17*(16+5*i) for i in range(5)]
tabEpsilon = [epsilonr(shooting_poisson(W,tabN[i],E0)) for i in range(len(tabN))]
tabElfield = [shooting_poisson(W,tabN[i],E0) for i in range(len(tabN))]

sizetabvg = 1000
Vgmax = 200.
stepVg = 2*Vgmax/sizetabvg


tabV = [-Vgmax + stepVg*i for i in range(sizetabvg)]
Length = 0.5*10**(-3)

#displacement (converted in cm-2)
def displacement(E):
	return epsilon0*epsilonr(E)*E/el

def deltaN(Vg,length_sample,electric_field):
	Egate = -Vg/length_sample
	return displacement(electric_field+Egate)-displacement(electric_field)-displacement(Egate)

a = 3.9*10**(-10)

tabdensity = [[a**2*(tabN[j]+deltaN(tabV[i],Length,tabElfield[j])) for i in range(len(tabV))] for j in range(len(tabN))]

for k in range(1):
	plot(tabV,tabdensity[k],label = tabN[k]*a**2)

xlabel("Backgate voltage (V)")
ylabel("Density of the 2DEG (eletron per UC)")

#legend()

plt.show()

