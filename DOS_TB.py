import os

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src/")

execfile("Header.py")

from scipy.special import ellipk
import mpmath

tt = [-1+0.002*i for i in range(1000)]
ttt = [ellipk(pp) for pp in tt]

#return the correct value of
def redE(gamfactor,E):
	if abs(E) < 1-gamfactor:
		return sqrt(4*gamfactor/(abs((1+gamfactor)**2-E**2)))
	else:
		return sqrt(abs(((1+gamfactor)**2-E**2))/(4*gamfactor))
		
def DOS(gamfactor,E):
	kred = redE(gamfactor,E)
	if abs(E) < 1-gamfactor:
		return 1/sqrt(gamfactor) * kred * (mpmath.ellipk(kred+10**(-3)*1j)).real
	elif abs(E) < 1+gamfactor:
		return sqrt(gamfactor)*(mpmath.ellipk(kred+10**(-3)*1j)).real
	else:
		return 0
	
gamfactor = 1

tabE = [-2.000+0.0002*i for i in range(20000)]
tabK = [0 for ene in tabE]

for i in range(len(tabE)):
	tabK[i] = DOS(gamfactor,tabE[i])


plot(tabE,tabK)

print(sum(tabK))


plt.show()