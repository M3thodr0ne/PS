import os

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src/")

execfile("Constantes.py")
execfile("DOS.py")
execfile("Files.py")


def band_densities(well_energies,band_energies,mx,my,Ef):
	densities = [[dens_single_band(well_energies[i][j]+band_energies[i],Ef,mx[i],my[i]) for j in range(len(well_energies[0]))] for i in range(len(well_energies))]
	return densities
	
def density_iterations(file_number):
	T = readIteration(str(file_number))
	well_energies = T[0]
	Ef = T[1]
	BD = band_densities(well_energies,band_energies,mx,my,Ef)
	return [sum(BD[0]),sum(BD[1])]

numberfiles = 18

tabindexfiles = [i for i in range(numberfiles)]
densityxz = [0 for i in range(numberfiles)]
densityyz = [0 for i in range(numberfiles)]
densityxy = [0 for i in range(numberfiles)]

	
os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/Results25/")

for i in range(numberfiles):
	BBD = density_iterations(i)
	densityxz[i] = BBD[1]/(BBD[0]+2*BBD[1])
	densityyz[i] = BBD[1]/(BBD[0]+2*BBD[1])
	densityxy[i] = BBD[0]/(BBD[0]+2*BBD[1])
	
#plot(tabindexfiles,densityxz,label = "Heavy")
#plot(tabindexfiles,densityxy,label = "Light")

tabden = [0.5,1.0,1.6,3.3]
tabxy = [0.24,0.35,0.42,0.53]

tabxz = [0.38,0.32,0.29,0.23]

tabxym = [tabxy[i]*tabden[i] for i in range(len(tabden))]
tabxzm = [tabxz[i]*tabden[i] for i in range(len(tabden))]

tef = [0.0384,0.1258,0.284,0.973]
t0l = [0.0167,0.0648,0.160,0.616]
t1l = [0.0255,0.0951,0.229,0.843]
t2l = [0.0296,0.1084,0.257,0.923]
t3l = [0.0319,0.1145,0.268,0.947]
t0h = [0.0285,0.1054,0.251,0.909]


#loglog(tabden,[tef[i]-t3l[i] for i in range(len(tabden))],label = "Fermi energy")
#plot(tabden,tabxzm,label = "Heavy")

#legend()
#plt.show()

print(h**2/(2*0.7*(3.9*10**(-10))**2*el)/m)
print(h**2/(2*7*(3.9*10**(-10))**2*el)/m)
