import os
chemin = "/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src"
os.chdir(chemin)

execfile("Schrodinger.py")

#compute the DOS of the dispersion relation energy_func taking into account the well energies, the band energies, up to an energy Emax, with the step DE in energy
def DOS_monteCarlo(energy_func,well_energies,band_energies,Emax,DE,mx,my):
	a = 3.9*10**(-10)
    #sampling k
	stepk = 0.0001
	klim = int(1./stepk)
	kf = pi/a
	Nsample = 100000
	countingE = [0 for i in range(int(Emax/DE))]
	for i in range(Nsample):
		kx = randint(-klim,klim)*kf*stepk
		ky = randint(-klim,klim)*kf*stepk
        #tab of energy because there can be several energies if multiple bands
		tabE = energy_func(kx,ky,mx,my,a,well_energies,band_energies)
		for energ in tabE:
			if energ < Emax:
        	#factor 2 for spin, factor 4 for taking into account the full BZ
				countingE[int(energ/DE)] += 1./(Nsample) * (1/a)**2
	return countingE
    
#energy of a tight binding model with bands of energies well energies + band_energies
def energy_TB(kx,ky,mx,my,a,we,be):
	tabE = [1000 for i in range(len(we[0])*len(we))]
	h = 1.06*10**(-34)
	el = 1.6*10**(-19)
	for band_index in range(len(we)):
		for well_index in range(len(we[0])):
			tabE[len(we)*well_index+band_index] = we[band_index][well_index] + be[band_index] + (h/a)**2 *(1-cos(kx*a))/(mx[band_index]*el) + (h/a)**2 *(1-cos(ky*a))/(my[band_index]*el)
	return tabE
	
def energy_2DEG(kx,ky,mx,my,a,we,be):
	tabE = [1000 for i in range(len(we[0])*len(we))]
	for band_index in range(len(we)):
		for well_index in range(len(we[0])):
			tabE[len(we)*well_index+band_index] = we[band_index][well_index] + be[band_index] + (h*kx)**2/(2*mx[band_index]*el) + (h*ky)**2/(2*my[band_index]*el)
	return tabE

#seeking the Fermi energy by dichotomy
def MC_fermiEnergy(energy_func,well_energies,band_energies,Emax,DE,mx,my,nsurf):
	tabE = DOS_monteCarlo(energy_func,well_energies,band_energies,Emax,DE,mx,my)
	tabint = [0 for i in range(len(tabE))]
	jmil = 0
	jmin = 0
	jmax = len(tabE)-1
	s = 0
	for i in range(len(tabE)):
		if i > 0:
			s+=tabE[i-1]
		tabint[i] = s
	nmin = tabint[jmin]
	while abs(jmin-jmax) > 2:
		jmil = int((jmin+jmax)/2)
		nmil = tabint[jmil]
		if nmil > nsurf:
			jmax = jmil
		else:
			jmin = jmil
	return jmil*DE
	
we = [[0.1*i for i in range(8)] for j in range(3)]
be = [0 for k in range(3)]
Emax = 5
DE = 0.0001
nsurf = 1.6*10**18
s=0
tabE = DOS_monteCarlo(energy_TB,we,be,Emax,DE,mx,my)
tabint = [0 for i in range(len(tabE))]
tabEplot = [DE*i for i in range(len(tabint))]
for i in range(len(tabE)):
	if i > 0:
		s+=tabE[i-1]
	tabint[i] = s


plot(tabEplot,tabint)

DE2 = 0.00001
tabE2 = DOS_monteCarlo(energy_TB,we,be,Emax,DE2,mx,my)
tabint2 = [0 for i in range(len(tabE2))]

s = 0

for i in range(len(tabE2)):
	if i > 0:
		s+=tabE2[i-1]
	tabint2[i] = s

tabEplot2 = [DE2*i for i in range(len(tabE2))]

plot(tabEplot2,tabint2,color = "red")

#plt.show()


print(MC_fermiEnergy(energy_TB,we,be,Emax,DE2,mx,my,nsurf))
print(MC_fermiEnergy(energy_TB,we,be,Emax,DE,mx,my,nsurf))
