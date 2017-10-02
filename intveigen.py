import os

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src/")

execfile("Schrodinger.py")
execfile("Poisson.py")
execfile("DOS.py")
execfile("Interaction.py")
execfile("Band_structure.py")

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src/")

execfile("Files.py")


#to plot the potential and the eigenfunctions below the Fermi energy of a file stored under the name name
#what legend to put on the axis ? title ?
#save the file under the name name_vfig in a folder vfig
def plotVeigen(name):
	T = readIterationint(name)
	well_energies = T[0]
	Ef = T[1]
	n2d = T[2]
	z = T[3]
	
	V = T[4]
	#utile ?
	Vnew = T[5]
	intenerg = T[6]
	zrev = [10**9*z[len(z)-1-i] for i in range(len(z))]
	zred = [10**9*z[i] for i in range(len(z))]
	for i in range(2):
		for j in range(len(well_energies[i])):
			if well_energies[i][j] < Ef:
				#psi = renorm(eigenstate(z,well_energies[i][j],V,mz[i]))
				psi = ith_eigenfuncacc(well_energies[i][j],mz[i],V,z)
				psi = [well_energies[i][j]+0*intenerg[i][j] + 0.1*p for p in psi]
				if i == 0:
					plt.plot(zred,psi,label = well_energies[i][j],color = "red")
				if i == 1:
					plt.plot(zred,psi,label = well_energies[i][j],color = "blue")
	plt.plot(zred,V,color = "Black")

	plt.xlabel("z (nm)")
	plt.ylabel("E (eV)")
	#xlim(0,50)
	#ylim(0,0.4)
	plt.title("Eigenstates of the well n = "+str(n2d*10**(-4))+"cm-2")
	#plt.savefig('Well_eigen.pdf')
	plt.show()
	
def plotDensity(name):
	T = readIterationint(name)
	well_energies = T[0]
	Ef = T[1]
	n2d = T[2]
	z = T[3]
	zrev = [10**9*z[i] for i in range(len(z))]
	V = T[4]
	intenerg = T[6]
	
	P = rho3d_multiband_int(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d,intenerg)
	
	Vtild = newpotential(P,n2d,z)
	
	plt.plot(zrev,P,color = "Black")

	plt.xlabel("z (nm)")
	plt.ylabel("LDOS ($m^{-2}$)")
	xlim(0,50)
	plt.title("Density near the interface")
	plt.savefig('Density.pdf')
	plt.show()
	
def computeError(name):
	T = readIterationint(name)
	Vanc = T[4]
	Vnew = T[5]
	return potential_error(Vanc,Vnew)

def plotEnerg(file_number):
	for i in range(file_number+1):
		T = readIterationint(str(i))
		well_energies = T[0]
		intenerg = T[6]
		totalenerg = sum_tab_energies(well_energies,intenerg)
		Ef = T[1]
		tabtab = [i for j in range(len(well_energies[0]))]
		tabEf = [Ef for j in range(len(well_energies[0]))]
		toto = [totalenerg[0][k] - Ef for k in range(len(totalenerg[0]))]
		plt.plot(tabtab,toto,color = 'red', linestyle="None",marker="x")
		#plt.plot(tabtab,tabEf,color = 'black', linestyle = "None",marker = "+")
		toto = [totalenerg[1][k] - Ef for k in range(len(totalenerg[0]))]
		plt.plot(tabtab,toto,color = 'blue', linestyle="None",marker="x")
	plt.show()
	return 1

#to plot the band structure of a given file in the case of interaction
def plotBandsFile(file_number):
	T = readIterationint(str(file_number))
	totalenerg = sum_tab_energies(T[0],T[6])
	Ef = T[1]
	plot_bands(totalenerg,band_energies,Ef,mx,my)

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/IResults16/")

i = 19

for i in range(i+1):
	print(computeError(str(i)))


plotBandsFile(i)
#plotEnerg(i)
plotVeigen(str(i))
plotDensity(str(i))