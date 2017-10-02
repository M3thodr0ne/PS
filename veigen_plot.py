import os

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src/")

execfile("Schrodinger.py")
execfile("Poisson.py")
execfile("DOS.py")
execfile("Files.py")


#to plot the potential and the eigenfunctions below the Fermi energy of a file stored under the name name
#what legend to put on the axis ? title ?
#save the file under the name name_vfig in a folder vfig
def plotVeigen2(name):
	T = readIteration(name)
	well_energies = T[0]
	Ef = T[1]
	n2d = T[2]
	z = T[3]
	
	V = T[4]
	#utile ?
	Vnew = T[5]
	zrev = [10**9*z[len(z)-1-i] for i in range(len(z))]
	zred = [10**9*z[i] for i in range(len(z))]
	for i in range(2):
		for j in range(len(well_energies[i])):
			if well_energies[i][j] < Ef:
				#psi = renorm(eigenstate(z,well_energies[i][j],V,mz[i]))
				psi = ith_eigenfuncacc(well_energies[i][j],mz[i],V,z)
				psi = [well_energies[i][j] + 0.1*p for p in psi]
				if i == 0:
					plt.plot(zred,psi,label = well_energies[i][j],color = "red")
				if i == 1:
					plt.plot(zred,psi,label = well_energies[i][j],color = "blue")
	plt.plot(zred,V,color = "Black")

	plt.xlabel("z (nm)")
	plt.ylabel("E (eV)")
	xlim(0,50)
	ylim(0.,0.45)
	plt.title("Eigenstates of the well n = "+str(n2d*10**(-4))+"cm-2")
	plt.savefig('Well_eigen_lim.pdf')
	plt.show()
	
def plotVeigen(name):
	T = readIteration(name)
	well_energies = T[0]
	Ef = T[1]
	n2d = T[2]
	z = T[3]
	
	V = T[4]
	#utile ?
	Vnew = T[5]
	zrev = [10**9*z[len(z)-1-i] for i in range(len(z))]
	zred = [10**9*z[i] for i in range(len(z))]
	for i in range(2):
		for j in range(len(well_energies[i])):
			if well_energies[i][j] < Ef:
				#psi = renorm(eigenstate(z,well_energies[i][j],V,mz[i]))
				psi = ith_eigenfuncacc(well_energies[i][j],mz[i],V,z)
				psi = [well_energies[i][j] + p**2 for p in psi]
				if i == 0:
					plt.plot(zred,psi,label = well_energies[i][j],color = "red")
				if i == 1:
					plt.plot(zred,psi,label = well_energies[i][j],color = "blue")
	plt.plot(zred,V,color = "Black")

	plt.xlabel("z (nm)")
	plt.ylabel("E (eV)")
	xlim(0,50)
	ylim(0.,0.45)
	plt.title("Eigenstates of the well n = "+str(n2d*10**(-4))+"cm-2")
	plt.savefig('Well_eigen_lim.pdf')
	plt.show()

#plot the total density of the file with the name name	
def plotDensity(name):
	T = readIteration(name)
	well_energies = T[0]
	Ef = T[1]
	n2d = T[2]
	z = T[3]
	zrev = [10**9*z[i] for i in range(len(z))]
	V = T[4]
	
	P = rho3d_multiband(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d)
	
	Vtild = newpotential(P,n2d,z)
	
	stepzzz = z[1]-z[0]
	
	Pscaled = [p/stepzzz * 10**(-27) for p in P]
	
	plt.plot(zrev,Pscaled,color = "Black")

	plt.xlabel("z (nm)",fontsize = 20)
	plt.ylabel("Density ($nm^{-3}$)",fontsize = 20)
	xlim(0,50)
	#ylim(0,1)
	plt.title("Density near the interface",fontsize = 16)
	xlim(0,50)
	#ylim(0,1)
	
	plt.savefig('Density.pdf')
	plt.show()

#plot the density of only the light band
def plotdxyDensityMulti(name):
	T = readIteration(name)
	well_energies = T[0]
	Ef = T[1]
	n2d = T[2]
	z = T[3]
	zrev = [10**9*z[i] for i in range(len(z))]
	V = T[4]
	
	P = rho3d_band(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d,0)
	
	print(sum(z[i]*P[i] for i in range(len(z)))/sum(P))
	
	Vtild = newpotential(P,n2d,z)
	
	stepzzz = z[1]-z[0]
	
	Pscaled = [p/stepzzz * 10**(-27) for p in P]
	
	plt.plot(zrev,Pscaled,label = name)

	plt.xlabel("z (nm)")
	plt.ylabel("Density ($nm^{-3}$)")
	xlim(0,50)
	#ylim(0,1)
	plt.title("Density of light bands near the interface")
	
#plot the density of only the heavy band
def plotdxzDensityMulti(name):
	T = readIteration(name)
	well_energies = T[0]
	Ef = T[1]
	n2d = T[2]
	z = T[3]
	zrev = [10**9*z[i] for i in range(len(z))]
	V = T[4]
	
	P = rho3d_band(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d,1)
	
	print(sum(z[i]*P[i] for i in range(len(z)))/sum(P))
	
	Vtild = newpotential(P,n2d,z)
	
	stepzzz = z[1]-z[0]
	
	Pscaled = [p/stepzzz * 10**(-27) for p in P]
	
	plt.plot(zrev,Pscaled,label = name)

	plt.xlabel("z (nm)")
	plt.ylabel("Density ($nm^{-3}$)")
	xlim(0,50)
	#ylim(0,1)
	plt.title("Density of heavy bands near the interface")

def plotdxzMulti(state_number,name):
	T = readIteration(name)
	well_energies = T[0]
	Ef = T[1]
	n2d = T[2]
	z = T[3]
	
	a = 3.9*10**(-10)
	
	V = T[4]
	#utile ?
	Vnew = T[5]
	zred = [10**9*z[i] for i in range(len(z))]
	
	psi = ith_eigenfuncacc(well_energies[1][state_number],mz[1],V,z)
	psicar = [abs(ps)**2 for ps in psi]
	plt.plot(zred,psicar,label = str(name))
	plt.xlabel("z (nm)")
	plt.ylabel('$|\psi|^2$')
	#xlim(0,50)
	#ylim(0,0.4)
	plt.title("Eigenfunctions of the well, n2d = "+str(round(a**2*n2d,3))+" e per unit cell")

def plotdxyMulti(state_number,name):
	T = readIteration(name)
	well_energies = T[0]
	Ef = T[1]
	n2d = T[2]
	z = T[3]
	
	a = 3.9*10**(-10)

	
	V = T[4]
	#utile ?
	Vnew = T[5]
	zred = [10**9*z[i] for i in range(len(z))]
	
	psi = ith_eigenfuncacc(well_energies[0][state_number],mz[0],V,z)
	psicar = [abs(ps)**2 for ps in psi]
	plt.plot(zred,psicar,label = "xy"+str(state_number))
	plt.xlabel("z (nm)")
	plt.ylabel('$|\psi|^2$')
	#xlim(0,50)
	#ylim(0,0.4)
	plt.title("Eigenfunctions of the well, n2d = "+str(round(a**2*n2d,3))+" e per unit cell")

	
def plotDensityMulti(name):
	T = readIteration(name)
	well_energies = T[0]
	Ef = T[1]
	n2d = T[2]
	z = T[3]
	zrev = [10**9*z[i] for i in range(len(z))]
	V = T[4]
	
	P = rho3d_multiband(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d)
	
	Vtild = newpotential(P,n2d,z)
	
	stepzzz = z[1]-z[0]
	
	Pscaled = [p/stepzzz * 10**(-27)/n2d for p in P]
	
	plt.plot(zrev,Pscaled,label = str(n2d*(3.9*10**(-10))**2))

	plt.xlabel("z (nm)",fontsize = 20)
	plt.ylabel("Density ($nm^{-3}$)",fontsize = 20)
	xlim(0,50)
	#ylim(0,1)
	plt.title("Density near the interface",fontsize = 16)
	
def plotElectronicDensity(name):
	T = readIteration(name)
	well_energies = T[0]
	Ef = T[1]
	n2d = T[2]
	z = T[3]
	V = T[4]
	
	P = rho3d_multiband(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d)
	
	
	a = 3.9*10**(-10)
	
	cell_size = int(a/stepz)
	
	N_cell = int(len(z)/cell_size)
	
	zred = [i for i in range(N_cell)]
	
	P2 = [sum(P[cell+cell_size*ncell]*a**2 for cell in range(cell_size)) for ncell in range(N_cell)]
	
	print(sum(z[i]*P[i] for i in range(len(z)))/sum(P))
	
	plt.plot(zred,P2,color = "Black")
	
	plt.xlabel("z (UC)")
	plt.ylabel("Density (electron per layer)")
	xlim(0,20)
	#ylim(10**(-4),1)
	plt.title("Density near the interface")
	plt.savefig('DensityE.pdf')
	plt.show()
	
def computeError(name):
	T = readIteration(name)
	Vanc = T[4]
	Vnew = T[5]
	return potential_error(Vanc,Vnew)

def plotEnerg(file_number):
	for i in range(file_number+1):
		T = readIteration(str(i))
		well_energies = T[0]
		Ef = T[1]
		tabtab = [i for j in range(len(well_energies[0]))]
		tabEf = [Ef for j in range(len(well_energies[0]))]
		plt.plot(tabtab,[well_energies[0][k]-Ef for k in range(len(well_energies[0]))],color = 'red', linestyle="None",marker="x")
		#plt.plot(tabtab,tabEf,color = 'black', linestyle = "None",marker = "+")
		plt.plot(tabtab,[well_energies[1][k]-Ef for k in range(len(well_energies[0]))],color = 'blue', linestyle="None",marker="x")
	plt.show()
	return 1

def plotVanc(file_number):
	
	T = readIteration(str(file_number))
	Ef = T[1]
	z = T[3]
	zrev = [10**9*z[i] for i in range(len(z))]
	zef = [0.01*10**9*z[-1]*i for i in range(100)]
	TEF = [Ef for i in range(len(zef))]
	V = T[4]
	V2 = T[5]
	plt.plot(zrev,V)
	plt.plot(zrev,V2)
	plt.plot(zef,TEF,'--',color = "black")
	plt.show()
	return 1

def plotepsilon(file_number):
	T = readIteration(str(file_number))
	well_energies = T[0]
	Ef = T[1]
	n2d = T[2]
	z = T[3]
	
	V = T[4]
	
	zred = [zz * 10**9 for zz in z]
	
	P = rho3d_multiband(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d)
	
	Pscaled = [p/stepz * 10**(-27) for p in P]
	
	tabE = [0 for i in range(len(z))]
	
	for i in range(len(z)-1):
		tabE[i] = (V[i+1]-V[i])/stepz
		
	print(sum(tabE[i]*P[i] for i in range(len(tabE)))/sum(P ))
	print(sum(z[i]*P[i] for i in range(len(z)))/sum(P ))
	print(sum(z[i]**2*P[i] for i in range(len(z)))/sum(P ))
	
	a = 3.9*10**(-10)
	

	fig, ax1 = plt.subplots()

	ax2 = ax1.twinx()
	ax1.plot(zred, Pscaled, 'black')
	ax2.plot(zred, tabE, 'b-')

	ax1.set_xlabel('z (nm)')
	ax1.set_ylabel("Density ($nm^{-3}$)", color='black')
	ax2.set_ylabel("E ($V.m^{-1}$)", color='b')

	plt.show()
	
def Efinal(file_number):
	T = readIteration(str(file_number))
	Ef = T[1]
	z = T[3]
	V = T[4]
	return [Ef,(V[-1]-V[-2])/(z[1]-z[0])*5*10**(-4)]
	

numberfile = 17

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/Results16111accu/")

tabEf = [0 for i in range(numberfile+1)]
tabV = [0 for i in range(numberfile+1)] 

for j in range(numberfile+1):
	print(computeError(str(j)))
	VVV = Efinal(j)
	tabEf[j] = VVV[0]
	tabV[j] = VVV[1]

#plot(tabV,tabEf)
#plt.show()

#for j in range(8):
	#plotdxyDensityMulti(str(4*j+3))
	#plotdxzMulti(0,str(4*j))

#plotdxyMulti(0,str(numberfile))
#plotdxyMulti(1,str(numberfile))
#plotdxzMulti(0,str(numberfile))
#plotdxyMulti(6,str(numberfile))
#plotdxzMulti(1,str(numberfile))

#plt.legend()
#plt.show()

#legend()
#plt.show()


#plotVanc(i)
plotEnerg(i)
plotVeigen(str(i))

#plotepsilon(i)
plotDensity(str(i))