import os

execfile("Poisson.py")
execfile("Schrodinger.py")
execfile("DOS.py")
execfile("Files.py")

def ksiDyn(Vanc,Vnew):
	if potential_error(Vanc,Vnew) < 10**(-1):
		return 0.2
	return 0.4

def autocoherentLoopInit(Nspectrumz,acc,surf_density):
	V_intensity = shooting_poisson(W,surf_density,E0/2)

	V = [V_intensity*z[i] for i in range(N_unitcell*L)]
	Vtild = [V[i] for i in range(len(V))]

	ksi = 0.4

	well_energies = [[1000*i for i in range(Nspectrumz)] for j in range(len(band_energies))]
	energ = 0
	energz = 0
	V = potential_corrected(V,Vtild,ksi)
	for i in range(Nspectrumz):
		#well_energies[0][i] = finding_eigenval_i(z,energz,V,mz[0],i,acc)
		well_energies[0][i] = finding_accurate_eigen(z,energz,V,mz[0],i,acc)
		energz = well_energies[0][i]
		#well_energies[1][i] = finding_eigenval_i(z,energ,V,mz[1],i,acc)
	well_energies[1][0] = finding_accurate_eigen(z,energ,V,mz[1],0,acc)
	energ = well_energies[1][0]
	well_energies[2][0] = well_energies[1][0]

	Ef = findingFermiEnergy2(well_energies,band_energies,mx,my,surf_density)



	P = rho3d_multiband2(z,Ef,well_energies,band_energies,mx,my,mz,V,surf_density)
	plot(zr,P,label=j)

	Vtild = newpotential(P,surf_density,z)

	saveIteration(0,well_energies,Ef,z,V,Vtild,surf_density)
	return 1

#read the file number j
def autocoherentLoopOneStep(Nspectrumz,acc,file_number):
	TT = readIteration(str(file_number))
	well_energies = [[1000 for i in range(Nspectrumz)] for j in range(3)]
	Ef = TT[1]
	surf_density = TT[2]
	z = TT[3]
	V = TT[4]
	Vtild = TT[5]
	ksi = ksiDyn(V,Vtild)
	energ = 0
	energz = 0

	V = potential_corrected(V,Vtild,ksi)
	for i in range(Nspectrumz):
		#well_energies[0][i] = finding_eigenval_i(z,energz,V,mz[0],i,acc)
		well_energies[0][i] = finding_accurate_eigen(z,energz,V,mz[0],i,acc)
		energz = well_energies[0][i]
		#well_energies[1][i] = finding_eigenval_i(z,energ,V,mz[1],i,acc)
	well_energies[1][0] = finding_accurate_eigen(z,0,V,mz[1],0,acc)
	energ = well_energies[1][0]
	well_energies[2][0] = well_energies[1][0]
	Ef = findingFermiEnergy2(well_energies,band_energies,mx,my,surf_density)
	#revenir a rho3Dmultiband si ca marche pas
	P = rho3d_multiband2(z,Ef,well_energies,band_energies,mx,my,mz,V,surf_density)
	plot(zr,P,label=j)
	Vtild = newpotential(P,surf_density,z)
	saveIteration(file_number+1,well_energies,Ef,z,V,Vtild,surf_density)
	return 1

def nextLoop(Nspectrumz,acc):
	k = nextFileNumber()
	print(k)
	surf_density = 1.6*10**18
	if k == 0:
		autocoherentLoopInit(Nspectrumz,acc,surf_density)
	else:
		autocoherentLoopOneStep(Nspectrumz,acc,k-1)
	return 1

chemin = os.getcwd()

makeFold("Results16",chemin)

Nspectrumz = 4
acc = 0.00005
for l in range(20):
	nextLoop(Nspectrumz,acc)
