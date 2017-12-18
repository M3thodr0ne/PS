import os

execfile("Poisson.py")
execfile("Schrodinger.py")
execfile("DOS_alter_interact.py")
execfile("interact_alter.py")
execfile("Files.py")

def ksiDyn(Vanc,Vnew):
	if potential_error(Vanc,Vnew) < 10**(-1):
		return 0.2
	return 0.4

def autocoherentLoopInit(Nspectrumz,acc,surf_density,param):
	a = 3.9*10**(-10)
	V_intensity = shooting_poisson(W,surf_density/a**2,E0/2)
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
		well_energies[1][i] = energz
		well_energies[2][i] = energz

	Ef = fermiEnergyDefault(well_energies,param,surf_density)
	Mpop = defaultDensity(well_energies,param,Ef)

	U = 3
	size = 15
	Nstep = 20
	T = interaction_loop_alter(well_energies,surf_density,Ef,U,Nstep,param,size,a,acc)
	Efcorr = T[1]
	Mpopbands = T[0]

	savePop(0,Mpopbands)
	P = rho3d_alter(z,Efcorr,well_energies,V,Mpopbands)
	Vtild = newpotential(P,surf_density/a**2,z)
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
	param = 0
	a = 3.9*10**(-10)
	V = potential_corrected(V,Vtild,ksi)
	for i in range(Nspectrumz):
		#well_energies[0][i] = finding_eigenval_i(z,energz,V,mz[0],i,acc)
		well_energies[0][i] = finding_accurate_eigen(z,energz,V,mz[0],i,acc)
		energz = well_energies[0][i]
		#well_energies[1][i] = finding_eigenval_i(z,energ,V,mz[1],i,acc)
		well_energies[1][i] = energz
		well_energies[2][i] = energz

	Ef = fermiEnergyDefault(well_energies,param,surf_density)
	Mpop = defaultDensity(well_energies,param,Ef)

	U = 3
	size = 15
	Nstep = 20

	T = interaction_loop_alter(well_energies,surf_density,Ef,U,Nstep,param,size,a,acc)

	Efcorr = T[1]
	Mpopbands = T[0]
	savePop(file_number+1,Mpopbands)
	P = rho3d_alter(z,Efcorr,well_energies,V,Mpopbands)
	Vtild = newpotential(P,surf_density/a**2,z)
	saveIteration(file_number+1,well_energies,Ef,z,V,Vtild,surf_density)

	return 1

def nextLoop(Nspectrumz,acc):
    k = nextFileNumber()
    print(k)
    surf_density = 0.25
    param = 0
    if k == 0:
        autocoherentLoopInit(Nspectrumz,acc,surf_density,param)
    else:
        autocoherentLoopOneStep(Nspectrumz,acc,k-1)
    return 1

chemin = os.getcwd()
makeFold("111",chemin)
chemin = os.getcwd()
makeFold("Kmesh15testU3pop2",chemin)

Nspectrumz = 4
acc = 0.0001
for l in range(10):
	nextLoop(Nspectrumz,acc)
