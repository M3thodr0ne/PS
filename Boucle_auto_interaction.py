import os
chemin = "/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src"
os.chdir(chemin)
execfile("Poisson.py")
execfile("Schrodinger.py")
execfile("DOS.py")
execfile("Interaction.py")
execfile("Files.py")

def ksiDyn(Vanc,Vnew):
	if potential_error(Vanc,Vnew) < 10**(-1):
		return 0.2
	return 0.4

def autocoherentLoopInitInt(Nspectrumz,acc,surf_density,U,U1):
	V_intensity = shooting_poisson(W,surf_density,E0/2)

	V = [V_intensity*z[i] for i in range(N_unitcell*L)]
	Vtild = [V[i] for i in range(len(V))]

	ksi = 0.4

	well_energies = [[1000*i for i in range(Nspectrumz)] for j in range(len(band_energies))]
	energ = 0
	energz = 0
	V = potential_corrected(V,Vtild,ksi)
	for i in range(Nspectrumz):
		well_energies[0][i] = finding_eigenval_i(z,energz,V,mz[0],i,acc)
		energz = well_energies[0][i]
		well_energies[1][i] = finding_eigenval_i(z,energ,V,mz[1],i,acc)
		energ = well_energies[1][i]
		well_energies[2][i] = well_energies[1][i]
    
	Ef = findingFermiEnergy2(well_energies,band_energies,mx,my,surf_density)
	
	Tab2 = interaction_loop(well_energies,band_energies,mx,my,surf_density,Ef,U,U1,40)
	
	Ef2 = Tab2[1]
	intenerg = Tab2[2]
	
	total_energy = sum_tab_energies(well_energies,intenerg)
	
	P = rho3d_multiband_int(z,Ef2,well_energies,band_energies,mx,my,mz,V,surf_density,intenerg)
	#plot(zr,P,label=j)
    
	Vtild = newpotential(P,surf_density,z)
    
	saveIterationint(0,well_energies,Ef2,z,V,Vtild,surf_density,intenerg,U,U1)
    
	return 1
    
#read the file number j
def autocoherentLoopOneStepint(Nspectrumz,acc,file_number,U,U1):
	TT = readIterationint(str(file_number))
	well_energies = [[0 for i in range(Nspectrumz)] for j in range(3)]
	Ef = TT[1]
	surf_density = TT[2]
	z = TT[3]
	V = TT[4]
	Vtild = TT[5]
	U = TT[7]
	U1 = TT[8]
	ksi = ksiDyn(V,Vtild)
	energ = 0
	energz = 0
	
	V = potential_corrected(V,Vtild,ksi)
	for i in range(Nspectrumz):
		well_energies[0][i] = finding_accurate_eigen(z,energz,V,mz[0],i,acc)
		energz = well_energies[0][i]
		well_energies[1][i] = finding_accurate_eigen(z,energ,V,mz[1],i,acc)
		energ = well_energies[1][i]
		well_energies[2][i] = well_energies[1][i]
	Ef = findingFermiEnergy2(well_energies,band_energies,mx,my,surf_density)
	
	Tab2 = interaction_loop(well_energies,band_energies,mx,my,surf_density,Ef,U,U1,40)
	
	Ef2 = Tab2[1]
	intenerg = Tab2[2]
	
	P = rho3d_multiband_int(z,Ef2,well_energies,band_energies,mx,my,mz,V,surf_density,intenerg)
	
	#plot(zr,P,label=j)
	Vtild = newpotential(P,surf_density,z)
	saveIterationint(file_number+1,well_energies,Ef2,z,V,Vtild,surf_density,intenerg,U,U1)
	return 1
 
def nextLoop(Nspectrumz,acc):
	k = nextFileNumber()
	U = 2
	U1 = 0
	print(k)
	surf_density = 1.2*10**18
	if k == 0:
		autocoherentLoopInitInt(Nspectrumz,acc,surf_density,U,U1)
	else:
		autocoherentLoopOneStepint(Nspectrumz,acc,k-1,U,U1)
	return 1
 

makeFold("IResults12",chemin)
   
Nspectrumz = 8
acc = 0.0001
for l in range(10):
	nextLoop(Nspectrumz,acc)