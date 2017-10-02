import os
chemin = "/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src"
os.chdir(chemin)
execfile("DOS.py")


def Kronecker(p,q):
	if p == q:
		return 1
	return 0

def interaction_energies(well_densities,Hubbard,Kanamori):
	ntot = sum(well_densities)
	wd = [[well_densities[i][j]/ntot for j in range(len(well_densities[i]))] for i in range(len(well_densities))]
	intEnergies = [[0 for i in range(len(well_densities[0]))] for j in range(len(well_densities))]
	for i in range(len(well_densities)):
		for j in range(len(well_densities[0])):
		#divide by 4 for spin degeneracy of each band for Hubbard, and by 2 for Kanamori because 
			intEnergies[i][j] += Hubbard * sum(wd[i][j] for j in range(len(well_densities[i])))
	return intEnergies
	
def band_densities(well_energies,band_energies,mx,my,Ef):
	densities = [[dens_single_band(well_energies[i][j]+band_energies[i],Ef,mx[i],my[i]) for j in range(len(well_energies[0]))] for i in range(len(well_energies))]
	return densities

#sommer les energies du puits avec l'energie d'interaction
def sum_tab_energies(well_energies,int_energies):
	return [[well_energies[i][j] + int_energies[i][j] for j in range(len(well_energies[0]))] for i in range(len(well_energies))]

#renvoie les energies d'interaction associees a l'occupation actuelles des bandes, ainsi que la nouvelle energie de Fermi.
#Basiquement il s'agit de la brique de base de la boucle d'autocoherence des interactions
def next_energies(well_energies,int_energies,band_energies,mx,my,ntotal,Ef,U,U1):
	BD = band_densities(sum_tab_energies(well_energies,int_energies),band_energies,mx,my,Ef)
	for i in range(len(BD)):
		for j in range(len(BD[i])):
			BD[i][j]/=ntotal
	tabIntE = interaction_energies(BD,U,U1)
	well_energies_corrected = sum_tab_energies(well_energies,tabIntE)
	Ef2 = findingFermiEnergy2(well_energies_corrected,band_energies,mx,my,ntotal)
	return [tabIntE,Ef2]
	
def next_densities(well_energies,int_energies,band_energies,densities,mx,my,ntotal,Ef,U,U1):
	BD = band_densities(sum_tab_energies(well_energies,int_energies),band_energies,mx,my,Ef)
	for i in range(len(BD)):
		for j in range(len(BD[i])):
			BD[i][j] = 0.05*BD[i][j]+0.95*densities[i][j]
	tabIntE = interaction_energies(BD,U,U1)
	well_energies_corrected = sum_tab_energies(well_energies,tabIntE)
	Ef2 = findingFermiEnergy2(well_energies_corrected,band_energies,mx,my,ntotal)
	return [BD,Ef2]
	
def interaction_loop(well_energies,band_energies,mx,my,ntotal,Ef,U,U1,Nstep):
	BD = band_densities(well_energies,band_energies,mx,my,Ef)
	Ef2 = Ef
	for i in range(Nstep):
		int_energies = interaction_energies(BD,U,U1)
		MM = next_densities(well_energies,int_energies,band_energies,BD,mx,my,ntotal,Ef2,U,U1)
		BD = MM[0]
		Ef2 = MM[1]
		int_energies = interaction_energies(BD,U,U1)
	return [BD,Ef2,int_energies]