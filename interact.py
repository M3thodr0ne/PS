import os

execfile("DOS.py")



#sommer les energies du puits avec l'energie d'interaction
def sum_tab_energies(well_energies,int_energies):
	return [[well_energies[i][j] + int_energies[i][j] for j in range(len(well_energies[0]))] for i in range(len(well_energies))]


def band_densities(well_energies,band_energies,mx,my,Ef):
	densities = [[dens_single_band(well_energies[i][j]+band_energies[i],Ef,mx[i],my[i]) for j in range(len(well_energies[0]))] for i in range(len(well_energies))]
	return densities

def interaction_energies(well_densities,U):
	ntot = sum(well_densities)
	wd = [[well_densities[i][j]/ntot for j in range(len(well_densities[i]))] for i in range(len(well_densities))]
	intEnergies = [[0 for i in range(len(well_densities[0]))] for j in range(len(well_densities))]
	for i in range(len(well_densities)):
		for j in range(len(well_densities[0])):
		#divide by 4 for spin degeneracy of each band for Hubbard, and by 2 for Kanamori because
			intEnergies[i][j] += U/2 * sum(wd[i][j] for j in range(len(well_densities[i])))
	return intEnergies



def interaction_loop(well_energies,band_energies,mx,my,ntotal,Ef,U,Nstep):

	int_energies = [[0 for j in range(len(well_energies[0]))] for i in range(len(well_energies))]

	BD = band_densities(well_energies,band_energies,mx,my,Ef)

	Efcorr = Ef

	for i in range(Nstep):
		int_energies = interaction_energies(BD,U)

		total_energies = sum_tab_energies(well_energies,int_energies)

		Efcorr = findingFermiEnergy2(total_energies,band_energies,mx,my,ntotal)

		BD2 = band_densities(total_energies,band_energies,mx,my,Efcorr)

		BD = [[0.99*BD[i][j] + 0.01*BD2[i][j] for j in range(len(well_energies[i]))] for i in range(len(well_energies))]

	return [BD,Efcorr,int_energies]
