execfile("Schrodinger.py")


#compute the LDOS in the qwell given the Fermi energy, the well energies, the band energies, the effective masses. V is to compute the eigenfunctions
def rho3d_multiband(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d):
	rhoz = [0 for i in range(len(z))]
	el = 1.6*10**(-19)
	for i in range(len(band_energies)):
		for j in range(len(well_energies[i])):
			e_nalpha = well_energies[i][j]+band_energies[i]
			if e_nalpha < Ef:
				ksi_nalpha = ith_eigenfunction(well_energies[i][j],j,mz[i],V,z)
				#ksi_nalpha = renorm(eigenstate(z,,V,mz[i]))
				for k in range(len(z)):
			#multiply by el to have energies in the correct unit as they are in eV
					rhoz[k] += sqrt(mx[i] * my[i])/(pi * h**2) * el*(Ef - e_nalpha) * abs(ksi_nalpha[k])**2
	s = sum(rhoz[k] for k in range(len(z)))
	rhoz = [rhoz[k]*n2d/s for k in range(len(z))]
	return rhoz
	
def rho3d_multiband2(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d):
	rhoz = [0 for i in range(len(z))]
	el = 1.6*10**(-19)
	for i in range(len(band_energies)):
		for j in range(len(well_energies[i])):
			e_nalpha = well_energies[i][j]+band_energies[i]
			if e_nalpha < Ef:
				ksi_nalpha = ith_eigenfuncacc(well_energies[i][j],mz[i],V,z)
				#ksi_nalpha = renorm(eigenstate(z,,V,mz[i]))
				for k in range(len(z)):
			#multiply by el to have energies in the correct unit as they are in eV
					rhoz[k] += sqrt(mx[i] * my[i])/(pi * h**2) * el*(Ef - e_nalpha) * abs(ksi_nalpha[k])**2
	s = sum(rhoz[k] for k in range(len(z)))
	rhoz = [rhoz[k]*n2d/s for k in range(len(z))]
	return rhoz
	
def rho3d_band(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d,band_number):
	rhoz = [0 for i in range(len(z))]
	el = 1.6*10**(-19)
	for j in range(len(well_energies[band_number])):
		e_nalpha = well_energies[band_number][j]+band_energies[band_number]
		if e_nalpha < Ef:
			ksi_nalpha = ith_eigenfunction(well_energies[band_number][j],j,mz[band_number],V,z)
				#ksi_nalpha = renorm(eigenstate(z,,V,mz[i]))
			for k in range(len(z)):
			#multiply by el to have energies in the correct unit as they are in eV
				rhoz[k] += sqrt(mx[band_number] * my[band_number])/(pi * h**2) * el*(Ef - e_nalpha) * abs(ksi_nalpha[k])**2
	s = sum(rhoz[k] for k in range(len(z)))
	#rhoz = [rhoz[k]*n2d/s for k in range(len(z))]
	return rhoz
	
def rho3d_multiband2(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d):
	rhoz = [0 for i in range(len(z))]
	el = 1.6*10**(-19)
	for i in range(len(band_energies)):
		for j in range(len(well_energies[i])):
			e_nalpha = well_energies[i][j]+band_energies[i]
			if e_nalpha < Ef:
				ksi_nalpha = ith_eigenfuncacc(well_energies[i][j],mz[i],V,z)
				#ksi_nalpha = renorm(eigenstate(z,,V,mz[i]))
				for k in range(len(z)):
			#multiply by el to have energies in the correct unit as they are in eV
					rhoz[k] += sqrt(mx[i] * my[i])/(pi * h**2) * el*(Ef - e_nalpha) * abs(ksi_nalpha[k])**2
	s = sum(rhoz[k] for k in range(len(z)))
	rhoz = [rhoz[k]*n2d/s for k in range(len(z))]
	return rhoz
	
#compute the LDOS in the qwell given the Fermi energy, the well energies, the band energies, the effective masses. V is to compute the eigenfunctions
def rho3d_multiband_int(z,Ef,well_energies,band_energies,mx,my,mz,V,n2d,intenerg):
	rhoz = [0 for i in range(len(z))]
	el = 1.6*10**(-19)
	for i in range(len(band_energies)):
		for j in range(len(well_energies[i])):
			e_nalpha = well_energies[i][j]+band_energies[i]
			if e_nalpha + intenerg[i][j] < Ef:
				ksi_nalpha = ith_eigenfuncacc(well_energies[i][j],mz[i],V,z)
				#ksi_nalpha = renorm(eigenstate(z,,V,mz[i]))
				for k in range(len(z)):
			#multiply by el to have energies in the correct unit as they are in eV
					rhoz[k] += sqrt(mx[i] * my[i])/(pi * h**2) * el*(Ef - e_nalpha - intenerg[i][j]) * abs(ksi_nalpha[k])**2
	s = sum(rhoz[k] for k in range(len(z)))
	rhoz = [rhoz[k]*n2d/s for k in range(len(z))]
	return rhoz

#compute the 2D density given the Fermi energy Ef, the bottom of the bands and the effective masses            
def density2d(Ef,well_energies,mxt,myt):
	el = 1.6*10**(-19)
	n2 = 0
	for well_state in range(len(well_energies)):
		e_nalpha = well_energies[well_state]
		if e_nalpha < Ef:
			n2 += sqrt(mxt * myt) / (pi * h**2) * el * (Ef - e_nalpha)
	return n2

#density of carriers of a band whose bottom at energy, with effective masses mx and my and filled up to Ef
def dens_single_band(energy,Ef,mx,my):
	if energy < Ef:
		return sqrt(mx*my)/(pi*h**2) * el*(Ef-energy)
	return 0

def dens_multi_band(well_energies,band_energies,Ef,mx,my):
	Ndens = 0
	for i in range(len(well_energies)):
		for j in range(len(well_energies[0])):
			energ = well_energies[i][j] + band_energies[i]
			if energ < Ef:
				Ndens += dens_single_band(energ,Ef,mx[i],my[i])
	return Ndens

#compute the 2D density with multiple bands up to the energy E
def ddens_multib(well_energies,band_energies,E,mx,my):
	Nguess = 0
	for i in range(len(well_energies)):
		#bottom of the band
		we = [well_energies[i][j]+band_energies[i] for j in range(len(well_energies[i]))]
		#density of a 2DEG
		Nguess += density2d(E,we,mx[i],my[i])
	return Nguess

#find the Fermi energy by a dichotomy method
def findingFermiEnergy2(well_energies,band_energies,mx,my,n):
	acc = n/10000.
	Emin = 0
	Emax = 0.5
	nmax = dens_multi_band(well_energies,band_energies,Emax,mx,my)
	while nmax < n:
		Emin = Emax
		Emax = Emax+0.5
		nmax = dens_multi_band(well_energies,band_energies,Emax,mx,my)
	nmil = nmax
	while abs(nmil-n) > acc:
		Emil = (Emin+Emax)/2
		nmil = dens_multi_band(well_energies,band_energies,Emil,mx,my)
		if nmil > n:
			Emax = Emil
		else:
			Emin = Emil
	return Emil