execfile("DOS_alter.py")


#computes the eigenvalues of the problem
def band_struct_Inter(well_energ,param,k,a,U,Mpopbands):
	M = matrixBandsInteract(well_energ,param,k,a,U,Mpopbands)
	lambd = eig(M)
	return lambd[0]

def Ucorr(level):
	if level == 0:
		return 1
	return 0.7

#construct the matrix of the Hamiltonian,param = [U,Mpopbands]
def matrixBandsInteract(well_energ,param,k,a,U,Mpopbands):
	kx = k[0]*a
	ky = k[1]*a
	t1 = 0.05
	t2 = 0.07
	t3 = 1.6
	trigfield = -0.02/3.
	pot = 0.

	t2nxy = - 2 * t1 * cos(  sqrt(3)    * kx             )
	t2nyz = - 2 * t1 * cos( -sqrt(3)/2. * kx + 3./2. * ky)
	t2nxz = - 2 * t1 * cos(  sqrt(3)/2. * kx + 3./2. * ky)

	T2N = [t2nxy , t2nyz , t2nxz]

	tnnxy = - t2 \
	        - 2 * t3 * cos(sqrt(3) /2. *kx) * cmath.exp(-1j * 3/2. * ky)
	tnnyz = - t2 * cmath.exp(-1j * (sqrt(3)/2. * kx + 3./2. * ky)) \
	        - t3 * (1 + cmath.exp( 1j * (sqrt(3)/2.*kx - 3./2. * ky ) ) )
	tnnxz = - t2 * cmath.exp( 1j * (sqrt(3)/2. * kx - 3./2. * ky)) \
	        - t3 * (1 + cmath.exp(-1j * (sqrt(3)/2.*kx + 3./2. * ky ) ) )

	TNN = [tnnxy , tnnyz , tnnxz]

	Nbands = len(well_energ[0])
	Norb = len(well_energ)

	M = [[0 for i in range(4*Norb*Nbands)] for j in range(4*Norb*Nbands)]
	HTRIG = [[0 for i in range(4*Norb*Nbands)] for j in range(4*Norb*Nbands)]
	HSOC = [[0 for i in range(4*Norb*Nbands)] for j in range(4*Norb*Nbands)]

	for i in range(len(M)):
		orb_index = i     / (4 * Nbands)
		i_orb     = i     % (4 * Nbands)
		sub_index = i_orb / 4
		Utild = Ucorr(sub_index)*U
		M[i][i] = well_energ[orb_index][sub_index] \
		        + T2N[orb_index] \
		        + pot * (i % 2) \
				+ Utild*sum(Mpopbands)/3.
				#+ 2**(-sub_index)*0.7*U*sum(Mpopbands)*2./3.
		        #+ 2**(-sub_index)*U*sum(Mpopbands[orb_index]) \
				#+ 2**(-sub_index)*0.7*U*sum(Mpopbands[(orb_index+1)%3]+Mpopbands[(orb_index+2)%3])

	for i_orb in range(Norb):
		for i_sub in range(Nbands):
			for i in range(2):
				M[4*Nbands*i_orb + 4*i_sub + 2*i + 1][4*Nbands*i_orb + 4*i_sub + 2*i] =           TNN[i_orb]
				M[4*Nbands*i_orb + 4*i_sub + 2*i][4*Nbands*i_orb + 4*i_sub + 2*i + 1] = conjugate(TNN[i_orb])
	for i in range(4*Nbands):
		HTRIG[i][4*Nbands+i] = 1
		HTRIG[4*Nbands+i][i] = 1
		HTRIG[i][8*Nbands+i] = 1
		HTRIG[8*Nbands+i][i] = 1
		HTRIG[4*Nbands+i][8*Nbands+i] = 1
		HTRIG[8*Nbands+i][4*Nbands+i] = 1

	for i in range(Nbands):
		HSOC[0*Nbands+4*i  ][4*Nbands+4*i+2] =  1
		HSOC[0*Nbands+4*i+1][4*Nbands+4*i+3] =  1
		HSOC[0*Nbands+4*i+2][4*Nbands+4*i  ] = -1
		HSOC[0*Nbands+4*i+3][4*Nbands+4*i+1] = -1

		HSOC[4*Nbands+4*i  ][0*Nbands+4*i+2] = -1
		HSOC[4*Nbands+4*i+1][0*Nbands+4*i+3] = -1
		HSOC[4*Nbands+4*i+2][0*Nbands+4*i  ] =  1
		HSOC[4*Nbands+4*i+3][0*Nbands+4*i+1] =  1

		HSOC[0*Nbands+4*i  ][8*Nbands+4*i+2] = -1j
		HSOC[0*Nbands+4*i+1][8*Nbands+4*i+3] = -1j
		HSOC[0*Nbands+4*i+2][8*Nbands+4*i  ] = -1j
		HSOC[0*Nbands+4*i+3][8*Nbands+4*i+1] = -1j

		HSOC[8*Nbands+4*i  ][0*Nbands+4*i+2] =  1j
		HSOC[8*Nbands+4*i+1][0*Nbands+4*i+3] =  1j
		HSOC[8*Nbands+4*i+2][0*Nbands+4*i  ] =  1j
		HSOC[8*Nbands+4*i+3][0*Nbands+4*i+1] =  1j

			#yz and xz
		HSOC[4*Nbands+4*i  ][8*Nbands+4*i  ] =  1j
		HSOC[4*Nbands+4*i+1][8*Nbands+4*i+1] =  1j
		HSOC[4*Nbands+4*i+2][8*Nbands+4*i+2] = -1j
		HSOC[4*Nbands+4*i+3][8*Nbands+4*i+3] = -1j

		HSOC[8*Nbands+4*i  ][4*Nbands+4*i  ] = -1j
		HSOC[8*Nbands+4*i+1][4*Nbands+4*i+1] = -1j
		HSOC[8*Nbands+4*i+2][4*Nbands+4*i+2] =  1j
		HSOC[8*Nbands+4*i+3][4*Nbands+4*i+3] =  1j

	Mtot = [[M[i][j] + trigfield * HTRIG[i][j] + param/3*HSOC[i][j] for i in range(len(M))] for j in range(len(M))]
	return Mtot



#compute the density for a given fermi energy
def compute_density_Inter(well_energ,param,Ef,size,a,U,Mpopbands):
	ntot = 0
	for i in range(2*size):
		for j in range(2*size):
			k = Kmesh(a,2*size,size-i,size-j)
			B = band_struct_Inter(well_energ,param,k,a,U,Mpopbands)
			for b in B:
				if real(b) < Ef:
					ntot+=volKmesh(a,2*size)
	return ntot

def rho3d_alter(z,Ef,well_energies,V,Mpop):
	#print(V)
	Mpopbands = [[sum(Mpop[j][4*l+i] for i in range(4)) for l in range(len(well_energies[1]))] for j in range(3)]
	ntot = sum(Mpopbands)
	rhoz = [0 for i in range(len(z))]
	el = 1.6*10**(-19)
	a = 3.9*10**(-10)
	for i in range(len(band_energies)):
		for j in range(len(well_energies[i])):
			e_nalpha = well_energies[i][j]
			if Mpop[i][j] > 10**(-5):
				print(well_energies[0][j])
				ksi_nalpha = ith_eigenfuncacc(well_energies[0][j],mz[0],V,z)
					#ksi_nalpha = renorm(eigenstate(z,,V,mz[i]))
				for k in range(len(z)):
					#multiply by el to have energies in the correct unit as they are in eV
					rhoz[k] += Mpopbands[i][j]*abs(ksi_nalpha[k])**2
	#print(rhoz)
	s = sum(rhoz[k] for k in range(len(z)))
	rhoz = [rhoz[k]*ntot/(s*a**2) for k in range(len(z))]
	return rhoz



def fermiEnergyComputerInter(well_energ,param,size,a,acc,ntot,Efmin0,Efmax0,U,Mpopbands):
	#computation of the initial majoration of the Fermi energy, which is taken at k=0.5/sqrt(pi)kmax in the x direction
	Efmin = Efmin0
	Efmax = Efmax0
	while Efmax-Efmin > acc:
		Efmoy = (Efmin+Efmax)/2
		n = compute_density_Inter(well_energ,param,Efmoy,size,a,U,Mpopbands)
		if n < ntot:
			Efmin = Efmoy
		else:
			Efmax = Efmoy
			print(n)
	return Efmoy


def bandPopComputerInter2(well_energies,param,size,a,Ef,U,Mpopbands):
	Nbands = len(well_energies[0])
	P = [[0 for i in range(4*Nbands)] for j in range(len(well_energies))]
	for i in range(size):
		for j in range(size):
			k = Kmesh(a,2*size,(-size+i),(-size+j))
			B = eig(matrixBandsInteract(well_energies,param,k,a,U,Mpopbands))
			eigenval = B[0]
			eigenvect = B[1]
			for index in range(len(B[0])):
				if real(eigenval[index]) < Ef:
					for index2 in range(len(eigenvect)):
						P[index2/(4*Nbands)][index2%(4*Nbands)] += volKmesh(a,2*size)*abs(eigenvect[index2,index])**2
	return P

def bandPopComputerInter(well_energies,param,size,a,Ef,U,Mpopbands):
	Nbands = len(well_energies[0])
	P = [[0 for i in range(4*Nbands)] for j in range(len(well_energies))]
	for i in range(2*size):
		for j in range(2*size):
			k = Kmesh(a,2*size,(-size+i),(-size+j))
			B = eig(matrixBandsInteract(well_energies,param,k,a,U,Mpopbands))
			eigenval = B[0]
			eigenvect = B[1]
			for index in range(len(B[0])):
				if real(eigenval[index]) < Ef:
					for index2 in range(len(eigenvect)):
						P[index2/(4*Nbands)][index2%(4*Nbands)] += volKmesh(a,2*size)*abs(eigenvect[index2,index])**2
	return P


def energyTotalComputer(well_energies,param,size,a,Ef):
	totalE = 0
	for i in range(2*size):
		for j in range(2*size):
			k = Kmesh(a,2*size,(-size+i),(-size+j))
