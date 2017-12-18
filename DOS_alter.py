
execfile("Schrodinger.py")

import cmath

#computes the eigenvalues of the problem
def band_struct(well_energ,param,k,a):
	M = matrixBands(well_energ,param,k,a)
	lambd = eig(M)
	return lambd[0]


#construct the matrix of the Hamiltonian,param = [U,Mpopbands]
#matrix per block, 3 big blocks for the orbitals of size 4*Nbands
#inside which Nbands blocks for each subband of size 4
#then site 1 up site 2 up site 1 down site 2 down
def matrixBands(well_energ,param,k,a):
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

		M[i][i] = well_energ[orb_index][sub_index] + T2N[orb_index] + pot * (i % 2)

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

	Mtot = [[M[i][j] + trigfield * HTRIG[i][j] + param/3 * HSOC[i][j] for i in range(len(M))] for j in range(len(M))]
	return Mtot

#band dispersion of free electrons on a 2DEG, in eV
def free_elec(mx,my,a,k):
	return h**2 * (k[0]**2 / mx + k[1]**2 / my)/(2*el)


#element of the Kmesh of coordinates x,y < size
def Kmesh(a,size,x,y):
	deltak = pi/(a*size)
	k1x = deltak*sqrt(3)/3
	k1y = deltak
	return [(x+y)*k1x,x*k1y-y*k1y]
	#return [deltak*x,deltak*y]

#elementary volume of the Kmesh
def volKmesh(a,size):
	return 1. * 1/size**2
	#return 2. * 1/size**2 *sqrt(3)/2

#Fermi Dirac statistics, T in eV
def fermiDirac(Ef,T,E):
	return 1/(exp((E-Ef)/T))



#compute the density for a given fermi energy
def compute_density(well_energ,param,Ef,size,a):
	ntot = 0
	for i in range(2*size):
		for j in range(2*size):
			k = Kmesh(a,2*size,size-i,size-j)
			B = band_struct(well_energ,param,k,a)
			for b in B:
				if real(b) < Ef:
					ntot+=volKmesh(a,2*size)
	return ntot

#values per default for the size of the unit cell and the number of points in the k mesh
def fermiEnergyDefault(well_energ,param,ntot):
	a = 3.9*10**(-10)
	size = 15
	acc = 0.001
	Ef = fermiEnergyComputer(well_energ,param,size,a,acc,ntot)
	return Ef

#return the occupation of the differnet bands and subbands in the default case
def defaultDensity(well_energies,param,Ef):
	a = 3.9*10**(-10)
	size = 15
	acc = 0.001

	P = bandPopComputer(well_energies,param,size,a,Ef)
	return P


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
			ksi_nalpha = ith_eigenfuncacc(well_energies[i][j],mz[i],V,z)
				#ksi_nalpha = renorm(eigenstate(z,,V,mz[i]))
			for k in range(len(z)):
			#multiply by el to have energies in the correct unit as they are in eV
				rhoz[k] += Mpopbands[i][j]*abs(ksi_nalpha[k])**2
	#print(rhoz)
	s = sum(rhoz[k] for k in range(len(z)))
	rhoz = [rhoz[k]*ntot/(s*a**2) for k in range(len(z))]
	return rhoz



#compute the Fermi energy of the system by dichotomy
def fermiEnergyComputer(well_energ,param,size,a,acc,ntot):
	#computation of the initial majoration of the Fermi energy, which is taken at k=0.5/sqrt(pi)kmax in the x direction
	Ef0 = initialFermiEnergy(well_energ,param,size,a)
	Efmin = -Ef0
	Efmax = Ef0
	while Efmax-Efmin > acc:
		Efmoy = (Efmin+Efmax)/2
		n = compute_density(well_energ,param,Efmoy,size,a)
		print(n)
		if n < ntot:
			Efmin = Efmoy
		else:
			Efmax = Efmoy
	return Efmoy

def fermiEnergyComputerTot(well_energ,param,size,a,acc,ntot):
	#computation of the initial majoration of the Fermi energy, which is taken at k=0.5/sqrt(pi)kmax in the x direction
	tabE = []
	tabN = []
	tabAcc = []
	Ef0 = initialFermiEnergy(well_energ,param,size,a)
	Efmin = -Ef0
	Efmax = Ef0
	while Efmax-Efmin > acc:
		Efmoy = (Efmin+Efmax)/2
		n = compute_density(well_energ,param,Efmoy,size,a)
		print(n)
		if n < ntot:
			Efmin = Efmoy
		else:
			Efmax = Efmoy
		tabN.append(n)
		tabE.append(Efmoy)
		tabAcc.append(Efmax-Efmin)
	return [tabAcc,tabE,tabN]



#to compute the initial Fermi energy, we compute the highest energy for k in the diagonal at a certain radius for the bands we have computed
def initialFermiEnergy(well_energ,param,size,a):
	k = Kmesh(a,size,size/(2*sqrt(pi)),size/(2*sqrt(pi)))
	B = band_struct(well_energ,param,k,a)
	return max(max(real(B))+1,max(real(-B))+1)

#result is in electron per unit cell
def bandPopComputer2(well_energies,param,size,a,Ef):
	Nbands = len(well_energies[0])
	P = [[0 for i in range(4*Nbands)] for j in range(len(well_energies))]
	for i in range(size):
		for j in range(size):
			k = Kmesh(a,2*size,(-size+i),(-size+j))
			B = eig(matrixBands(well_energies,param,k,a))
			eigenval = B[0]
			eigenvect = B[1]
			for index in range(len(B[0])):
				if real(eigenval[index]) < Ef:
					for index2 in range(len(eigenvect)):
						P[index2/(4*Nbands)][index2%(4*Nbands)] += volKmesh(a,2*size)*abs(eigenvect[index2,index])**2
	return P

def bandPopComputer(well_energies,param,size,a,Ef):
	Nbands = len(well_energies[0])
	P = [[0 for i in range(4*Nbands)] for j in range(len(well_energies))]
	for i in range(2*size):
		for j in range(2*size):
			k = Kmesh(a,2*size,(-size+i),(-size+j))
			B = eig(matrixBands(well_energies,param,k,a))
			eigenval = B[0]
			eigenvect = B[1]
			for index in range(len(B[0])):
				if real(eigenval[index]) < Ef:
					for index2 in range(len(eigenvect)):
						P[index2/(4*Nbands)][index2%(4*Nbands)] += volKmesh(a,2*size)*abs(eigenvect[index2,index])**2
	return P
