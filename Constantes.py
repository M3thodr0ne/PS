#constants for the 111 orientation

el = 1.6*10**(-19)
epsilon0 = 8.85 * 10**(-12)
B = 25462.
E0 = 82213.
h = 1.054*10**(-34)

#mass of the electron
m = 9.1*10**(-31)

#size 100nm
L = 50
N_unitcell = 160
step = 0.5*10**(-10)
stepz = 0.1*10**(-10)
step2 = 1*10**(-10)
z = [float(i)*step for i in range(N_unitcell*L)]
for i in range(L*N_unitcell/2):
	z[i] = stepz*float(i)
	z[L*N_unitcell/2+i] = stepz*N_unitcell*L/2 + step2*i
z = [float(i)*step for i in range(N_unitcell*L)]
zr = [z[len(z)-i-1] for i in range(len(z))]

#effective masses and band energies in the basis xy, yz, xz
mx = [8*m,8*m,8*m]
my = [0.3*m,0.3*m,0.3*m]
mz = [m,m,m]
band_energies = [0,0,0]