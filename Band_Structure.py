import os

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src/")

execfile("Constantes.py")
execfile("Files.py")

chemin = "/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations"

#makeFold("Results",chemin)	


def plot_bands(well_energies,band_energies,Ef,mx,my):
	a = 3.9*10**(-10)
	N = 1000
	kvectab = [-0.5*pi/a + 0.5*pi/(a*N)*i for i in range(2*N)]
	kplot = [kvectab[i]*a/pi for i in range(len(kvectab))]
	tabE = [0 for i in range(len(kvectab))]
	for i in range(len(well_energies)):
		for j in range(len(well_energies[i])):
			if well_energies[i][j] < 10:
				for k in range(N):
					tabE[k] = well_energies[i][j]+band_energies[i] - Ef + (h*kvectab[k])**2/(2*mx[i]*el)
					tabE[N+k] = well_energies[i][j]+band_energies[i] - Ef + (h*kvectab[N+k-1])**2/(4*el)*(1/mx[i] + 1/my[i])
				if i == 0:
					if j == 0:
						plot(kplot,tabE,color = 'r',label = "xy")
					else:
						plot(kplot,tabE,color = 'r')
				if i == 1:
					if j == 0:
						plot(kplot,tabE,color = 'b',label = "xz")
					else:
						plot(kplot,tabE,color = 'b')
				if i == 2:
					if j == 0:
						plot(kplot,tabE,color = 'g',label = "yz")
					else:
						plot(kplot,tabE,color = 'g')
	efall = [0 for i in range(len(kvectab))]
	legend()
	plot(kplot,efall,color = 'black')
	xlim(-0.3,0.3)
	ylim(-0.25,0.05)
	title("X - $\Gamma$ - M band structure",fontsize = 20)
	xlabel("k ($\pi/a$)",fontsize = 20)
	ylabel("E (eV)",fontsize = 20)
	plt.savefig('Band_struct.pdf')
	plt.show()
	return 1
	
def getEnergies(file_number):
	T = readIteration(str(file_number))
	well_energies = T[0]
	Ef = T[1]
	plot_bands(well_energies,band_energies,Ef,mx,my)
	return 1


import matplotlib.pyplot as plt

def band_structure_BZ(well_energies,band_energies,mx,my,Ef):
	a = 3.9*10**(-10)
	N = 50
	tabkX = [0 for i in range(N)]
	tabkY = [0 for i in range(N)]
	tabkY2 = [0 for i in range(N)]
	tabkM = [0 for i in range(N)]
	tabTot = [0 for i in range(4*N)]
	ticke = ['$\Gamma$','X','M','Y','$\Gamma$']
	xtick = [0,N,2*N,3*N,4*N]
	kf = pi/a
	for band in range(len(well_energies)):
		for subband in range(len(well_energies[band])):
			energ = well_energies[band][subband] + band_energies[band]
			for k in range(N):
				tabkX[k] = energy_TB(k*kf/N,0,mx[band],my[band],a) + energ
				tabkM[k] = energy_TB(kf,k*kf/N,mx[band],my[band],a) + energ
				tabkY[k] = energy_TB(kf-k*kf/N,kf,mx[band],my[band],a) + energ
				tabkY2[k] = energy_TB(0,kf-k*kf/N,mx[band],my[band],a) + energ
				tabTot[k] = tabkX[k] - Ef
				tabTot[N+k] = tabkM[k] - Ef
				tabTot[2*N+k] = tabkY[k] - Ef
				tabTot[3*N+k] = tabkY2[k] - Ef
			if band == 0:
				plot(tabTot,color = "red")
			if band == 1:
				plot(tabTot,color = "blue")
			if band == 2:
				plot(tabTot,color = "green")
				
	plt.xticks(xtick, ticke, rotation='horizontal')
	plt.show()
	
def band_structure_BZ2(well_energies,band_energies,mx,my,Ef):
	a = 3.9*10**(-10)
	N = 200
	tabkX = [0 for i in range(N)]
	tabkM = [0 for i in range(N)]
	tabTot = [0 for i in range(2*N)]
	ticke = ['X','X/2','$\Gamma$','M/2','M']
	xtick = [0,N/2,N,3*N/2,2*N]
	kf = pi/a
	for band in range(len(well_energies)):
		for subband in range(len(well_energies[band])):
			energ = well_energies[band][subband] + band_energies[band]
			for k in range(N):
				tabkX[k] = energy_TB(kf-k*kf/N,0,mx[band],my[band],a) + energ
				tabkM[k] = energy_TB(k*kf/N,k*kf/N,mx[band],my[band],a) + energ
				tabTot[k] = tabkX[k] - Ef
				tabTot[N+k] = tabkM[k] - Ef
			if band == 0:
				plot(tabTot,color = "red")
			if band == 1:
				plot(tabTot,color = "blue")
			if band == 2:
				plot(tabTot,color = "green")
				
	plt.xticks(xtick, ticke, rotation='horizontal')
	plt.show()

def energy_2DEG(kx,ky,mx,my):
	return h**2 * (kx**2/(2*mx) + ky**2/(2*my))/(el)

def energy_TB(kx,ky,mx,my,a):
	return (h/a)**2*((1-cos(kx*a))/mx + (1-cos(ky*a))/my)/el

def band_structure(file_number):
	T = readIteration(str(file_number))
	well_energies = T[0]
	Ef = T[1]
	band_structure_BZ2(well_energies,band_energies,mx,my,Ef)

def findK(theta,energy,mx,my,a,Ef):
	kmin = 0
	kf = pi/a
	kmax = kf
	kmil = 0
	
	while abs(kmax-kmin)/kf > 0.0001:
		kmil = (kmax+kmin)/2.
		energy2 = energy_TB(kmil*cos(theta),kmil*sin(theta),mx,my,a)
		if energy+energy2-Ef > 0:
			kmax = kmil
		else:
			kmin = kmil
	return kmil
		
			

def Fermi_surface(file_number):
	a = 3.9*10**(-10)
	kf = pi/a
	
	T = readIteration(str(file_number))
	well_energies = T[0]
	Ef = T[1]
	
	N = 40
	
	tabtheta = [i*pi/N for i in range(N+1)]
	tabKx = [0 for i in range(2*N+1)]
	tabKy = [0 for i in range(2*N+1)]
	

	for band in range(len(well_energies)):
		for subband in range(len(well_energies[band])):
			energy = well_energies[band][subband] + band_energies[band]
			if energy < Ef:
				for angle in range(len(tabtheta)):
				
					kfound = findK(tabtheta[angle],energy,mx[band],my[band],a,Ef)
					tabKx[angle] = kfound*cos(tabtheta[angle])/kf
					tabKy[angle] = kfound*sin(tabtheta[angle])/kf
					tabKx[N+angle] = -kfound*cos(tabtheta[angle])/kf
					tabKy[N+angle] = -kfound*sin(tabtheta[angle])/kf
					tabKx[-1] = tabKx[0]
					tabKy[-1] = tabKy[0]
				if band == 0:
					plot(tabKx,tabKy,color = 'red')
				if band == 1:
					plot(tabKx,tabKy,color = 'blue')
				if band == 2:
					plot(tabKx,tabKy,color = 'green')
	plt.show()

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/Results16/")	

#getEnergies(18)

#band_structure(26)