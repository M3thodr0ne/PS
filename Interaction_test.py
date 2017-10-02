import os

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src/")

execfile("Interaction.py")
execfile("Files.py")

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
					tabE[N+k] = well_energies[i][j]+band_energies[i] - Ef + (h*kvectab[N+k-1])**2/(2*my[i]*el)
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
	title("Y - $\Gamma$ - X band structure")
	xlabel("k ($\pi/a$)")
	ylabel("E (eV)")
	#plt.savefig('Band_struct.pdf')
	plt.show()
	return 1

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/Results16/")	
	
U = 1
U1 = 0

T = readIteration(str(14))
well_energies = T[0]
Ef = T[1]
ntotal = T[2]

BD = band_densities(well_energies,band_energies,mx,my,Ef)

print(sum(BD[1]))

Nstep = 100

for i in range(Nstep):
	TTT = interaction_loop(well_energies,band_energies,mx,my,ntotal,Ef,U,U1,i+1)

	Ef2 = TTT[1]

	BD2 = TTT[0]
	
	intEnergies = TTT[2]
	
	print(BD2[0][1])

#print(intEnergies)

energ = sum_tab_energies(well_energies,intEnergies)

#print(BD)
#print(BD2)

plot_bands(well_energies,band_energies,Ef,mx,my)
plot_bands(energ,band_energies,Ef2,mx,my)