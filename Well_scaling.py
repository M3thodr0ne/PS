import os

os.chdir("/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations/src/")

execfile("Header.py")
execfile("Constantes.py")
execfile("Poisson.py")

def nth_energy(m,E,n):
	el = 1.6*10**(-19)
	hbar = 1.054*10**(-34)
	return (el*E*hbar)**(2./3.)/((2*m)**(1./3.)) * (3*pi/2. * abs(n-1./4.))**(1./3.)/el

#extension of the nth Airy function in nm
def nth_extension(m,E,n):
	el = 1.6*10**(-19)
	hbar = 1.054*10**(-34)
	return 10**9*(hbar)**(2./3.)/((2*m*el*E)**(1./3.)) * (3*pi/2. * abs(n-1./4.))**(1./3.)
	
def density_energy_converter(n2d,n,mel):
	E = shooting_poisson(W,n2d,E0)
	return nth_energy(mel,E,n)


def densityComputer(n2d,Ef):
	NN = 10
	el = 1.6*10**(-19)
	hbar = 1.054*10**(-34)
	el_field = shooting_poisson(W,n2d,E0)
	#i+1 for not having 0 state
	tabEQwellHeavy = [nth_energy(mz[0],el_field,i+1) + 0.00 for i in range(NN)]
	tabEQwellLight = [nth_energy(mz[1],el_field,i+1) for i in range(NN)]
	res = 0
	for i in range(NN):
		if tabEQwellHeavy[i] < Ef:
			res += sqrt(mx[0]*my[0])/(pi*hbar**2) * el * (Ef - tabEQwellHeavy[i])
		if tabEQwellLight[i] < Ef:
		#yz and xz bands
			res += 2*sqrt(mx[1]*my[1])/(pi*hbar**2) * el * (Ef - tabEQwellLight[i])
	return res

def fermiEnergySolver(n2d):
	acc = n2d/10000.
	Emin = 0
	Emax = 0.5
	nmax = densityComputer(n2d,Emax)
	while nmax < n2d:
		Emin = Emax
		Emax = Emax+0.5
		nmax = densityComputer(n2d,Emax)
	nmil = nmax
	while abs(nmil-n2d) > acc:
		Emil = (Emin+Emax)/2
		nmil = densityComputer(n2d,Emil)
		if nmil > n2d:
			Emax = Emil
		else:
			Emin = Emil
	return Emax

def density_bands(n2d,Ef):
	NN = 10
	el = 1.6*10**(-19)
	hbar = 1.054*10**(-34)
	el_field = shooting_poisson(W,n2d,E0)
	#i+1 for not having 0 state
	tabEQwellHeavy = [nth_energy(mz[0],el_field,i+1) + 0.002 for i in range(NN)]
	tabEQwellLight = [nth_energy(mz[1],el_field,i+1) for i in range(NN)]
	resheavy = 0
	reslight = 0
	for i in range(NN):
		if tabEQwellHeavy[i] < Ef:
			resheavy += sqrt(mx[0]*my[0])/(pi*hbar**2) * el * (Ef - tabEQwellHeavy[i])
		if tabEQwellLight[i] < Ef:
		#yz and xz bands
			reslight += 2*sqrt(mx[1]*my[1])/(pi*hbar**2) * el * (Ef - tabEQwellLight[i])
	return [resheavy,reslight]

tabN = [10**(16+0.1*i) for i in range(40)]
tabextheavy = [nth_extension(mz[1],shooting_poisson(W,tabN[i],E0),1) for i in range(len(tabN))]
tabextlight = [nth_extension(mz[0],shooting_poisson(W,tabN[i],E0),1) for i in range(len(tabN))]

tabEheavy = [0 for i in range(len(tabN))]
tabElight = [0 for i in range(len(tabN))]
tabEf = [0 for i in range(len(tabN))]
dens_light = [0 for i in range(len(tabN))]
dens_heavy = [0 for i in range(len(tabN))]

for i in range(len(tabN)):
	tabEheavy[i] = density_energy_converter(tabN[i],1,mz[0])+0.00
	tabElight[i] = density_energy_converter(tabN[i],1,mz[1])
	tabEf[i] = fermiEnergySolver(tabN[i])
	T = density_bands(tabN[i],tabEf[i])
	dens_light[i] = T[1]
	dens_heavy[i] = T[0]

for i in range(15):
	print(density_energy_converter(1.6*10**18,i+1,0.7*m))

tabNplot = [tabN[i] *10**(-4) for i in range(len(tabN))]

#loglog(tabNplot,tabEheavy,label = "xy")
#loglog(tabNplot,tabElight,label = "xz or yz")
#loglog(tabNplot,tabEf,label = "Fermi energy")
plot(tabNplot,dens_light,label = "Heavy (xz/yz)")
plot(tabNplot,dens_heavy,label = "Light (xy)")
plot(tabNplot,tabN,label = "total")



#loglog(tabNplot,tabextheavy,label = "Heavy")
#loglog(tabNplot,tabextlight,label = "Light")

xlabel("$n_{2d}$ ($cm^{-2}$)",fontsize = 16)
#ylabel("Spreading length (nm)",fontsize = 16)
#title("Spatial extension of the wavefunction of the ground state of the band")

legend(loc=2, borderaxespad=0.)

#plt.savefig('extension.pdf')

plt.show()
