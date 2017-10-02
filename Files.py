execfile("Header.py")

chemin = "/Users/Pierre/Desktop/Stage_4A/Poisson_Schrodinger/Simulations"



#creer un dossier avec le nom name dans le chemin chemin, et se placer dans le dossier en question
def makeFold(name,chemin):
    os.chdir(chemin)
    chemin2 = chemin+"/"+name
    #ifnot permet de ne pas creer si ca existe deja, obligatoire de le mettre sinon erreur possible
    if not os.path.exists(chemin2):
        os.mkdir(name)
    #on se place dans le nouveau dossier cree
    os.chdir(chemin2)
    
    #sauvegarder le resultat d'une iteration dans un fichier
#Format eigenvalues, Fermi energy, z, Vancien, Vnew
def saveIteration(l,eigen,Ef,z,Vanc,Vnew,n2d):
    #ouvrir en lecture qui ecrase tout, attention, initialisation du fichier
    mon_fichier = open(str(l),"w")
    #on print le tableau d'eigenvalues
    for i in range(len(eigen)):
        for j in range(len(eigen[i])):
            mon_fichier.write(str(eigen[i][j]))
            if j < len(eigen[i])-1:
                mon_fichier.write(",")
        mon_fichier.write("\n")
    mon_fichier.write(str(Ef)+"\n")
    mon_fichier.write(str(n2d)+"\n")
    #separateur entre les eigenenergies
    mon_fichier.write("ZV")
    for i in range(len(z)):
        mon_fichier.write(str(z[i])+","+str(Vanc[i])+","+str(Vnew[i])+"\n")
    mon_fichier.close()
    return 1
    
#save a file in the case of the presence of interaction. Saves first the energies of the well, then the energy of the interaction part, then the Fermi energy in the presence of the interaction part (anyway the Fermi energy can be computed from the eigenergies), the 2d density of the gas, the two interaction constants U and U1, and then the z component and the potential
def saveIterationint(l,eigen,Ef,z,Vanc,Vnew,n2d,intenerg,U,U1):
    #ouvrir en lecture qui ecrase tout, attention, initialisation du fichier
	mon_fichier = open(str(l),"w")
	#on print le tableau d'eigenvalues
	for i in range(len(eigen)):
		for j in range(len(eigen[i])):
			mon_fichier.write(str(eigen[i][j]))
			if j < len(eigen[i])-1:
				mon_fichier.write(",")
		mon_fichier.write("\n")
	for i in range(len(intenerg)):
		for j in range(len(intenerg[i])):
			mon_fichier.write(str(intenerg[i][j]))
			if j < len(intenerg[i])-1:
				mon_fichier.write(",")
		mon_fichier.write("\n")
	mon_fichier.write(str(Ef)+"\n")
	mon_fichier.write(str(n2d)+"\n")
	mon_fichier.write(str(U)+"\n")
	mon_fichier.write(str(U1)+"\n")
	#separateur entre les eigenenergies
	mon_fichier.write("ZV")
	for i in range(len(z)):
		mon_fichier.write(str(z[i])+","+str(Vanc[i])+","+str(Vnew[i])+"\n")
	mon_fichier.close()
	return 1

#read a file with the name name, and return in a tab the well energies, the Fermi energy, the z tab, Vanc and Vnew, Vanc being the potential used to compute the well energies and Vnew the computed potential 
def readIteration(name):
    mon_fichier = open(name,"r")

    c = mon_fichier.read()

    t = c.split("ZV")
    
    
    eigenval = t[0].split("\n")
    #counting the number of energy levels
    N = len(eigenval[0].split(","))
    
    well_energies = [[1000*i for i in range(N)] for i in range(3)]
    
    for i in range(3):
        band_eigenval = eigenval[i].split(",")
        for j in range(len(band_eigenval)):
            well_energies[i][j] = float(band_eigenval[j])
    
    Ef = float(eigenval[3])
    n2d = float(eigenval[4])
    
    #on recupere la partie z et V, et on enleve une ligne (la derniere ligne est vide)
    potz = t[1].split("\n")
    L = len(potz)-1
    
    z = [0 for i in range(L)]
    Vanc = [0 for i in range(L)]
    Vnew = [0 for i in range(L)]
    
    for i in range(L):
        s = potz[i].split(",")
        z[i] = float(s[0])
        Vanc[i] = float(s[1])
        Vnew[i] = float(s[2])
    
    mon_fichier.close()

    return [well_energies,Ef,n2d,z,Vanc,Vnew]
    
def readIteration111(name):
    mon_fichier = open(name,"r")

    c = mon_fichier.read()

    t = c.split("ZV")
    
    
    eigenval = t[0].split("\n")
    #counting the number of energy levels
    N = len(eigenval[0].split(","))
    
    well_energies = [[1000*i for i in range(N)] for i in range(6)]
    
    for i in range(6):
        band_eigenval = eigenval[i].split(",")
        for j in range(len(band_eigenval)):
            well_energies[i][j] = float(band_eigenval[j])
    
    Ef = float(eigenval[6])
    n2d = float(eigenval[7])
    
    #on recupere la partie z et V, et on enleve une ligne (la derniere ligne est vide)
    potz = t[1].split("\n")
    L = len(potz)-1
    
    z = [0 for i in range(L)]
    Vanc = [0 for i in range(L)]
    Vnew = [0 for i in range(L)]
    
    for i in range(L):
        s = potz[i].split(",")
        z[i] = float(s[0])
        Vanc[i] = float(s[1])
        Vnew[i] = float(s[2])
    
    mon_fichier.close()

    return [well_energies,Ef,n2d,z,Vanc,Vnew]


#read the file in the interacting case, more things to read
def readIterationint(name):
	mon_fichier = open(name,"r")
	c = mon_fichier.read()

	t = c.split("ZV")
	
	eigenval = t[0].split("\n")
    #counting the number of energy levels
	N = len(eigenval[0].split(","))
    
	well_energies = [[1000*i for i in range(N)] for i in range(3)]
	intenerg = [[0 for i in range(N)] for i in range(3)]
    
	for i in range(3):
		band_eigenval = eigenval[i].split(",")
		inteigenval = eigenval[i+3].split(",")
		for j in range(len(band_eigenval)):
			well_energies[i][j] = float(band_eigenval[j])
			intenerg[i][j] = float(inteigenval[j])
    
	Ef = float(eigenval[6])
	n2d = float(eigenval[7])
	U = float(eigenval[8])
	U1 = float(eigenval[9])
    
	#on recupere la partie z et V, et on enleve une ligne (la derniere ligne est vide)	
	potz = t[1].split("\n")
	L = len(potz)-1
	z = [0 for i in range(L)]
	Vanc = [0 for i in range(L)]
	Vnew = [0 for i in range(L)]
    
	for i in range(L):
		s = potz[i].split(",")
		z[i] = float(s[0])
		Vanc[i] = float(s[1])
		Vnew[i] = float(s[2])
    
	mon_fichier.close()

	return [well_energies,Ef,n2d,z,Vanc,Vnew,intenerg,U,U1]
   
#to determine at which number one should continue to produce a file 
def nextFileNumber():
	i=0
	while os.path.isfile(str(i)):
		i+=1
	return i
		
makeFold("Results",chemin)	