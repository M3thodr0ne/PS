execfile("Header.py")

def Kmesh(a,size,x,y):
	deltak = pi/(a*size)
	k1x = deltak*sqrt(3)/3
	k1y = deltak
	return [(x+y)*k1x,(x-y)*k1y]
	#return [deltak*x,deltak*y]

def sampling(a,size):
    tabK = [Kmesh(a,2*size,float(-size+i+0.5),float(-size+j+0.5)) for i in range(2*size) for j in range(2*resol)]
    tabKx = [tabK[i][0] for i in range(len(tabK))]
    tabKy = [tabK[i][1] for i in range(len(tabK))]
    return [tabKx,tabKy]

def hexagmesh(a,size):
    v1 = [0,1]
    v2 = [1,0]
    v3 = [1,1]
    deltak = pi/(a*size)
    k1x = deltak*sqrt(3)/3
    k1y = deltak
    tabKx = [0 for i in range(6*size)]
    tabKy = [0 for i in range(6*size)]
    for i in range(size):
        k = Kmesh(a,2*size,i*v1[0],i*v1[1])
        tabKx[6*i] = k[0]
        tabKy[6*i] = k[1]
        k = Kmesh(a,2*size,-i*v1[0],-i*v1[1])
        tabKx[6*i+1] = k[0]
        tabKy[6*i+1] = k[1]
        k = Kmesh(a,2*size,i*v2[0],i*v2[1])
        tabKx[6*i+2] = k[0]
        tabKy[6*i+2] = k[1]
        k = Kmesh(a,2*size,-i*v2[0],-i*v2[1])
        tabKx[6*i+3] = k[0]
        tabKy[6*i+3] = k[1]
        k = Kmesh(a,2*size,i*v3[0],i*v3[1])
        tabKx[6*i+4] = k[0]
        tabKy[6*i+4] = k[1]
        k = Kmesh(a,2*size,-i*v3[0],-i*v3[1])
        tabKx[6*i+5] = k[0]
        tabKy[6*i+5] = k[1]
    return [tabKx,tabKy]


a = 1.
resol = 5

T = hexagmesh(a/2,resol)

tabKx = T[0]
tabKy = T[1]

plot(tabKx,tabKy,'o')

T = sampling(a,resol)

tabKx = T[0]
tabKy = T[1]

plot(tabKx,tabKy,'x')


plt.show()
