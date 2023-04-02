# Kart over Bølgefunksjon**2 eller sannsynlighet for å finne elektron et bestemt sted i z-x planet
# Etter å ha laget tetthetsplottet simuleres et antall 'eksperimenter' som gir posisjon av et elektron
# i rommet i henhold til sannsynlighetsfunksjonen  
# På grunn av hurtig varierende bølgefunksjon, så må man ofte gjøre mange 'forsøk'
# for å få et plott med et rimelig antall elektroner i skyen 
import numpy as np
import scipy
import scipy.special as sp
import scipy.misc as sm
import matplotlib.pyplot as plt

def psi(x,n,l):         # Define your function here
    # Den radielle bølgefunksjonen til hydrogenatomet
    #npr = n
    #nom = -4*sm.factorial(npr-l-1)*(-1)**(2*l+1)
    #denom = (a**3)*(npr**4)*(sm.factorial(npr+l)**3)
    #N = np.sqrt(nom/denom)
    # Laguerre polynomet  i Python er annerledes definert enn i Hemmer.
    # (Se kommentar i lærebok i QM av Amit Goswami)
    # Har ikke funnet ut av normeringen med denne definisjonen
    # Men man er i nærheten slik
    N = 1/sp.factorial(n)
    a = 0.529  # Bohr-radien i Ångstrøm)
    rho = 2*x/(n*a)
    y = N*(rho**l)*np.exp(-.5*rho)*sp.assoc_laguerre(rho,n-l-1,2*l+1)
    return y
# Her velges kvantetall
values = []
line = input("tast n,l og m :")
values = line.split()
n = int(values[0])
l = int(values[1])
m = int(values[2])
print("Lukk tetthetsplottet for å begynne simulering av elektronsky! (SOM KAN TA LITT TID)")
while n>0 :
    nbins = 200
    # Skala defineres utfra forventet midlere avstand
    #  som er .529*(3*n**2 - l**2-l)
    rm = 1.5*n**2
    
    x = []
    z = []
    w = []
    bw2 = .5*rm/float(nbins)
    wmax = 0.
    for i in range(nbins):
        zz = -rm + 2*i*rm/float(nbins)+bw2
        
        for j in range(nbins):
            phi = 0.  # Vi skal alltid kvadrere Ylm, så det er OK å bruke verdien for phi=0
            xx = -rm + 2*j*rm/float(nbins)+bw2
            r = np.sqrt(xx**2+zz**2)
            theta = abs(np.arccos(zz/r))
            wr = psi(r,n,l)
            ww = sp.sph_harm(int(m),int(l),phi,theta)
            weight = (r*wr*ww.real)**2
            if weight > wmax:
               wmax = weight 	
            w.append(weight)
            x.append(xx)
            z.append(zz)
    plt.hist2d(x,z,bins=nbins,weights=w)    
    plt.show()
    print("Maksimumverdi i sannsynlighetsfunksjon", wmax) # maksimumsverdi i tetthetsplottet, trengs for MC simulering
    
    # her begynner MonteCarlo-simulering
    line = input("Skriv inn antall forsøk i simulering (typisk noen hundre tusen): ")
    values = line.split()
    ntries = int(values[0]) 
    npinarray = 0
    xdata = []
    ydata = []
    zdata = []
    
    trz = np.random.random(ntries)
    trx = np.random.random(ntries)
    tryc = np.random.random(ntries)
    trtry = np.random.random(ntries)*wmax
    
    for k in range(ntries):
       # lage vilkårlig punkt i rommet
       zz = -rm + 2*trz[k]*rm
       xx = -rm + 2*trx[k]*rm
       yy = -rm + 2*tryc[k]*rm
       r = np.sqrt(xx**2 +yy**2 + zz**2)
       theta = np.arccos(zz/r)            
       wr = psi(r,n,l)
       ww = sp.sph_harm(int(m),int(l),0.,theta)
       weight = (r*wr*ww.real)**2  
       # sjekke sannsynligheten av dette punktet, akseptere med frekvens i henhold til sannsynlighet         
       if weight > trtry[k]:
          npinarray = npinarray + 1
          zdata.append(zz) 
          xdata.append(xx) 
          ydata.append(yy) 
    ax = plt.axes(projection='3d')
    ax.scatter3D(xdata,ydata,zdata,s=.01) 
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    print(" Antall punkter i skyen ",npinarray) 
    plt.show()   # Show the plot (Without this you don't see anything on the screen)
    line = input("tast n,l og m :")
    values = line.split()
    n = int(values[0])
    l = int(values[1])
    m = int(values[2])   
exit()
