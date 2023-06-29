# Todimensjonalt bilde av hydrogenatomets bølgefunksjon eller sannsynlighet i z-x planet
# Man velger kvantetall n,l,m og produserer plott med funksjonsverdier i passende fargeskala (automatisk valgt av plottefunksjonen)
import numpy as np
import scipy
import scipy.special as sp
import scipy.misc as sm
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
def psi(x,n,l):         # Define your function here
    # Den radielle bølgefunksjonen til hydrogenatomet
    #npr = n
    #nom = -4*sm.factorial(npr-l-1)*(-1)**(2*l+1)
    #denom = (a**3)*(npr**4)*(sm.factorial(npr+l)**3)
    #N = np.sqrt(nom/denom)
    # Laguerre polynomet  i Python er annerledes definert enn i Hemmer.
    # (Se kommentar i lærebok i QM av Amit Goswami)
    # Har ikke funnet ut av normeringen med denne definisjonen
    # men setter inn noe som ikke er helt på trynet
    N = 1/sp.factorial(n)
    a = 0.529  # Bohr-radien i Ångstrøm)
    rho = 2*x/(n*a)
    y = N*(rho**l)*np.exp(-.5*rho)*sp.assoc_laguerre(rho,n-l-1,2*l+1)
    return y
    
# Her velges kvantetall
# I denne versjonen kan man også plotte psi med korrekt fortegn
#Programmet avsluttes ved å velge valg = 0
# Du må lukke plottvinduet for å velge nye kvantetall

values = []
print("Valg = 1: psi**2  Valg 2: psi Valg 3: abs(psi)  Valg >= 4:  (r*psi)**2")
line = input("Tast valg (0=avslutt): ")
valg = int(line)
if valg > 0:
   line = input("tast n,l og m (MELLOMROM mellom tallene) ,n=0 betyr avslutt:")
   values = line.split()
   n = int(values[0])
   l = int(values[1])
   m = int(values[2])

while valg>0 :
    costh = -1.
    dcosth = .01
    bins = int(2/dcosth)
    norm = 0.
    for i in range(bins):
       th = np.arccos(costh)
       wttst = sp.sph_harm(int(m),int(l),0.,th)
       wrtst = (wttst.real)**2
       #print("wttst ",wttst," wrtst ",wrtst)
       norm = norm + wrtst*dcosth
       costh = costh+dcosth
    print("norm ",6.28*norm)
    lowest = 1.0e20
    highest = -1.0e20
    nbins = 200
    # Skala defineres utfra forventet midlere avstand
    #  som er .529*(3*n**2 - l**2-l)
    rm =  n**2
    #rm = 2.*n
    x = []
    z = []
    w = []
    bw2 = .5*rm/float(nbins)
    for i in range(nbins):
        zz = -rm + 2*i*rm/float(nbins)+bw2
        for j in range(nbins):
            phi = 0.  # Vi bruker verdien for phi=0
            xx = -rm + 2*j*rm/float(nbins)+bw2
            r = np.sqrt(xx**2+zz**2)
            theta = abs(np.arccos(zz/r))
            wr = psi(r,n,l)
            ww = sp.sph_harm(int(m),int(l),phi,theta)
        
            if valg == 1:
                weight = (wr*ww.real)**2
            elif valg == 2:
                weight = wr*ww.real
                if weight > highest:
                    highest = weight
                if weight < lowest:
                    lowest = weight
            elif valg == 3:
                weight = np.abs(wr*ww.real)
                if weight > highest:
                    highest = weight
                if weight < lowest:
                    lowest = weight           
            else:
                weight = (r*wr*ww.real)**2
                if weight > highest:
                    highest = weight
                if weight < lowest:
                    lowest = weight               	
            w.append(weight)
            x.append(xx)
            z.append(zz)
    fig,ax = plt.subplots()
    h = ax.hist2d(x,z,bins=nbins,weights=w)
    fig.colorbar(h[3],ax=ax)
    if valg == 1:
        title = "psi**2,  n = "+str(n)+" l = "+str(l)+" m = "+str(m) 
    elif valg == 2:
        title = "psi,  n = "+str(n)+" l = "+str(l)+" m = "+str(m)
    elif valg == 3:
        title = "|psi|,  n = "+str(n)+" l = "+str(l)+" m = "+str(m)
    else:
        title = "P(x,z)  n = "+str(n)+" l = "+str(l)+" m = "+str(m)
    plt.title(title)            
    plt.xlabel(" x/$10^{-10}$ m ")
    plt.ylabel(" z/$10^{-10}$ m ")
    #plt.hist2d(x,z,bins=nbins,weights=w,norm=LogNorm())
    print(" Lowest ", lowest)
    print(" Highest",highest)
    plt.show()   # Show the plot (Without this you don't see anything on the screen)
    line = input("Tast valg (0 betyr avslutt): ")
    valg = int(line)
    if valg > 0:
       line = input("tast n,l og m :")
       values = line.split()
       n = int(values[0])
       l = int(values[1])
       m = int(values[2])   
exit()
