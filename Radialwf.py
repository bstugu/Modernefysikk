# Plotte Radielle bølgefunksjoner for Hydrogenatomet
# R(r) og P(r)=r**2*R/r)**2 plottes
# Bølgefunksjoner for opp til fire verdier av l kan plottes i samme bilde.
# Ønskede verdier av n og l settes rett inn i koden
# Det beregnes også forventningsverdi av r (med 'ad hoc' øvre grense for r)

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
    # Har ikke funnet ut av normeringen med denne definisjonen, setter noe som
    # ikke er helt på trynet
    N = 1/sp.factorial(n)
    a = 0.529  # Bohr-radien i Ångstrøm)
    rho = 2*x/(n*a)
    y = N*(rho**l)*np.exp(-.5*rho)*sp.assoc_laguerre(rho,n-l-1,2*l+1)
    return y

#Sett verdier inn her
# Hovedkvantetall
n = 10
# Ønskede verdier av l
nl = 4
ll = [6,7,8,9]
clor = ['r','b','g','k']
lstyle = ['solid','dashed','dotted','dashdot']

fig, (ax1,ax2) = plt.subplots(nrows=2,ncols=1)

for i in range(0,nl):
   l = ll[i]
   xl = 0.
   # Skala kan defineres utfra forventet midlere avstand
   #  som er .5*(3*n**2 - l**2-l)
   #xh = 3.*n**2
   xh = 150.  # ..eller manuelt
   x = np.arange(xl,xh,0.01)
   y = psi(x,n,l)    # Bølgefunksjonen

   y1 = (x**2)*(y**2)
   # normering
   intv = np.trapz(y1,x)
   y = y/np.sqrt(intv)
   y1 = y1/intv
   
   ax1.plot(x,y,color=clor[i],linestyle=lstyle[i]) # 'make' the plot (array of points connected with lines)
   plt.ylabel("                               P(r)                                            $\psi (r)$")
   plt.xlabel(" $r/10^{-10}$m ")
   ax2.plot(x,y1,color=clor[i],linestyle=lstyle[i])
   

   y2 = x*y1
   intv = np.trapz(y2,x)
   intv = intv/.529
   print ("n= ",n," l= ",l," Forventningsverdi av avstand, <x>= ","%9.3f"%intv," Bohrradier")
plt.show()   # Show the plot (Without this you don't see anything on the screen)
