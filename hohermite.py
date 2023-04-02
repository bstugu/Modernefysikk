# Plotte egenfunksjoner for Harmonisk oscillatorer
# I det nederste plottet sammenlignes klassisk 
# og kvantemekanisk sannsynlighet for å finne partikkel i x
# vi setter hbar=1

import numpy as np
import scipy.special as sp
#import scipy.misc as sm
import matplotlib.pyplot as plt

def f(x,m,omega,n,p):        
    # Egenfunksjon som løser Schrødingerligningen for harmonisk
    # oscillator med kvantetall n
    N = 1./np.sqrt((2.**n)*sp.factorial(n))
    xx = x*np.sqrt(m*omega)
    y = (sp.eval_hermite(n,xx)*np.exp(-xx*xx/2.)*N)**p
    return y
def cl(x,m,omega,n):
    # Klassisk sannsynlighet for partikkelens posisjon
    # (ikke normert)
    E = (n+.5)*omega
    a = np.sqrt(2*E/m)/omega
    y = 1/(omega*np.sqrt(a**2-x**2))  # Python lar oss fortsette selv om vi iblant har x>a
    return y
# 
#  Hovedprogram
#
n = 20   # Kvantetall som gir energinivå (1/2+n)*hbar*omega
p = 1   #Bølgefunksjon: p=1, eller sannsynlighetsfordeling (psi**2): p=2
N = 1./np.sqrt((2.**n)*sp.factorial(n))  # Normeringskonstant
# Tilfeldig valgt masse og karakteristisk frekvens
m = 1.
omega = 1.
A = np.sqrt(2*(n+.5)/(m*omega))
print(" N ",N, " Amplitude ", A)
fig, (ax1,ax2) = plt.subplots(nrows=2,ncols=1)
xl = -np.sqrt((2*n+1)/omega)
xh = np.sqrt((2*n+1)/omega)
x1 = np.arange(xl*2,xh*2,0.001)
x = np.arange(xl*2,xh*2,0.001)
y = f(x,m,omega,n,p)    

ax1.plot(x,y) 

y1 = cl(x1,m,omega,n)
y2 = f(x,m,omega,n,2)
ax2.plot(x1,y1,'r--')
ax2.plot(x1,y2)
# Vi må finne et fornuftig maksimum på y-akseskala
maxv = np.max(y2)
maxy =1.1*maxv
ax2.axis([xl*2,xh*2,0.,maxy])

plt.show()   # Vise plottet
