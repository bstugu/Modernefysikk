#
# Diffraksjon/interferens plott
# Vi får et plott av intensitetsfordeling når parallelt lys går gjennom en
# eller to spalter.
# Beregner også standardavvik i gitt intervall, når intensitetsfunksjonen 
# tolkes som en statistisk fordeling av fotoner.
# Dette sammenlignes med usikkeren i transvers komponent av fotonimpuls som forutsagt 
# av Heisenberg. Beregningene har bare mening for énspalte tilfellet  
 
import numpy as np
import scipy
from scipy.special import *
from scipy.stats import *
import matplotlib.pyplot as plt

def f(x,a,b,bl):
# f gir intensitet for vinkel x  
# a: avstand mellom spalter
# b: Bredden til hver spalte
# bl: Bølgelengde
   arg1 = 3.14159*b*np.sin(x)/bl
   diftrm = np.sin(arg1)/arg1
   arg2 = 3.14159*a*np.sin(x)/bl
   inttrm = np.cos(arg2)
   return (diftrm*inttrm)**2
  
# hovedprogram
#
b = 40. # Bredden til hver spalte
a = 0. # Avstand mellom spaltene, når a=0, så har vi bare en spalte med bredde b
bl = 10. # Bølgelengde
# Selve diffraksjonsmønsteret avhenger bare av forholdet mellom bl og b
hbar = 1.
c = 1.
# Lage x-verdier (vinkler)
xl = -.25
xh = .25
stp = 0.0001
th = np.arange(xl,xh,stp)
npfu = int((xh-xl)/stp) + 1
plt.axis([xl,xh,-.0,1.3]) #nedre og øvre grenser i plott


y = f(th,a,b,bl)
#plt.plot(th,y)
plt.plot(th,y)
plt.title("Intensitetsfordeling")
plt.xlabel('theta')
plt.ylabel('intensitet')
 #   Integrering er en lek!
intv = scipy.integrate.trapz(y,th)
print(" Verdi av funksjonens integral fra ",xl," til ",xh,f" er: {intv:5.3f}" )
print(" (trapes-regel med ", npfu,"  punkter)")
middelv = scipy.integrate.trapz(y*th,th)/intv
intv2 = scipy.integrate.trapz(y*th*th,th)/intv
sd = np.sqrt(intv2-middelv**2)
print(f" Verdi av middelverdi og standardavvik {middelv:.3f}   {sd:.3f}")
E = hbar*2*3.14159*c/bl
p = E/c
pt = p*np.sin(sd)
Dx = b/np.sqrt(12)
DpDx=pt*Dx
hbarCn = 197. #(eV nm)
hbarC = hbarCn/(hbar*c)
# Vi gjør om til elektronvolt og meter
print("=================================")
print(f" Bølgelengden {bl:.3f} nm svarer til en foton-bevegelsesmengde: {hbarC*E/c:.3f} eV/c")
print(f" Usikkerhet i bevegelsesmengde (fra standardavvik i intensitetsplott) : {hbarC*pt:.3f} eV/c")
print(f" Bredde/sqrt(12) av spalt, Dx: {Dx:.3f}  nm ")
print(f" DpDx: {hbarC*pt*Dx:.3f} eVnm ")
print(f" verdien av hbar/2: {hbarCn/2.:.3f} eV nm")
print("=================================")

plt.show()
