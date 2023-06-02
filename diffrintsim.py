#   Simulere enkeltfotoner som passerer gjennom dobbeltspalt og treffer 
#   en skjerm i henhold til forventet intensitetsfordeling til en elektromagnetisk bølge 
#   Fotonene genereres ved Monte-Carlo simulering
import numpy as np
import scipy
from scipy.special import *
from scipy.stats import *
import matplotlib.pyplot as plt

def f(x,a,b,bl):
   arg1 = 3.14159*b*np.sin(x)/bl
   diftrm = np.sin(arg1)/arg1
   arg2 = 3.14159*a*np.sin(x)/bl
   inttrm = np.cos(arg2)
   return (diftrm*inttrm)**2
  
# hovedprogram
#
# Lage x-verdier (vinkler)
xl = -1.
xh = 1.
stp = 0.0001
#  Antall forsøk (dette er IKKE antall fotoner)  
ntries = 100000
# thetas er en array med ntries vinkler, uniformt fordelt 
thetas = xl + (xh-xl)*np.random.random(ntries)
# yl er en array med ntries sannsynligheter
yl = np.random.random(ntries)
ys = []
xs = []
b = .25  # Bredden til hver spalte
a = 1.  # Avstand mellom spaltene
bl = .25 # Bølgelengde
hbar = 1.
c = 1.
l = 1. # Avstand fra spalt til skjerm
nfotoner = 0
for i in range(ntries):
   th = thetas[i]
   
   y = f(th,a,b,bl)
   yll = yl[i]
   if y > yll:
      # nå har vi generert et foton som har passert spaltene med vinkel th
      xs.append(th)
      # Vi generer et vilkårlig punkt langs spaltenes lengderetning
      ypos = np.random.random(1)
      ys.append(ypos)
      nfotoner = nfotoner + 1
print(" Antall fotoner som treffer skjerm: ",nfotoner)
plt.plot(xs,ys,'o',markersize=1.)
title0 = " fotoner"
tit = str(nfotoner)+title0
plt.title(tit)
plt.show()
