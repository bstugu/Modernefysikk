#   Simulere magnetisering for system (paramagnet) som har spontane spinnflipp
#   Sammeligne med den deterministiske formelen ('teorien').
#   Utvikling av entropi vises også
#   Første del starter med fullstendig magnetisert 'stav' med magnetfelt avslått

import numpy as np
import scipy
from scipy.special import *
from scipy.stats import *
import matplotlib.pyplot as plt


#  Antall spinnpartikler er npos+nneg  
npos = 1000.
nneg = 0.
initnpos = npos-nneg
initnneg = nneg-npos
sm = int(npos+nneg)
ppos =npos/(npos+nneg)
# Sannsynlighet for splinnflipp pr. sekund
p = 0.05
# Antall skritt i simulering uten magnetfelt
duration = 100
# Antall skritt etter magnetfelt er slått på
durationB = 100
# Energiforskjell mellom de to spinntilstandene når magnetfelt er på 
dEkT = .2
Magtest = np.random.random(sm) 
Mag = []
for i in range(sm):
   if Magtest[i] < ppos:
      Mag.append(1)
   else:
      Mag.append(-1)

ys = []
xs = []
entr = []
TMag = npos-nneg
initTMag = TMag
thMag = []
xsth = []  # vi stopper teorien før vi slår på B-felt
nthpos = npos
for i in range(duration):
   xs.append(i)
   ys.append(TMag/(npos+nneg))
   #   Teoretisk forløp
   if initTMag > 0:
      nthpos= 0.5*initnpos*(1+np.exp(-2*p*i)) 
      nthneg = sm - nthpos
   else:
      nthneg= np.abs(0.5*initnneg*(1+np.exp(-2*p*i)))
      nthpos = sm - nthneg
   thMag.append((nthpos-nthneg)/sm)
   xsth.append(i)


   fliptest = np.random.random(sm)
   TMag = 0.
   nNpos = 0.

   for j in range(sm):
      if fliptest[j] < p:
         if Mag[j] == 1:
            Mag[j] = -1
         else:
            Mag[j] = 1
      if Mag[j] > 0.:
         nNpos = nNpos+1.      
      TMag = TMag + Mag[j]
   if nNpos > 0.:
      # entropi ved Stirlings formel 
      entr.append(sm*(np.log(sm)-1) - (nNpos*(np.log(nNpos)-1)+(sm-nNpos)*(np.log(sm-nNpos)-1)))
   else:
      entr.append(0.)
  

   

nBpos = np.exp(-dEkT)/(np.exp(dEkT)+np.exp(-dEkT))
nBneg = np.exp(dEkT)/(np.exp(dEkT)+np.exp(-dEkT))
  
popp = p*(1-(nBneg-nBpos)/(nBpos+nBneg))
pned =  p*(1+(nBneg-nBpos)/(nBpos+nBneg))
dur2 = duration + durationB 
for i in range(duration,dur2):   
   xs.append(i)   
   fliptest = np.random.random(sm)
   TMag = 0.
   nNpos = 0.
   nNneg = 0.
   for j in range(sm):
      if Mag[j] == 1:
         if fliptest[j] < pned:
            Mag[j] = -1
      else:
         if fliptest[j] < popp:
            Mag[j] = 1
      if Mag[j] > 0.:
         nNpos = nNpos + 1.
      TMag = TMag + Mag[j]  
   ys.append(TMag/(npos+nneg))
   entr.append(sm*(np.log(sm)-1) - (nNpos*(np.log(nNpos)-1)+(sm-nNpos)*(np.log(sm-nNpos)-1)))
fig, (ax1,ax2) = plt.subplots(nrows=2,ncols=1,layout="constrained")
ax1.set_title('Magnetisering')
ax1.plot(xs,ys)
ax1.plot(xsth,thMag)
ax2.set_title('Entropi/k')
ax2.plot(xs,entr)
plt.show()
