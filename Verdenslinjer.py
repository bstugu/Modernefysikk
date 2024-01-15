import scipy
import numpy as np
import matplotlib.pyplot as plt
import math
# Her kan man plotte en verdenslinje for en relativistisk tur 
# og finne tider for mottak av meldinger for turen
# Vi setter lyshastigheten lik 1 i dette programmet
# Hastighet til objekt som reiser ut og hjem
beta = 0.5
gamma = 1/np.sqrt(1.-beta*beta)
print(" beta ",beta," gamma ",gamma)
# Returnering ved første glimt etter dette tidspunktet
rt = 5.
# Signalfrekvens i rakettens hvilesystem til hjemsendte signaler
f0 = 1.5



# Dopplerforskjøvede frekvenser
f1 = np.sqrt((1+beta)/(1-beta))*f0
f2 = np.sqrt((1-beta)/(1+beta))*f0
dT0 = 1/f0
dT1 = 1/f1
dT2 = 1/f2

# Tid mellom hvert signal sendt fra romskip, som observert hjemme
dTlab = dT0*gamma
Tlab = 0.

# Plottegreier
# Ytterkanter av figuren
XHI = 6.
XLO = -1.
YLO = 0.
YHI = 12.
plt.axis([XLO,XHI,YLO,YHI])
plt.xlabel('x')
plt.ylabel('ct')
#linja ct = x
x1=[0.,10.]
y1=[0.,10.]
plt.plot(x1,y1,'ro--',linewidth=0.5,markersize=.5)
slope = beta
# xtick er tilbakelagt reiselengde i som sett i lab
# ytick er ct i lab
xtick = 0.
ytick = 0.
# tilbakelagt reiselengde mellom hvert blink
dxtick = dTlab*beta
dytick = dTlab

dprimex = []
dprimey = []
# Dopplerformeltid (doptx velges for plassering i plottet
dopt = []
dopx = []
while ytick < rt :
   xtick = xtick + dxtick
   ytick = ytick + dytick
   dprimex.append(xtick)
   dprimey.append(ytick)
   Tlab = Tlab + dT2
   dopt.append(Tlab)
   dopx.append(-.3)
   x1 = [xtick,0.]
   y1 = [ytick,ytick+xtick]
   plt.plot(x1,y1,'ro--',linewidth=0.5,markersize=.5)
snux = slope*ytick
x1=[0.,snux]
snuy = ytick
y1=[0.,snuy]
plt.plot(x1,y1,'k-',linewidth=1.,markersize=1)   
while ytick < 2*rt :
   Tlab = Tlab + dT1
   xtick = xtick - dxtick
   ytick = ytick + dytick
   dprimex.append(xtick)
   dprimey.append(ytick)
   dopt.append(Tlab)
   dopx.append(-.3)
   x1 = [xtick,0.]
   y1 = [ytick,ytick+xtick]
   plt.plot(x1,y1,'ro--',linewidth=0.5,markersize=.5)  
plt.plot(dprimex,dprimey,'ko',markersize=3.)
plt.plot(dopx,dopt,'kx',markersize=3.)
x1=[snux,0.]
y1=[snuy,snuy*2]
plt.plot(x1,y1,'k-',linewidth=1.,markersize=1)
#plt.arrow(1.,1.,4.,0.,head_width=.2,head_length=.5,fc='k',ec='k')
#plt.text(4.5,.5,'x',fontsize=sz)


  
plt.title("Verdenslinjer og Doppler-tider")
plt.grid() 
#plt.axis('off')
plt.show()
