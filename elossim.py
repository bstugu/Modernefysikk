# Simulere og plotte energitap for ioniserende
# partikler som går gjennom et materiale.
# Tre valg:
#
# 1: Simulere differensielt energitap (dE/dx) og plotte fordeling av tapene i henhold til Landaufluktuasjoner
# 2: Plotte Bethe-Bloch formelen som funksjon av prosjektilenergi
# 3: Plotte energitap når et prosjektil trenger inn i et tykt materiale, helt til det stopper 
#    Dette vil gi et omtrentlig bilde av hvordan energitapet fordeler seg innover i materialet.
#    Det er lagt på Landau-fluktuasjoner i hvert skritt
#    Resultatet vil bl.a. avhenge av valgt skrittlengde. 
#    Derfor må skrittlengde 'kalibreres' mot kjente resultater.
#    Bethe-Bloch resultatet vises også 
#    Fra Google:
#    "The range of a 160 MeV proton is 17.6 cm in water. In air it is 16700 cm (167 m)".
#    Passende skrittlengde er minst 1/1000 av forventet rekkevidde 
#  en alpha partikkel på 5 MeV skal ha en rekkevidde på 3,6 cm i luft
import numpy as np
import scipy.special as sp
#import scipy.misc as sm
import matplotlib.pyplot as plt

def landau(x,N):         # Definer funksjonen
    # Tilnærming av landaufordeling som nevnt i wikipedia den 8/6-2022
    # N er normeringsfaktor, lik 1/sqrt(2*pi)= .399 om integralet skal være 1
    # N=1.6 gir maxverdi omtrent lik 1
    y = N*np.exp(-0.5*(x+np.exp(-x)))
    return y

def BB(beta,q,Z,A,rho,I,dx):
    # Bethe-Bloch formelen for gjennomsnittlig differensielt energitap 
    # som går gjennom et stykke dx i et materiale 
    # beta: Partikkelens hastighet  ( = v/c)
    #    q: Partikkelens ladning
    #    Z: Mediets kjerneladning
    #    A: Mediets nukleontall
    #  rho: Mediets tetthet
    #    I: Ioniseringspotensial    
    gamma = np.sqrt(1./(1.-beta**2))
    eta = beta*gamma
    Wmax = 1022000*eta
    # konstanten K = .3070 MeVcm**2mol**-1 Nedenfor står K/2
    y = dx*.1535*rho*(Z/A)*q**2*(np.log(1022000*gamma**2*beta**2*Wmax/(I**2))-2*beta**2)/beta**2
    return y

def landausim(N,dx,q,beta,rho,Z,A,I):
    # Generere energitap i materiale av tykkelse dx, med Landaufluktuasjon
    # N events genereres
    #Kan brukes f.eks. for simulering av energitap i tynne sensorer (f.eks. av silisium)
    Deltam = BB(beta,q,Z,A,rho,I,dx)
    xsi = 0.1535*q**2*Z*rho*dx/(A*beta**2)
    gamma = 1./np.sqrt(1.-beta**2)
    eta = beta*gamma
    Wmax =  1.022000*eta
    k = xsi/Wmax
    logeprime = np.log(xsi*1000000/I)
    #logeprime = np.log(I**2/(gamma**2*beta**2*1022000))+beta**2
    Deltamp = xsi*(np.log(xsi)-logeprime +.198)*dx
    #print("Deltam ",Deltam," Deltamp ",Deltamp," xsi ",xsi," k ",k," logeprime ",logeprime)
    Deltash = []
    for itr in range(N):

       xrndm = -2.8 +25.*np.random.random()
       yrndm = landau(xrndm,1.6)
       yrndm2 = np.random.random()
       while yrndm2 > yrndm:
          xrndm = -2.8 +25.*np.random.random()
          yrndm = landau(xrndm,1.6)
          yrndm2 = np.random.random()     
 
       Deltas = Deltam+ xsi*(xrndm + beta**2 + np.log(k)  +1 - 0.577)
       Deltash.append(Deltas)
    plt.title("Fordeling av energitap")
    plt.xlabel("energitap (MeV)")
    plt.ylabel("Antall")
    plt.hist(Deltash,bins=100)
    plt.show()
    return

def braggrule(ai,fZ,fA,N):
    # Finne gjennomsnitt av elementene i materialet
    # og gjennomsnittlig ioniseringspotensial
    # ai er relativt antall av hvert element
    Am = 0.
    Zm = 0.
    sai = 0.
    Ii = []
    for i in range(N):
       sai = sai + ai[i]
    Zeff = 0
    Aeff = 0
    Iw = 0
    for i in range(N):
       Zeff = Zeff + ai[i]*fZ[i]
       Aeff = Aeff + ai[i]*fA[i]
    Zeff = Zeff/sai
    Aeff = Aeff/sai
    for i in range(N): 
       # Ioniseringspotensial fra Leo
       if fZ[i] < 13:
          Ii.append((12 + 7./fZ[i])*fZ[i])
       else:
          Ii.append((9.76 + 58.8*fZ[i]**(-1.19))*fZ[i])
       Iw =  Iw + np.log(Ii[i])*ai[i]
       print(" I ",Ii[i])
    Iw = Iw/sai
    I = np.exp(Iw)
    print( "Ieff ",I," Zeff ",Zeff," Aeff ",Aeff)
    return Zeff,Aeff,I 


def rangesim(mp,Ekstart,q,Nsim,dt,ibbl):
    # simulere energitap til Nsim partikler inntil de stopper
    # mp = partikkelmasse, Ekstart = kinetisk energi (E-masse),q=ladning
    # dt = skrittlengde
    # hvis ibbl =1 så legges BetheBloch prediksjon til slutt
    # Først vil alle energitapene vises, så vil forløp av gjennomsnittlig energitap vises. 
    tt = []
    ddE = []
    Rng = []
    Maxiter = 0
    for isim in range(Nsim+ibbl):
        Ek = Ekstart
        E = Ekstart+mp
        P = np.sqrt(E**2-mp**2)
        beta = P/E
        tn = 0.
        t = []
        dE = []
        deBs = []
        Iter = 0
        de = 0.
        while Ek > de:
            Iter = Iter + 1 
            dx = dt
            xrndm = -2. +25.*np.random.random()
            yrndm = landau(xrndm,1.6)
            yrndm2 = np.random.random()
            while yrndm2 > yrndm:
                xrndm = -2. +25.*np.random.random()
                yrndm = landau(xrndm,1.6)
                yrndm2 = np.random.random()      
            Deltam = BB(beta,q,Z,A,rho,I,dx)
            xsi = 0.1535*q**2*Z*rho*dx/(A*beta**2)
            gamma = 1./np.sqrt(1.-beta**2)
            eta = beta*gamma
            Wmax =  1.022000*eta
            k = xsi/Wmax
            Deltas = Deltam+ xsi*(xrndm + beta**2 + np.log(k)  +1 - 0.577)      
            de = Deltas
            #print("Deltas ",Deltas," Deltam ",Deltam," xsi ",xsi)
            if isim == Nsim:
                de = Deltam
            if Ek > de:
                t.append(tn+dt)
                dE.append(de)
                deBs.append(Deltam)
                P = np.sqrt(E**2-mp**2)
                beta = P/E
                Ek = Ek - de
                E = Ek + mp
                tn = tn+dt
            else:
                t.append(tn+dt)
                dE.append(Ek)
                deBs.append(Ek)
                tn = tn + dt
                Ek = 0.
        #print("Ek ",Ek)
        Rng.append(tn)
        t.append(tn+dt)
        dE.append(0.) 
        deBs.append(0.)
        if Iter > Maxiter:
            Maxiter = Iter
        print("ISIM ", isim," Antall skritt ",Iter, "Rekkevidde ",tn+2*dt)
        tt.append(t)
        ddE.append(dE)
        plt.plot(tt[isim],ddE[isim])
 
    for isim in range(Nsim+ibbl):
        last = len(tt[isim])
        ttn = Rng[isim]
        if last <= Maxiter:
            for ind in range(last,Maxiter+1):
                ttn = ttn + dt
                tt[isim].append(ttn)
                ddE[isim].append(0.)
   
    plt.title("Energitap som funksjon av dybde")
    plt.ylabel("energitap pr. skritt (MeV)")
    plt.xlabel("Dybde (cm)")  
    plt.show()
    avg = []
    for ind in range(Maxiter+1):
        avg.append(0)
    for simind in range(Nsim):
        for ind in range(Maxiter+1):
            avg[ind] = avg[ind] + ddE[simind][ind]
    for ind in range(Maxiter+1):
        avg[ind] = avg[ind]/Nsim      
    plt.plot(tt[0],avg)   
    if ibbl == 1:
        plt.plot(tt[0],ddE[Nsim])    
    plt.title("Gjennomsnittlig energitap som funksjon av dybde")
    plt.ylabel("energitap pr. skritt (MeV)")
    plt.xlabel("Dybde (cm)") 
    plt.show()
    return
    
def safeinput(txt,nvar):
    # Lese tall uten kræsj
    # Antar at korrekt input er nvar tall separert med blanke
    OK = False
    line = input(txt)
    values = []
    floatvalues =[]
    while OK == False:
        values = line.split()
        if len(values) == int(nvar):
            try:              
                for i in range(0,nvar):
                    floatvalues.append(float(values[i]))
                OK = True
            except:
                print("Tast ",nvar," TALL, NB desimalskilletegn er punktum")
                line = input(txt)
        else:
            print("Tast ",nvar," tall separert med mellomrom: ")
            line = input(txt)
    return floatvalues
    
# HOVEDPROGRAM
# Valg: = 1: Simulere differensielt energitap (dE/dx) og plotte fordeling av tapene for NN events 
#         2: Plotte Bethe-Bloch formelen som funksjon av prosjektilenergi
#         3: Plotte energitap når et prosjektil trenger inn i et tykt materiale, helt til det stopper



# Lage Landaufunksjon på intervall xl til xh i array y
xl = -4.
xh = 14.
x = np.arange(xl,xh,0.001)
y = landau(x,1.6)
# hjelpevariabel for innlesning
par = []
# Lese inn materialet: (en blanding defineres isteden rett i programmet) 
# (silisium har Z = 14, A = 28 og rho = 2.329 g/cm**3)
par = safeinput("Tast inn materialets Z, A og tetthet (0 0 0 -> egendefinert blanding) ",3)
Z = float(par[0])
A = float(par[1])
rho = float(par[2])
if Z < .5:
    #Egendefinert materialblanding
    # Vann
    
    Nel = 2
    fZ = [1,8]
    fA = [1,16]
    ai = [2,1]
    Nel = 2
    rho = 1.
    """
    #luft
    Nel = 2
    fZ = [7,8]
    fA = [14,16]
    ai = [78,22]
    rho = 1.2/1000.
    
    #ScintillatorPlast (fant materialet EJ-296 på nettet)
    Nel = 2
    fZ = [1,6]
    fA = [1,12]
    ai = [5.17,4.69]
    rho = 1.02
    """  
    Z,A,I = braggrule(ai,fZ,fA,Nel)  

print("Materialet som traverseres har Z= ",Z," A= ",A," Tetthet= ",rho, " g/cm**3")
# Ioniseringspotensial fra Leo
if Z < 13:
    I = 12 + 7./Z
else:
    I = 9.76 + 58.8*Z**(-1.19)
#print("I ",I)
par = safeinput("Tast valg (1,2 eller 3)",1)
Valg = int(par[0])
if Valg == 1:
    print(" Simulere energitap i tynt materiale")
    par = safeinput("Tast prosjektilmasse, impuls (MeV/c**2 og MeV/c) og ladning ",3)
    mp = float(par[0])
    p = float(par[1])
    q = float(par[2])
    par = safeinput("Tast materialtykkelse og antall events",2)
    Dx = float(par[0])
    Nsim = int(par[1])
    beta = p/np.sqrt(mp**2+p**2)
    landausim(Nsim,Dx,q,beta,rho,Z,A,I)

if Valg == 2:
    print(" Plotte (dE/dx)*dx i impulsintervall i henhold til Bethe-Bloch formel")
    par = safeinput("Tast prosjektilmasse (MeV/c**2, nedre og øvre impuls (MeV/c,MeV/c), og ladning  ",4)
    mp = float(par[0])
    pl = float(par[1])
    ph = float(par[2])
    q = float(par[3])
    par = safeinput("Tast materialtykkelse ",1)
    Dx = float(par[0])   
    p = np.arange(pl,ph,.005)
    Ex = np.sqrt(p**2+mp**2)
    Beta = p/Ex
    yBB = BB(Beta,q,Z,A,rho,I,Dx)
    plt.plot(p,yBB)
    plt.title("Energitap som funksjon av impuls")
    plt.xlabel("Impuls (MeV)")
    plt.ylabel("dE (MeV)")
    plt.show()   # Vise plottet

if Valg == 3:
    print("Simulere energitapet for et antall partikler som går gjennom materiale inntil de stopper")
    par = safeinput("Tast prosjektilmasse, Kinetisk energi (MeV/c**2, MeV) og ladning ",3)
    mp = float(par[0])
    Ekstart = float(par[1])
    q = float(par[2])
    
    par = safeinput("Tast skrittlengde i integrasjon, og antall events ",2)
    dt = float(par[0])
    Nsim = int(par[1])
    par = safeinput(" Vil du legge på deterministisk Bethe-Bloch prediksjon? (1/0)",1)
    ibbl = int(par[0])
    #  Skrittlengde i cm: Må avhenge av Ek
    # Rekkevidde til protoner i vann er 0,12 cm for 10 MeV kE, 
    # vil ha minst 1000 skritt her (Må iterere)
    rangesim(mp,Ekstart,q,Nsim,dt,ibbl)

