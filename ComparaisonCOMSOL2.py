#Auteur : Mathieu Goureau/ Arona Sottas
#mathieu.goureau@polymtl.ca

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

##Définition des constantes

S =  1e-8 #Terme source constant
k = 4e-9
D_eff = 1e-10 #Coefficient de diffusion constant
Ce = 10 #Concentration extérieure
D = 1 #Diamètre du pilier
R = D/2 #Rayon du pilier

Vn = np.array([501,1501,3001]) #Choix du nombre de points dans le maillage

dt = 1e6 #Base de temps en seconde
Vt = np.arange(0,1e10,dt) #Vecteur des temps t

##Lecture fichier COMSOL :

solCOMSOL = []

COMSOL_501 = np.loadtxt("Valeurs1D_S_avec_kCquadratique_header_501_1e6.txt",comments="%")
COMSOL_1501 = np.loadtxt("Valeurs1D_S_avec_kCquadratique_header_1501_1e6.txt",comments="%")
COMSOL_3001 = np.loadtxt("Valeurs1D_S_avec_kCquadratique_header_3001_1e6.txt",comments="%")

solCOMSOL.append(COMSOL_501[:,1])
solCOMSOL.append(COMSOL_1501[:,1])
solCOMSOL.append(COMSOL_3001[:,1])

##Construction des matrices

def A2(alpha,V,dr):
    '''Coefficients D de la matrice avec schémas à l'ordre 2'''
    T=np.ones_like(V)
    for i in range(len(V)-1):
        T[i]=V[i+1]
    return -alpha*((1/dr**2)-1/(2*T*dr))

def B2(alpha,dr):
    '''Coefficients B de la matrice avec schémas à l'ordre 2'''
    return 1+(k*dt)+alpha*((2/dr**2))

def D2(alpha,V,dr):
    '''Coefficients D de la matrice avec schémas à l'ordre 2'''
    T=np.ones_like(V)
    T[2:]=V[1:-1]
    return -alpha*((1/dr**2)+1/(2*T*dr))
    
def MatriceGear(V,Ntot):
    ''' Création de la matrice avec schéma de Gear'''
    dr = R/(Ntot-1)
    alpha = dt*D_eff
    A = -alpha/dr**2
    M = A2(alpha,V,dr)*np.eye(Ntot,k=-1) + B2(alpha,dr)*np.eye(Ntot) + D2(alpha,V,dr)*np.eye(Ntot,k=1)
    M[0],M[-1]=0,0
    M[0,0],M[0,1],M[0,2]=-3,4,-1 #Schéma de Gear
    M[-1,-1] = 1
    return M

def Euler_implicite_solve(Vr,Vt,M,Y0,Matrice,Ntot):
    """ Résolution du système linéaire à chaque t selon la méthode Euler implicite"""
    solu=[] 
    Y=Y0 #Iniatialisation de la solution
    solu.append(Y0.tolist())
    M=Matrice(Vr,Ntot)
    for t in Vt: #Résolution à chaque t du système AX=B
        C=VecC(Y,t,Vr,Ntot) #Vecteur B
        Y=np.linalg.solve(M,C) #Résolution du système
        solu.append(Y.tolist())
    return np.array(solu)

def VecC(Y,t,Vr,Ntot):
    '''Renvoie le vecteur colonne C'''
    #source = (D_eff*t*(-4*Vr*np.sin(Vr)+(R**2-Vr**2)*np.cos(Vr))+2*np.cos(Vr)+k*(Ce*(R**2)*dt+t*(R**2-Vr**2)*np.cos(Vr))+((R**2-Vr**2)*np.cos(Vr)))/((R**2)*dt)
    #source = (D_eff*t*((R**2)*Vr[1:-1]*np.cos(Vr[1:-1])+(R**2)*np.sin(Vr[1:-1])-(Vr[1:-1]**3)*np.cos(Vr[1:-1])-5*(Vr[1:-1]**2)*np.sin(Vr[1:-1])+4*Vr[1:-1]*np.cos(Vr[1:-1]))+Vr[1:-1]*(k*(Ce*(R**2)*dt+t*(R**2-Vr[1:-1]**2)*np.cos(Vr[1:-1]))+(R**2-Vr[1:-1]**2)*np.cos(Vr[1:-1])))/((R**2)*dt*Vr[1:-1])
    source = 0
    C=np.zeros_like(Y)
    C[1:-1]=Y[1:-1]+(dt*source)
    C[0]=0
    C[-1]=Ce
    return C

def EL1(analytique,numérique):
    '''Calcul l'erreur L1'''
    errL1=(numérique-analytique)/analytique*abs(numérique-analytique)
    L1=(1/len(analytique))*sum(errL1)
    return L1

def EL2(analytique,numérique):
    '''Calcul l'erreur L2'''
    errL2=(numérique-analytique)/analytique*(abs(numérique-analytique))**2
    L2=(abs((1/len(analytique))*sum(errL2)))**0.5
    return L2
    
def EL3(analytique,numérique):
    '''Calcul erreur L3'''
    errL3=abs(numérique-analytique)
    L3=max(errL3)
    return L3

def Maillage(Ntot):
    '''Calcul les solutions numériques en fonction du nombre de point Ntot'''
    Vr = np.linspace(0,R,Ntot)
    M=MatriceGear(Vr,Ntot)
    Y0 = np.zeros(Ntot) #Condition initiale : C=0 dans tout le pilier
    solution = Euler_implicite_solve(Vr,Vt,M,Y0,MatriceGear,Ntot)
    sol_numérique=solution[-1]
    return Ntot, Vr, sol_numérique

def C(Vr,t):
    '''Donne la solution analytique MMS'''
    return Ce + (t*(1-(Vr**2/R**2))*np.cos(Vr))/dt


## Affichage des résultats
plt.figure("Résultats",figsize=(12,5))
res1,err1=plt.subplot(1,2,1),plt.subplot(1,2,2)
plt.subplots_adjust(left=0.05, right=0.99, bottom=0.06, top=0.94, wspace=0.3, hspace=0.25)

ErreurL2=[]

#Affichage des résultats

for i in range(len(Vn)):
    print(i)
    Ntot,Vr, sol_numérique = Maillage(Vn[i])
    L2=EL2(solCOMSOL[i],sol_numérique)
    ErreurL2.append(L2)
    #Affiche des solutions
    res1.plot(Vr,sol_numérique,".-",label="Ordre 1 : Ntot = {}".format(Ntot))
    res1.plot(Vr,solCOMSOL[i],label="Solution Comsol : Ntot = {}".format(Vr[0]))
    
print(ErreurL2)

#Affichage des erreurs en échelle Log-Log

err1.loglog(R/(Vn-1),ErreurL2,'r1-',label="L2")
    
#Affichage de la solution Analytique

    

err1.loglog([0.1,0.01],[0.01,0.0001],'g*-',label="Référence ordre 2")

err1.legend()
err1.grid()
err1.set_title("Tracé de l'erreur L2 en fonction du pas de temps")
err1.set_xlabel('dt')
err1.set_ylabel('$Erreur L2$')
res1.legend()
res1.grid()
res1.set_title("Tracé de la concentration en fonction du rayon")
res1.set_xlabel('$r$')
res1.set_ylabel('$Concentration$')

err1.invert_xaxis()

plt.show()
##