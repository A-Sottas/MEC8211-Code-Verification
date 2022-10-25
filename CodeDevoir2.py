#Auteur : Mathieu Goureau/ Arona Sottas
#mathieu.goureau@polymtl.ca

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

##Définition des constantes

S =  1e-8 #Terme source constant
D_eff = 1e-10 #Coefficient de diffusion constant
Ce = 10 #Concentration extérieure
D = 1 #Diamètre du pilier
R = D/2 #Rayon du pilier

Vn=np.array([5,50]) #Choix du nombre de points dans le maillage

dt = 365*3600*24 #Base de temps : 1 an
Vt = np.arange(0,1e10,dt) #Vecteur des temps t

##Construction des matrices

def B(alpha,V,dr):
    '''Coefficients B de la matrice avec schémas à l'ordre 2'''
    T=np.ones_like(V)
    for i in range(len(V)-1):
        T[i]=V[i+1]
    return 1+alpha*((2/dr**2)+1/(2*T*dr))

def D(alpha,V,dr):
    '''Coefficients D de la matrice avec schémas à l'ordre 2'''
    T=np.ones_like(V)
    T[1:]=V[1:]
    return -alpha*((1/dr**2)+1/(2*T*dr))
    
def MatriceGear(V,Ntot):
    ''' Création de la matrice avec schéma de Gear'''
    dr = R/(Ntot-1)
    alpha = dt*D_eff
    A = -alpha/dr**2
    M = A*np.eye(Ntot,k=-1) + B(alpha,V,dr)*np.eye(Ntot) + D(alpha,V,dr)*np.eye(Ntot,k=1)
    M[0],M[-1]=0,0
    M[0,0],M[0,1],M[0,2]=-3,4,-1 #Schéma de Gear
    M[-1,-1] = 1
    return M

def Euler_implicite_solve(Vr,Vt,M,Y0,Ntot):
    """ Résolution du système linéaire à chaque t selon la méthode Euler implicite"""
    solu=[] 
    Y=Y0 #Iniatialisation de la solution
    solu.append(Y0.tolist())
    for t in Vt: #Résolution à chaque t du système AX=B
        C=VecC(Y,Ntot) #Vecteur B
        Y=np.linalg.solve(M,C) #Résolution du système
        solu.append(Y.tolist())
    return np.array(solu)

def VecC(Y,Ntot):
    '''Renvoie le vecteur colonne C'''
    source=-dt*S*np.ones(Ntot) #Terme source
    C=Y+source
    C[0] = 0
    C[-1] = Ce
    return C

def C(Vr):
    '''Donne la solution analytique'''
    return 0.25*(S/D_eff)*R**2*(Vr**2/R**2-1)+Ce


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
    Vr = np.linspace(0,R,Ntot) #Définition du Maillage
    M=MatriceGear(Vr,Ntot)
    Y0 = np.zeros(Ntot) #Condition initiale : C=0 dans tout le pilier
    solution = Euler_implicite_solve(Vr,Vt,M,Y0,Ntot)
    sol_analytique=C(Vr)
    sol_numérique=solution[-1]
    L1=EL1(sol_analytique,sol_numérique)
    L2=EL2(sol_analytique,sol_numérique)
    L3=EL3(sol_analytique,sol_numérique)
    return Ntot,abs(L1),abs(L2),abs(L3),Vr, sol_analytique, sol_numérique

## Affichage des résultats
res1,err1=plt.subplot(2,2,1),plt.subplot(2,2,2)
plt.subplots_adjust(left=0.05, right=0.99, bottom=0.06, top=0.94, wspace=0.3, hspace=0.25)

ErreurL1=[]
ErreurL2=[]
ErreurL3=[]

#Affichage des résultats
for N in Vn:
    Ntot, L1, L2, L3,Vr,sol_analytique, sol_numérique = Maillage(N)
    ErreurL1.append(L1)
    ErreurL2.append(L2)
    ErreurL3.append(L3)
    #Affiche des solutions
    res1.plot(Vr,sol_numérique,".-",label="Ordre 1 : Ntot = {}".format(Ntot))

#Affichage des erreurs en échelle Log-Log
err1.loglog(R/(Vn-1),ErreurL1,'b1-',label="L1")
err1.loglog(R/(Vn-1),ErreurL2,'r1-',label="L2")
err1.loglog(R/(Vn-1),ErreurL3,'k1-',label="Linf")
    
#Affichage de la solution Analytique
res1.plot(Vr,sol_analytique,"-",label="Solution analytique")

err1.loglog([0.1,0.01],[0.01,0.0001],'g*-',label="Référence ordre 2")

err1.legend()
err1.grid()
err1.set_title("Tracé de l'erreur L1 en fonction du maillage")
err1.set_xlabel('dr')
err1.set_ylabel('$Erreur L1$')
res1.legend()
res1.grid()
res1.set_title("Tracé de la concentration en fonction de la distance (méthode dérivée première)")
res1.set_xlabel('$r$')
res1.set_ylabel('$Concentration$')

err1.invert_xaxis()

plt.show()
##