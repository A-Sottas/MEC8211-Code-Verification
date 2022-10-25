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

Vn=np.array([5,50,500]) #Choix du nombre de points dans le maillage

dt = 365*3600*24 #Base de temps : 1 an
Vt = np.arange(0,1e10,dt) #Vecteur des temps t

##Construction des matrices

def B1(alpha,V,dr):
    '''Coefficients B de la matrice avec schémas à l'ordre 1'''
    T=np.ones_like(V)
    for i in range(len(V)-1):
        T[i]=V[i+1]
    return 1+dt*k+alpha*((2/dr**2)+1/(2*T*dr))

def D1(alpha,V,dr):
    '''Coefficients D de la matrice avec schémas à l'ordre 1'''
    T=np.ones_like(V)
    T[1:]=V[1:]
    return -alpha*((1/dr**2)+1/(2*T*dr))

def B2(alpha,V,dr):
    '''Coefficients B de la matrice avec schémas à l'ordre 2'''
    T=np.ones_like(V)
    for i in range(len(V)-1):
        T[i]=V[i+1]
    return 1+alpha*((2/dr**2)+1/(2*T*dr))

def D2(alpha,V,dr):
    '''Coefficients D de la matrice avec schémas à l'ordre 2'''
    T=np.ones_like(V)
    T[1:]=V[1:]
    return -alpha*((1/dr**2)+1/(2*T*dr))
    
def MatriceSansSource(V,Ntot):
    ''' Création de la matrice sans schéma de Gear'''
    dr = R/(Ntot-1) #Définition du dr
    alpha = dt*D_eff
    A = -alpha/dr**2 #Coefficient A indépendant de r
    M = A*np.eye(Ntot,k=-1) + B1(alpha,V,dr)*np.eye(Ntot) + D1(alpha,V,dr)*np.eye(Ntot,k=1) #Construction de la matrice
    M[0],M[-1] = 0,0
    M[0,0],M[0,1]=-1,1 #Schéma d'ordre 1
    M[-1,-1] = 1
    return M

def MatriceGear(V,Ntot):
    ''' Création de la matrice avec schéma de Gear'''
    dr = R/(Ntot-1)
    alpha = dt*D_eff
    A = -alpha/dr**2
    M = A*np.eye(Ntot,k=-1) + B1(alpha,V,dr)*np.eye(Ntot) + D1(alpha,V,dr)*np.eye(Ntot,k=1)
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
    M,M2=MatriceSansSource(Vr,Ntot),MatriceGear(Vr,Ntot)
    Y0 = np.zeros(Ntot) #Condition initiale : C=0 dans tout le pilier
    solution = Euler_implicite_solve(Vr,Vt,M,Y0,MatriceSansSource,Ntot)
    solution2 = Euler_implicite_solve(Vr,Vt,M2,Y0,MatriceGear,Ntot)
    sol_analytique=C(Vr)
    sol_numérique=solution[-1]
    sol_numérique2=solution2[-1]
    L1=EL1(sol_analytique,sol_numérique)
    L11=EL1(sol_analytique,sol_numérique2)
    L2=EL2(sol_analytique,sol_numérique)
    L22=EL2(sol_analytique,sol_numérique2)
    L3=EL3(sol_analytique,sol_numérique)
    L33=EL3(sol_analytique,sol_numérique2)
    return Ntot,abs(L1), abs(L11), abs(L2), abs(L22), abs(L3), abs(L33), Vr, sol_analytique, sol_numérique,sol_numérique2

## Affichage des résultats
res1,err1=plt.subplot(2,2,1),plt.subplot(2,2,2)
res2,err2=plt.subplot(2,2,3),plt.subplot(2,2,4)
plt.subplots_adjust(left=0.05, right=0.99, bottom=0.06, top=0.94, wspace=0.3, hspace=0.25)

ErreurL1=[]
ErreurL2=[]
ErreurL3=[]
ErreurL11=[]
ErreurL22=[]
ErreurL33=[]

#Affichage des résultats
for N in Vn:
    Ntot, L1, L11, L2, L22, L3, L33,Vr,sol_analytique, sol_numérique,sol_numérique2 = Maillage(N)
    ErreurL1.append(L1)
    ErreurL11.append(L11)
    ErreurL2.append(L2)
    ErreurL22.append(L22)
    ErreurL3.append(L3)
    ErreurL33.append(L33)
    #Affiche des solutions
    res1.plot(Vr,sol_numérique,".-",label="Ordre 1 : Ntot = {}".format(Ntot))
    res2.plot(Vr,sol_numérique2,".--",label="Ordre 2 : Ntot = {}".format(Ntot))

#Affichage des erreurs en échelle Log-Log
err1.loglog(R/(Vn-1),ErreurL1,'b1-',label="L1")
err2.loglog(R/(Vn-1),ErreurL11,'b2-',label="L2")
err1.loglog(R/(Vn-1),ErreurL2,'r1-',label="L2")
err2.loglog(R/(Vn-1),ErreurL22,'r2-',label="L2 Gear")
err1.loglog(R/(Vn-1),ErreurL3,'k1-',label="Linf")
err2.loglog(R/(Vn-1),ErreurL33,'k2-',label="Linf Gear")
    
#Affichage de la solution Analytique
res1.plot(Vr,sol_analytique,"-",label="Solution analytique")
res2.plot(Vr,sol_analytique,"-",label="Solution analytique")

err1.loglog([0.1,0.01],[0.01,0.0001],'g*-',label="Référence ordre 2")
err2.loglog([0.1,0.01],[0.01,0.0001],'g*-',label="Référence ordre 2")

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
err2.legend()
err2.grid()
err2.set_title("Tracé de l'erreur L2 en fonction du maillage")
err2.set_xlabel('dr')
err2.set_ylabel('$Erreur L2$')
res2.legend()
res2.grid()
res2.set_title("Tracé de la concentration en fonction de la distance (méthode Gear)")
res2.set_xlabel('$r$')
res2.set_ylabel('$Concentration$')

err1.invert_xaxis()
err2.invert_xaxis()

plt.show()
##
