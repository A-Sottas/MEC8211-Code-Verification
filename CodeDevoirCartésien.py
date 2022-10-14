#Auteur : Mathieu Goureau/Arona Sottas
#mathieu.goureau@polymtl.ca

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

#Définition des constantes

S =  1e-8 #Terme source constant
D_eff = 1e-10 #Coefficient de diffusion constant
Ce = 10 #Concentration extérieure
D = 1 #Diamètre du pilier
R = D/2 #Rayon du pilier

Vn=np.array([5,50,500]) #Choix du nombre de points dans le maillage

dt = 365*24*3600 #Base de temps : 1 an
Vt = np.arange(0,5e9,dt) #Vecteur des temps t
    
def Matrice(V,Ntot):
    ''' Création de la matricec sans schéma de Gear'''
    dr = R/(Ntot-1)
    alpha = dt*D_eff/(dr**2)
    A = [-alpha]*(Ntot-1)
    B = [1+2*alpha]*Ntot
    D = [-alpha]*(Ntot-1)
    M = np.diag(A,k=-1) + np.diag(B) + np.diag(D,k=1)
    M[0],M[-1] = 0,0
    M[0,0],M[0,1]=-1,1
    M[-1,-1] = 1
    return M

def MatriceGear(V,Ntot):
    ''' Création de la matrice avec schéma de Gear'''
    dr = R/(Ntot-1)
    alpha = dt*D_eff/(dr**2)
    A = [-alpha]*(Ntot-1)
    B = [1+2*alpha]*Ntot
    D = [-alpha]*(Ntot-1)
    M = np.diag(A,k=-1) + np.diag(B) + np.diag(D,k=1)
    M = np.diag(A,k=-1) + np.diag(B) + np.diag(D,k=1) #Matrice tridiagonale
    M[0],M[-1]=0,0
    M[0,0],M[0,1],M[0,2]=-3,4,-1
    M[-1,-1] = 1
    return M

def VecC(Y,Ntot):
    '''Renvoie le vecteur colonne C'''
    source=-dt*S*np.ones(Ntot)
    C=Y+source
    C[0] = 0
    C[-1] = Ce
    return C

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

def C(Vr): #Solution analytique
    return 0.5*(S/D_eff)*(R**2)*((Vr**2/R**2)-1)+Ce


def EL1(analytique,numérique): #Calcul erreur L1
    errL1=(numérique-analytique)/analytique*abs(numérique-analytique)
    L1=(1/len(analytique))*sum(errL1)
    return errL1, L1

def EL2(analytique,numérique): #Calcul erreur L2
    errL2=(numérique-analytique)/analytique*(abs(numérique-analytique))**2
    L2=((1/len(analytique))*sum(errL2))*0.5
    return errL2, L2
    
def EL3(analytique,numérique): #Calcul erreur L3
    errL3=abs(numérique-analytique)
    L3=max(errL3)
    return errL3, L3

def Maillage(Ntot):
    Vr = np.linspace(0,R,Ntot)
    M,M2=Matrice(Vr,Ntot),MatriceGear(Vr,Ntot)
    Y0 = np.zeros(Ntot)
    solution = Euler_implicite_solve(Vr,Vt,M,Y0,Matrice,Ntot)
    solution2 = Euler_implicite_solve(Vr,Vt,M2,Y0,MatriceGear,Ntot)
    sol_analytique=C(Vr)
    sol_numérique=solution[-1]
    sol_numérique2=solution2[-1]
    errL2,L2=EL2(sol_analytique,sol_numérique)
    errL22,L22=EL2(sol_analytique,sol_numérique2)
    return Ntot, L2, L22, Vr, sol_analytique, sol_numérique,sol_numérique2

## Affichage des résultats
res,err=plt.subplot(1,2,1),plt.subplot(1,2,2)
plt.subplots_adjust(left=0.05, right=0.99, bottom=0.06, top=0.94, wspace=0.3)

#Affichage des résultats
for N in Vn:
    Ntot, L2, L22,Vr,sol_analytique, sol_numérique,sol_numérique2 = Maillage(N)
    #Affichage des erreurs L2 en échelle Log-Log
    err.loglog(Ntot,abs(L2),'b1-',label="Sol1")
    err.loglog(Ntot,abs(L22),'r2-',label="Sol2")
    #Affiche des solutions
    res.plot(Vr,sol_numérique,".-",label="Ordre 1 : Ntot = {}".format(Ntot))
    res.plot(Vr,sol_numérique2,".--",label="Ordre 2 : Ntot = {}".format(Ntot))
    
#Affichage Solution Analytique
res.plot(Vr,sol_analytique,"-",label="Solution analytique")

err.loglog([10,100],[0.001,0.00001],'g*-') #Référence ordre 2

err.legend()
err.grid()
res.legend()
res.grid()

plt.show()