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

Vn=np.array([5,50,500,1000]) #Choix du nombre de points dans le maillage

#dt = 365*3600*24 #Base de temps : 1 an
dt = 1e6
Vt = np.arange(0,1e9,dt) #Vecteur des temps t

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
    source = (D_eff*t*((R**2)*Vr[1:-1]*np.cos(Vr[1:-1])+(R**2)*np.sin(Vr[1:-1])-(Vr[1:-1]**3)*np.cos(Vr[1:-1])-5*(Vr[1:-1]**2)*np.sin(Vr[1:-1])+4*Vr[1:-1]*np.cos(Vr[1:-1]))+Vr[1:-1]*(k*(Ce*(R**2)*dt+t*(R**2-Vr[1:-1]**2)*np.cos(Vr[1:-1]))+(R**2-Vr[1:-1]**2)*np.cos(Vr[1:-1])))/((R**2)*dt*Vr[1:-1])
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
    Vr = np.linspace(0,R,Ntot) #Définition du Maillage
    M=MatriceGear(Vr,Ntot)
    Y0 = np.zeros(Ntot) #Condition initiale : C=0 dans tout le pilier
    solution = Euler_implicite_solve(Vr,Vt,M,Y0,MatriceGear,Ntot)
    sol_numérique=solution[-1]
    return Ntot, Vr, sol_numérique

def C(Vr,t):
    '''Donne la solution analytique MMS'''
    return Ce + (t*(1-(Vr**2/R**2))*np.cos(Vr))/dt


## Affichage des résultats
plt.figure("Résultats Espace",figsize=(12,5))
res1,err1=plt.subplot(1,2,1),plt.subplot(1,2,2)
plt.subplots_adjust(left=0.05, right=0.99, bottom=0.06, top=0.94, wspace=0.3, hspace=0.25)

ErreurL2=[]

#Affichage des résultats
for N in Vn:
    Ntot,Vr, sol_numérique = Maillage(N)
    sol_analytique = C(Vr,Vt[-1])
    L2=EL2(sol_analytique,sol_numérique)
    ErreurL2.append(L2)
    #Affiche des solutions
    res1.plot(Vr,sol_numérique,".-",label="Ordre 1 : Ntot = {}".format(Ntot))
    res1.plot(Vr,sol_analytique,label="solution analytique : Ntot = {}".format(Ntot))

print(ErreurL2)

#Affichage des erreurs en échelle Log-Log

err1.loglog(R/(Vn-1),ErreurL2,'r1-',label="L2")
    
#Affichage de la solution Analytique

    

err1.loglog([0.1,0.01],[0.01,0.0001],'g*-',label="Référence ordre 2")

err1.legend()
err1.grid()
err1.set_title("Tracé de l'erreur L1 en fonction du maillage")
err1.set_xlabel('dr')
err1.set_ylabel('$Erreur L1$')
res1.legend()
res1.grid()
res1.set_title("Tracé de la concentration en fonction de la distance")
res1.set_xlabel('$r$')
res1.set_ylabel('$Concentration$')

err1.invert_xaxis()

plt.show()
##
