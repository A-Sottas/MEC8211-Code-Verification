import matplotlib.pyplot as plt

X = 1
Y = -61.1939
Erreur1 = 38.521
X2 = 2
Erreur2 = 44.038

plt.figure("Estimation erreur du modèle")
plt.plot(X,Y,'b.',label="Erreur avec loi Gausienne")
plt.errorbar(X, Y, xerr = 0, yerr = Erreur1 , fmt = 'none', capsize = 10, ecolor = 'red', zorder = 1)

plt.plot(X2,Y,'g.',label="Erreur avec loi Log-Normale")
plt.errorbar(X2, Y, xerr = 0, yerr = Erreur2 , fmt = 'none', capsize = 10, ecolor = 'red', zorder = 1)

plt.xlim(0,3)
plt.ylim(-120,120)
plt.grid(axis = "y")
plt.xticks(color='w')
plt.legend()
plt.show()
