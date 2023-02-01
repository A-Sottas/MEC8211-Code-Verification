import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import lhsmdu

nb_echan = 100 #Nombre d'échantillons
np.random.seed(0)

#Taille des fibres
moy_df = 12.5
ecart_df = 2.85
serie_df = np.random.normal(moy_df,ecart_df,nb_echan)

#Porosité
moy_poro = 0.9
ecart_poro = 7.5e-3
serie_poro = np.random.normal(moy_poro,ecart_poro,nb_echan)


plt.figure("PDF")
#Tracé des histogrammes normalisés
plot_pdf_df = plt.subplot(2,2,1)
plot_pdf_poro = plt.subplot(2,2,2)
plot_cdf_df = plt.subplot(2,2,3)
plot_cdf_poro = plt.subplot(2,2,4)

plt.subplots_adjust(hspace=0.5)

plot_pdf_df.hist(serie_df,25,density=True)
plot_pdf_poro.hist(serie_poro,25,density=True)

plot_pdf_df.set_xlabel("Diamètre des fibres")
plot_pdf_df.set_ylabel("Nombre")
plot_pdf_poro.set_xlabel("Porosité")
plot_pdf_poro.set_ylabel("Nombre")

#Tracer et fitter les distributions normales correspondantes
xmin_df,xmax_df = moy_df-4*ecart_df, moy_df+4*ecart_df
xmin_poro,xmax_poro = moy_poro-4*ecart_poro, moy_poro+4*ecart_poro
lnspc_df = np.linspace(xmin_df,xmax_df,len(serie_df))
lnspc_poro = np.linspace(xmin_poro,xmax_poro,len(serie_poro))
fit_moy_df,fit_ecart_df = stats.norm.fit(serie_df)
fit_moy_poro,fit_ecart_poro = stats.norm.fit(serie_poro)

#Superposition des PDF

pdf_df = stats.norm.pdf(lnspc_df,fit_moy_df,fit_ecart_df)
label = "Moyenne ="+"{:.2f}".format(fit_moy_df)+'\n'+"Ecart-type ="+"{:.2f}".format(fit_ecart_df)
plot_pdf_df.plot(lnspc_df,pdf_df,label=label)

pdf_poro = stats.norm.pdf(lnspc_poro,fit_moy_poro,fit_ecart_poro)
label = "Moyenne ="+"{:.2f}".format(fit_moy_poro)+'\n'+"Ecart-type ="+"{:.4f}".format(fit_ecart_poro)
plot_pdf_poro.plot(lnspc_poro,pdf_poro,label=label)

#Tracé des CDF
plot_cdf_df.hist(serie_df,20,cumulative=True,density=True)
plot_cdf_poro.hist(serie_poro,20,cumulative=True,density=True)

plot_cdf_df.set_xlabel("Diamètre des fibres")
plot_cdf_df.set_ylabel("Probabilité df < Diamètre des fibres")
plot_cdf_poro.set_xlabel("Porosité")
plot_cdf_poro.set_ylabel("Probabilité p < porosité")

cdf_df = stats.norm.cdf(lnspc_df,fit_moy_df,fit_ecart_df)
cdf_poro = stats.norm.cdf(lnspc_poro,fit_moy_poro,fit_ecart_poro)

plot_cdf_df.plot(lnspc_df,cdf_df,label="Norm")
plot_cdf_poro.plot(lnspc_poro,cdf_poro,label="Norm")

#Légende et plot

plot_pdf_poro.set_title("PDF df")
plot_pdf_df.set_title("PDF porosité")
plot_cdf_poro.set_title("CDF df")
plot_cdf_df.set_title("CDF df")

plot_pdf_poro.legend()
plot_pdf_df.legend()
plot_cdf_poro.legend()
plot_cdf_df.legend()

plt.show()

#MonteCarlo

MC = lhsmdu.createRandomStandardUniformMatrix(2,10)

MC_df = stats.norm.ppf(MC[0],fit_moy_df,fit_ecart_df)
MC_poro = stats.norm.ppf(MC[1],fit_moy_poro,fit_ecart_poro)

print(MC_df)
print(MC_poro)