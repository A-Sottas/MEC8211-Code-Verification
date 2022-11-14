import lhsmdu
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

moy = 12.5
ecart = 2.85
serie = 

xmin,xmax = moy-4*ecart, moy+4*ecart

lnspc_df = np.linspace(xmin,xmax,len(serie))

fit_moy_df,fit_ecart_df = stats.norm.fit(serie)

