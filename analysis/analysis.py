#%% codecell
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show

from settings import MEASUREMENTS, TAXONOMY


#%% codecell
measurementset = MEASUREMENTS / 'bacteriahalocyanin.csv'

taxonomyset = TAXONOMY / 'halocyanin_bacteria.csv'

properties = pd.read_csv(measurementset)
taxi = pd.read_csv(taxonomyset)
pt = properties.merge(taxi, on='id')


#%%codecell

plt.figure()
plt.plot(pt.dpm, pt.hpm)

# %%
