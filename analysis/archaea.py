# %%codecell
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show
import seaborn as sns

from settings import MEASUREMENTS, TAXONOMY, ANALYSIS

# %%codecell paths archaea
measurementset = MEASUREMENTS / 'archaeahalocyanintrim.csv'
taxonomyset = TAXONOMY / 'halocyanin_archaea.csv'
analysisfolder = ANALYSIS/'archaeahalocyanin'

# %%codecell Read csv files and merge for archaea
properties = pd.read_csv(measurementset)
taxi = pd.read_csv(taxonomyset)
pta = properties.merge(taxi, on='id')

# %%codecell specify amount of subgroups to view
topfourphyla = pta.value_counts(subset='phylum').index[:4]
topsixclass = pta.value_counts(subset= 'class').index[:6]
topfourorder = pta.value_counts(subset='order').index[:4]


#%% 

pta.value_counts('class')
# %%
sns.kdeplot(
        x= 'ncd',
        y= 'mwkda',
        hue= 'king',
        data = pta[pta['class'].isin(topsixclass)],
        palette = 'hot'
)
#plt.savefig(analysisfolder/'ncdmwkdakdeclass.png', bbox_inches='tight',dpi=200)
# %%codcell mwkda ncd
pta['logdpm'] = np.log(pta.dpm + 1)
pta['ncd1000'] = pta.ncd  * 1000
sns.kdeplot(
        x='ncd1000',
        y='logdpm',
        hue = 'class', 
        data = pta[pta['class'].isin(topsixclass)],
        palette = 'hls')
# %%

sns.relplot(
        x='ncd',
        y='dpm',

        hue = 'class', 
        data = pta[pta['class'].isin(topsixclass)],
        palette = 'hls')

# %%
#pta['logdpm'] = np.log(pta.dpm + 1)
pta['loghpm'] = np.log(pta.hpm + 1)
sns.regplot(
    x='logdpm',
    y='loghpm',
    #hue = 'phylum', 
    data = pta[pta['phylum'].isin(topfourphyla)],
    #palette = 'hls',
    #legend = False,
    
)
# %%
pta.corr().loc['dpm', 'hpm']
# %%
