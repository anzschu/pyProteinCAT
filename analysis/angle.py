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
properties = properties[properties.length >= 80]
properties.angle = properties.angle.astype(float)
taxi = pd.read_csv(taxonomyset)
pta = properties.merge(taxi, on='id')

# %%codecell specify amount of subgroups to view
topfourphyla = pta.value_counts(subset='phylum').index[:4]
topsixclass = pta.value_counts(subset= 'class').index[:6]
topfourorder = pta.value_counts(subset='order').index[:4]

# %%codecell

def anglemap(row):
    if row.angle < 90:
        return 'a<90'
    else:
        return 'a>90'

pta['angletype'] = pta.apply(anglemap, axis =1)

# %%codecell
pta['logdpm'] = np.log(pta.dpm + 1)
pta['loghpm'] = np.log(pta.hpm + 1)
fig = plt.figure()
ax = fig.subplots(1, 1)
pta.plot.scatter(
    'logdpm',
    'loghpm',
    c='angle',
    cmap='vlag',
    ax=ax
)
fig.savefig('analysis/archaeahalocyanin/dpmhpmangle.png', dpi=200, bbox_inches='tight')
# %%codecell
