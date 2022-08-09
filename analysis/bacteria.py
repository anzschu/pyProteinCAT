# %%codecell
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show
import seaborn as sns

from settings import MEASUREMENTS, TAXONOMY, ANALYSIS

# %%codecell paths bacteria
measurementset = MEASUREMENTS / 'bacteriahalocyanintrim.csv'
taxonomyset = TAXONOMY / 'halocyanin_bacteria.csv'
analysisfolder = ANALYSIS/'bacteriahalocyanin'

# %%code cell Read csv files and merge for bacteria
properties = pd.read_csv(measurementset)
taxi = pd.read_csv(taxonomyset)
ptb = properties.merge(taxi, on= 'id')

# %%codecell specify amount of subgroups to view

topfourphyla = ptb.value_counts(subset='phylum').index[:4]
topfourclass = ptb.value_counts(subset= 'class').index[:4]
topfourorder = ptb.value_counts(subset='order').index[:4]

#%%
ptb['ncd1000'] = ptb.ncd  * 1000
# %% 
ptb.sort_values('ncd1000', ascending = False)[['ncd1000','phylum','class', 'order', 'family', 'genus']].head(10)

#%%
sns.relplot(
    x = 'ncd1000',
    y = 'length',
    data = ptb
    #[ptb['phylum'].isin(topfourphyla)],
    #hue = 'phylum'
)
# %%
ptb['fPos'].sort_values()
# %%
ptb.iloc[178]
# %%

ptb[ptb.fPos >0.307692 ]
# %%
