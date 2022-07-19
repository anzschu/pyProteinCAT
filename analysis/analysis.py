# %%codecell
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show
import seaborn as sns

from settings import MEASUREMENTS, TAXONOMY, ANALYSIS

# %%codecell paths
measurementset = MEASUREMENTS / 'bacteriahalocyanin.csv'
taxonomyset = TAXONOMY / 'halocyanin_bacteria.csv'
analysisfolder = ANALYSIS/'bacteriahalocyanin'

# %%codecell Read csv files and merge
properties = pd.read_csv(measurementset)
taxi = pd.read_csv(taxonomyset)
pt = properties.merge(taxi, on='id')

# %%codecell specify amount of subgroups to view
topfourphyla = pt.value_counts(subset='phylum').index[:4]
toptenclass = pt.value_counts(subset= 'class').index[:10]

# %%code cell plot all possible plots

for x in pt.columns[1:12]:
    for y in pt.columns[1:12]:
        sns.relplot(x= f"{x}", y= f"{y}", hue ='phylum', data = pt[pt.phylum.isin(topfourphyla)])
        plt.show()
        #plt.savefig(analysisfolder/f"{xname}{yname}.png",bbox_inches='tight',dpi=200)

# %%codecell plot one specific plot

sns.relplot(x= 'length', y= 'fNeg', hue ='phylum', data = pt[pt.phylum.isin(topfourphyla)])
plt.show()
# %%codecell plot all plots in one plot
listofcols = ['mwkda', 'fPos', 'fNeg', 'fFatty', 'ncd', 'dpm', 'guy', 'hpm']
sns.pairplot( hue ='phylum', data = pt[pt.phylum.isin(topfourphyla)], vars = listofcols)
plt.show()

# %%codecell plot all plots as kdeplots

for xi in pt.columns[1:12]:
    for yi in pt.columns [1:12]:
        if not xi == yi:
            sns.kdeplot(
                x= xi,
                y= yi,
                hue = 'class', 
                data = pt[pt['class'].isin(toptenclass)],
                legend = False
            )
            plt.show()
            
# %%
