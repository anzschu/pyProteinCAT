import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show
import seaborn as sns

from settings import MEASUREMENTS, TAXONOMY, ANALYSIS

measurementset = MEASUREMENTS / 'bacteriahalocyanin.csv'
taxonomyset = TAXONOMY / 'halocyanin_bacteria.csv'
analysisfolder = ANALYSIS/'bacteriahalocyanin'

# Read csv files and merge
properties = pd.read_csv(measurementset)
taxi = pd.read_csv(taxonomyset)
pt = properties.merge(taxi, on='id')

# slice pt for only toptenphyla 
toptenphyla = pt.value_counts(subset='phylum').index[:10]
# drop rows if not in list
#print(pt.phylum.values)
rowstoremove = []
#for i in pt.iterrows():
for i in pt.iterrows():
    if pt.phylum.values[i] in toptenphyla:
        print('yes')
    #else:
    #    rowstoremove.append()

#for x in pt.columns[1:12]:
#    xname = str(x)
#    for y in pt.columns[1:12]:
#        yname = str(y)
#        sns.relplot(x= f"{xname}", y= f"{yname}", hue ='phylum', data = pt)
#        plt.savefig(analysisfolder/f"{xname}{yname}.png",bbox_inches='tight',dpi=200)