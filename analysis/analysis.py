# %%codecell
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show
import seaborn as sns

from settings import MEASUREMENTS, TAXONOMY, ANALYSIS

# %%codecell paths archaea
measurementseta = MEASUREMENTS / 'archaeahalocyanintrim.csv'
taxonomyseta = TAXONOMY / 'halocyanin_archaea.csv'
analysisfoldera = ANALYSIS/'archaeahalocyanin'

# %%codecell paths bacteria
measurementsetb = MEASUREMENTS / 'bacteriahalocyanintrim.csv'
taxonomysetb = TAXONOMY / 'halocyanin_bacteria.csv'
analysisfolderb = ANALYSIS/'bacteriahalocyanin'

analysisfolderab = ANALYSIS/'archaeabacteria'
# %%codecell Read csv files and merge for archaea
propertiesa = pd.read_csv(measurementseta)
taxia = pd.read_csv(taxonomyseta)
pta = propertiesa.merge(taxia, on='id')

# %%code cell Read csv files and merge for bacteria
propertiesb = pd.read_csv(measurementsetb)
taxib = pd.read_csv(taxonomysetb)
ptb = propertiesb.merge(taxib, on= 'id')

# %%codecell specify amount of subgroups to view
topfourphylaa = pta.value_counts(subset='phylum').index[:4]
topfourclassa = pta.value_counts(subset= 'class').index[:4]
topfourordera = ptb.value_counts(subset='order').index[:4]

topfourphylab = ptb.value_counts(subset='phylum').index[:4]
topfourclassb = ptb.value_counts(subset= 'class').index[:4]
topfourorderb = ptb.value_counts(subset='order').index[:4]


# %%code Plot all possible plots

for x in pta.columns[1:12]:
    for y in pta.columns[1:12]:
        sns.relplot(x= f"{x}", y= f"{y}", hue ='kingdom', data = ptb)
        plt.show()
        #plt.savefig(analysisfolder/f"{xname}{yname}.png",bbox_inches='tight',dpi=200)

# %%code cell Plot both bacteria and archaea in one plot, all plots

d1 = pta
#[pta['phylum'].isin(topfourphylaa)]
d2 = ptb
#[ptb['phylum'].isin(topfourphylab)]
frames = [d1, d2]
da = pd.concat(frames).reset_index(drop=True)

for xi in da.columns[1:12]:
    for yi in da.columns [1:12]:
        if not xi == yi:
            sns.relplot(
                x = xi,
                y = yi,
                hue = 'kingdom',
                data = da,
                palette = 'coolwarm'
                )
            plt.show()

# %% Plot archaea and bacteria in one plot, plots only one plot
d1 = pta[pta['class'].isin(topfourclassa)]
#['phylum'].isin(topfourphylaa)
d2 = ptb[ptb['class'].isin(topfourclassb)]
#['phylum'].isin(topfourphylab)
frames = [d1, d2]
da = pd.concat(frames).reset_index(drop=True)
da['ncd1000'] = da.ncd * 1000
da['fc'] = da.fPos + da.fNeg
#print(da)

sns.relplot(
    x = 'ncd',
    y = 'mwkda',
    hue = 'kingdom',
    data = da,
    palette = 'coolwarm_r',
    )
#plt.savefig(analysisfolderab/'ncdmwkdaphylumrelplot.png', bbox_inches='tight',dpi=200)