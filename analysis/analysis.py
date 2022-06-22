import pandas as pd

from settings import MEASUREMENTS, TAXONOMY

measurementset = MEASUREMENTS / 'bacteriahalocyanin.csv'

taxonomyset = TAXONOMY / 'halocyanin_bacteria.csv'

properties = pd.read_csv(measurementset)
taxi = pd.read_csv(taxonomyset)

print(properties.merge(taxi, on='id', ))
