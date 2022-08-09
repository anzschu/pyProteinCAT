import pandas as pd
from pathlib import Path
from Bio.PDB import PDBParser

from datetime import datetime

from src.metrics import Builder
from settings import DATA, MEASUREMENTS


dataset = DATA / 'bacteriahalocyanin'
measurementset = MEASUREMENTS / 'bacteriahalocyanintrim.csv'

def read_data(fname):
    parser = PDBParser(QUIET=1, structure_builder=Builder(is_AF=True))
    s = parser.get_structure( fname.stem , fname)
    return s
#start = datetime.now()
results = []
for pdbfile in dataset.iterdir():
    #print(pdbfile)
    s = read_data(pdbfile)
    s.measure()
    results.append(s.serializer())

proteindata = pd.DataFrame(results)
proteindata.to_csv( measurementset, index = False)

#end = datetime.now()
#print(end-start)