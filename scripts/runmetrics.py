import pandas as pd
from src.metrics import Builder
from pathlib import Path
from Bio.PDB import PDBParser

dataset = Path('data/models')



def read_data(fname):
    parser = PDBParser(QUIET=1, structure_builder=Builder())
    s = parser.get_structure( fname.stem , fname)
    return s

results = []
for pdbfile in dataset.iterdir():
    s = read_data(pdbfile)
    s.measure()
    results.append(s.serializer())

proteindata = pd.DataFrame(results)
print(proteindata)