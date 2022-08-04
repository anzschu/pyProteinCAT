import pandas as pd
from src.metrics import Builder
from pathlib import Path
from Bio.PDB import PDBParser

dataset = Path('data/models')



def read_data(fname):
    '''
    Read file and build structure using Builder from metrics module.
    '''
    parser = PDBParser(QUIET=1, structure_builder=Builder())
    s = parser.get_structure( fname.stem , fname)
    return s

if __name__ == '__main__':
    results = []
    for pdbfile in dataset.iterdir():
        s = read_data(pdbfile)
        s.measure()
        results.append(s.serializer())
    proteindata = pd.DataFrame(results)
    #print(proteindata)