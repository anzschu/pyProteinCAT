import os
import pandas as pd
import matplotlib as mp
import matplotlib.pyplot as plt
from matplotlib.pyplot import plot, show
from pathlib import Path
from Bio.PDB import PDBParser
from Bio.PDB.PDBIO import PDBIO

from src.metrics import Builder
from settings import DATA, TRIM
from src.metrics.metrics import ModStructure


dataset = DATA / 'archaeahalocyanin'
trimset = TRIM 

def read_data(fname):
    parser = PDBParser(QUIET=1, structure_builder=Builder(is_AF=True))
    s = parser.get_structure(fname.stem, fname)
    return s

def save_pdb(fname):
    s = read_data(fname)
    io = PDBIO()
    io.set_structure(s)
    with open(TRIM / f"{fname.stem}_trim.pdb", mode='w') as f:
        io.save(f)

def findingthreshold(fname):
    s= read_data(fname)
    y_pLDDT = []
    x_id= []
    for residue in s.get_residues():
        for atom in residue.get_atoms():
            y_pLDDT.append(atom.get_bfactor())
            x_id.append(residue.id[1])
            break
    print(y_pLDDT)
    plt.figure()
    plt.plot(x_id, y_pLDDT, linestyle="-", marker = "o")
    plt.xlabel("Residue")
    plt.ylabel("pLDDT")
    plt.title(f"{fname.stem}")
    plt.savefig(TRIM/f"{fname.stem}2.png", bbox_inches='tight',dpi=200)
   
if __name__ =='__main__':
    #firstten = 0
    #for data in dataset.iterdir():
    #    firstten += 1
    #    save_pdb(data)
    #    findingthreshold(data)
    #    if firstten == 10:
    #        break
    # s.measure()
    data = dataset/'A0A0F7IFE1.pdb'
    save_pdb(data)
    findingthreshold(data)

