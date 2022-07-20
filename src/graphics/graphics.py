import os
import sys
import textwrap
from pathlib import Path
from Bio.PDB import PDBParser
from chimerax.core.commands import run

sys.path.append(Path(__file__).parent.parent.parent.as_posix())

from src.metrics import Builder


def drawDipole(session, proteinfile):
    dvectorfile = generateVectorFile(proteinfile)[0]
    run(session, f"open {proteinfile}")
    run(session, f"open {dvectorfile}")
    # clean up
    os.remove(dvectorfile)

def drawHydrophobe(session, proteinfile):
    hvectorfile = generateVectorFile(proteinfile)[1]
    run(session, f"open {proteinfile}")
    run(session, f"open {hvectorfile}")
    os.remove(hvectorfile)

def generateVectorFile(proteinfile): #maybe pass argument here
    proteinfile = Path(proteinfile)
    dvectorfile = proteinfile.stem + f"{'dipole'}"+ ".bild"
    hvectorfile = proteinfile.stem + f"{'hydrophobe'}" + ".bild"
    parser = PDBParser(QUIET=True, structure_builder=Builder())
    structure = parser.get_structure(proteinfile.stem, proteinfile)
    structure.measure()
    cx, cy, cz = structure.center_of_mass()
    dx, dy, dz = structure.dipolevector()
    hx, hy, hz = structure.hydrophobicvector()
    with open(dvectorfile, "w") as vector:
        body = f"""
                .color 1 1 1
                .sphere {cx} {cy} {cz} 2
                .color 1 0 0
                .arrow {cx} {cy} {cz} {dx} {dy} {dz} 1 3 0.75
                """
        vector.write(textwrap.dedent(body))
    with open(hvectorfile, "w") as vector:
        body = f"""
                .color 1 1 1
                .sphere {cx} {cy} {cz} 2
                .color 0 1 0
                .arrow {cx} {cy} {cz} {hx} {hy} {hz} 1 3 0.75
                """
        vector.write(textwrap.dedent(body))
    return dvectorfile, hvectorfile

drawDipole(session, "data/data1/1ris.pdb")
drawHydrophobe(session, "data/data1/1ris.pdb")
