import os
import sys
import textwrap
from pathlib import Path
from Bio.PDB import PDBParser
from chimerax.core.commands import run

sys.path.append(Path(__file__).parent.parent.parent.as_posix())

from src.metrics import Builder


def drawDipole(session, proteinfile):
    vectorfile = generateDipoleFile(proteinfile)
    run(session, f"open {proteinfile}")
    run(session, f"open {vectorfile}")
    # clean up
    os.remove(vectorfile)


def generateDipoleFile(proteinfile):
    proteinfile = Path(proteinfile)
    vectorfile = proteinfile.stem + ".bild"
    parser = PDBParser(QUIET=True, structure_builder=Builder())
    structure = parser.get_structure(proteinfile.stem, proteinfile)
    cx, cy, cz = structure.center_of_mass()
    dx, dy, dz = structure.dipolevector()
    with open(vectorfile, "w") as vector:
        body = f"""
                .sphere {cx} {cy} {cz} 2
                .color 1 0 0
                .arrow {cx} {cy} {cz} {dx} {dy} {dz} 1 3 0.75
                .color 1 1 0
                """
        vector.write(textwrap.dedent(body))
    return vectorfile


drawDipole(session, "data/data1/1ris.pdb")
