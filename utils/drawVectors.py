'''
Visualizes dipole and hydrophobic vectors with respect to the protein's structure.

Example: 
    $ chimera --script "drawVectors.py 1abc.pdb"

'''
import os
import sys
import textwrap
import argparse
import numpy as np
from pathlib import Path
from Bio.PDB import PDBParser
from chimerax.core.commands import run
from chimerax.markers import MarkerSet
from chimerax.struct_measure.tool import StructMeasureTool
from chimerax.label import label_create

sys.path.append(Path(__file__).parent.parent.as_posix())

from src.metrics import Builder  # noqa
from settings import TRIM

def labels(session, label: str):
    '''
    Creates labels for hydrophobic vector, dipole vector and protein name.
    '''
    label_create(
        session,
        name='hv',
        text='hydrophobic vector',
        color=(255, 140, 0, 255),
        xpos=0.5,
        ypos=0.9,
        size=21
    )
    label_create(
        session,
        name='dv',
        text='dipole vector',
        color=(0, 0, 255, 255),
        xpos=0.1,
        ypos=0.9,
        size=21
    )
    label_create(
        session,
        name=label,
        text=label,
        xpos=0.45,
        ypos=0.1,
        size=21
    )

def drawArc(session, structure):
    '''
    Draws arc between both vectors.
    '''
    com, dpv, hpv = createcoordinates(structure)
    dx, dy, dz = com + dpv / np.linalg.norm(dpv) * 10
    hx, hy, hz = com + hpv / np.linalg.norm(hpv) * 10
    body = f"""
            .color 0 0.5 0.5
            .transparency 1
            .dot {dx} {dy} {dz}
            .transparency 0
            .draw {hx} {hy} {hz}
            """
    with open('_arc.bild', "w") as arcfile:
        arcfile.write(textwrap.dedent(body))

    run(session, "open _arc.bild")
    # clean up
    os.remove('_arc.bild')

def angles(session, structure):
    '''
    Defines markers for the center of mass and the tip of the hydrophobic and
    dipole vector. Calculates the angle between the both vectors. Source for marker set:
    https://rbvi.github.io/chimerax-recipes/mark_blobs/mark_blobs.html
    '''
    c, d, h = createcoordinates(structure)
    marker_set = MarkerSet(session, name="markersforangle")
    marker_set.create_marker(c + d, (0, 0, 255, 255), 0.5)
    marker_set.create_marker(c, (255, 255, 255, 255), 0.5)
    marker_set.create_marker(c + h, (255, 140, 0, 255), 0.5)
    session.models.add([marker_set])
    run(session, "select add #4")
    anglemaker = StructMeasureTool(session)
    anglemaker._create_angle()

def drawDipole(session, structure):
    '''
    Creates dipole vector in ChimeraX session.
    '''
    dvectorfile = generateVectorFile(structure)[0]
    run(session, f"open {dvectorfile}")
    # clean up
    os.remove(dvectorfile)

def drawHydrophobe(session, structure):
    '''
    Creates hydrophobic vector in ChimeraX session.
    '''
    hvectorfile = generateVectorFile(structure)[1]
    run(session, f"open {hvectorfile}")
    # clean up
    os.remove(hvectorfile)

def drawProtein(session, proteinfile: Path):
    '''
    Displays protein structure in ChimeraX session.
    '''
    run(session, f"open {proteinfile}")


def generateVectorFile(structure):
    '''
    Creates vector files for dipole and hydrophobic vector.
    '''
    com, dpv, hpv = createcoordinates(structure)
    cx, cy, cz = com
    dx, dy, dz = com + dpv
    hx, hy, hz = com + hpv
    dvectorfile =  "_dipole.bild"
    hvectorfile =  "_hydrophobe.bild"
    with open(dvectorfile, "w") as vector:
        body = f"""
                .color 1 1 1
                .sphere {cx} {cy} {cz} 2
                .color 0 0 1
                .arrow {cx} {cy} {cz} {dx} {dy} {dz} 1 3 0.75
                """
        vector.write(textwrap.dedent(body))
    with open(hvectorfile, "w") as vector:
        body = f"""
                .color 1 1 1
                .sphere {cx} {cy} {cz} 2
                .color 1 {140/255} 0
                .arrow {cx} {cy} {cz} {hx} {hy} {hz} 1 3 0.75
                """
        vector.write(textwrap.dedent(body))
    return dvectorfile, hvectorfile

def createcoordinates(structure):
    '''
    Creates coordinates for the center of mass, the dipole vector and the hydrophobic vector.
    '''
    c = structure.center_of_mass()
    d = structure.dipolevector()
    h = structure.hydrophobicvector() / 200
    return c, d, h

def generateStructure(proteinfile: Path):
    '''
    Creates protein structure from Builder from metrics module.
    '''
    parser = PDBParser(QUIET=True, structure_builder=Builder(is_AF=True))
    structure = parser.get_structure(proteinfile.stem, proteinfile)
    structure.measure()
    return structure

def main(proteinfile):
    '''
    Reads protein file to create structure. Visualizes protein structure, dipole vector and hydrophobic vector in Chimera X.
    '''
    proteinfile = Path(proteinfile)
    trimmedprotein = TRIM/f"{proteinfile.stem}_trim.pdb"
    s = generateStructure(proteinfile)
    s.save_pdb(trimmedprotein)
    drawProtein(session, trimmedprotein)
    drawDipole(session, s)
    drawHydrophobe(session, s)
    angles(session, s)
    drawArc(session, s)
    labels(session, proteinfile.stem)

parser = argparse.ArgumentParser()
parser.add_argument("proteinfile", help="A protein structure in pdb format.")
args = parser.parse_args()
main(args.proteinfile)
