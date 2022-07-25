import os
import sys
import textwrap
import argparse
from pathlib import Path
from Bio.PDB import PDBParser
from chimerax.core.commands import run
from chimerax.markers import MarkerSet
from chimerax.struct_measure.tool import StructMeasureTool
from chimerax.label import label_create

sys.path.append(Path(__file__).parent.parent.parent.as_posix())

from src.metrics import Builder

def labels(session, proteinfile):
    '''
    Create labels in ChimeraX session for the hydrophobic vector, dipole vector and the file name of the pdb structure.
    '''
    label_create(session, name = 'hv', text = 'hydrophobic vector', color = (255, 140, 0, 255), xpos = 0.5, ypos = 0.9, size = 21)
    label_create(session, name = 'dv', text = 'dipole vector', color = (0, 0, 255,255), xpos = 0.1, ypos = 0.9, size = 21)
    label_create(session, name = proteinfile.stem, text = proteinfile.stem, xpos = 0.45, ypos = 0.1, size = 21)

def angles(session, proteinfile: Path):
    '''
    Defines markers for the center of mass and the tips of the hydrophobic and dipole vector. Calculates the angle
    between the both vectors.
    Source for marker set: https://rbvi.github.io/chimerax-recipes/mark_blobs/mark_blobs.html
    '''
    c, d, h = createcoordinates(proteinfile)
    marker_set = MarkerSet(session, name = "markersforangle")
    marker_set.create_marker(d, (0, 0, 255,255), 0.5)
    marker_set.create_marker(c, (255, 255, 255,255), 0.5)
    marker_set.create_marker(h, (255, 140, 0, 255), 0.5)
    session.models.add([marker_set])
    run(session, f"select add #4")
    anglemaker = StructMeasureTool(session)
    anglemaker._create_angle()
    
def drawDipole(session, proteinfile: Path):
    '''
    Creates dipole vector in ChimeraX session.
    '''
    dvectorfile = generateVectorFile(proteinfile)[0]
    run(session, f"open {dvectorfile}")
    # clean up
    os.remove(dvectorfile)

def drawHydrophobe(session, proteinfile: Path):
    '''
    Creates hydrophobic vector in ChimeraX session.
    '''
    hvectorfile = generateVectorFile(proteinfile)[1]
    run(session, f"open {hvectorfile}")
    os.remove(hvectorfile)

def drawProtein(session, proteinfile: Path):
    '''
    Displays protein structure in ChimeraX session.
    '''
    run(session, f"open {proteinfile}")

def generateVectorFile(proteinfile: Path):
    '''
    Generates vector files for dipole and hydrophobic vectors.
    '''
    structure = generateStructure(proteinfile)
    coordinates = createcoordinates(proteinfile)
    cx, cy, cz = coordinates[0]
    dx, dy, dz = coordinates[1]
    hx, hy, hz = coordinates[2]/200
    dvectorfile = proteinfile.stem + f"{'dipole'}"+ ".bild"
    hvectorfile = proteinfile.stem + f"{'hydrophobe'}" + ".bild"
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

def createcoordinates(proteinfile: Path):
    '''
    Creates coordinates for the center of mass, the dipole vector and the hydrophobic vector.
    '''
    structure = generateStructure(proteinfile)
    structure.measure()
    c = structure.center_of_mass()
    d = structure.dipolevector()
    h = structure.hydrophobicvector()
    return c, d, h

def generateStructure(proteinfile: Path):
    '''
    Creates protein structure from Builder from metrics module.
    '''
    parser = PDBParser(QUIET=True, structure_builder=Builder())
    structure = parser.get_structure(proteinfile.stem, proteinfile)
    return structure
    
def main(proteinfile):
    '''
    Reads protein file to create structure. Visualizes protein structure, dipole vector and hydrophobic vector in Chimera X.
    '''
    proteinfile = Path(proteinfile)
    drawProtein(session,proteinfile)
    drawDipole(session, proteinfile)
    drawHydrophobe(session, proteinfile)
    angles(session, proteinfile)
    labels(session, proteinfile)

parser = argparse.ArgumentParser()
parser.add_argument("proteinfile", help="pass filename to script")
args = parser.parse_args()
main(args.proteinfile)