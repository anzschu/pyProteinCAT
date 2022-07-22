from src.metrics import Builder
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

sys.path.append(Path(__file__).parent.parent.parent.as_posix())


def labels(session, proteinfile):
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
        name=proteinfile.stem,
        text=proteinfile.stem,
        xpos=0.45,
        ypos=0.1,
        size=21
    )


def angles(session, proteinfile: Path):
    '''
    Defines markers for the center of mass and the tip of the hydrophobic and
    dipole vector. Calculates the angle between the both vectors and displays
    this. Source for marker set:
    https://rbvi.github.io/chimerax-recipes/mark_blobs/mark_blobs.html
    '''
    c, d, h = createcoordinates(proteinfile)
    marker_set = MarkerSet(session, name="markersforangle")
    marker5 = marker_set.create_marker(c + d, (0, 0, 255, 255), 0.5)
    marker6 = marker_set.create_marker(c, (255, 255, 255, 255), 0.5)
    marker7 = marker_set.create_marker(c + h, (255, 140, 0, 255), 0.5)
    session.models.add([marker_set])
    run(session, f"select add #4")
    anglemaker = StructMeasureTool(session)
    anglemaker._create_angle()


def drawDipole(session, proteinfile: Path):
    dvectorfile = generateVectorFile(proteinfile)[0]
    run(session, f"open {dvectorfile}")
    # clean up
    os.remove(dvectorfile)


def drawHydrophobe(session, proteinfile: Path):
    hvectorfile = generateVectorFile(proteinfile)[1]
    run(session, f"open {hvectorfile}")
    # clean up
    os.remove(hvectorfile)


def drawProtein(session, proteinfile: Path):
    run(session, f"open {proteinfile}")


def generateVectorFile(proteinfile: Path):
    com, dpv, hpv = createcoordinates(proteinfile)
    cx, cy, cz = com
    dx, dy, dz = com + dpv
    hx, hy, hz = com + hpv
    dvectorfile = proteinfile.stem + f"{'dipole'}" + ".bild"
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
    structure = generateStructure(proteinfile)
    structure.measure()
    c = structure.center_of_mass()
    d = structure.dipolevector()
    h = structure.hydrophobicvector() / 200
    return c, d, h


def generateStructure(proteinfile: Path):
    parser = PDBParser(QUIET=True, structure_builder=Builder())
    structure = parser.get_structure(proteinfile.stem, proteinfile)
    structure.measure()
    return structure


def drawArc(session, proteinfile):
    com, dpv, hpv = createcoordinates(proteinfile)
    dx, dy, dz = com + dpv / np.linalg.norm(dpv) * 10
    hx, hy, hz = com + hpv / np.linalg.norm(hpv) * 10
    body = f"""
            .color 0 0.5 0.5
            .transparency 1
            .dot {dx} {dy} {dz}
            .transparency 0
            .draw {hz} {hy} {hz}
            """
    with open(f'{proteinfile.stem}_arc.bild', "w") as arcfile:
        arcfile.write(textwrap.dedent(body))

    run(session, f"open {proteinfile.stem}_arc.bild")
    # clean up
    os.remove(f'{proteinfile.stem}_arc.bild')


def main(proteinfile):
    proteinfile = Path(proteinfile)
    drawProtein(session, proteinfile)
    drawDipole(session, proteinfile)
    drawHydrophobe(session, proteinfile)
    angles(session, proteinfile)
    drawArc(session, proteinfile)
    labels(session, proteinfile)


parser = argparse.ArgumentParser()
parser.add_argument("proteinfile", help="pass filename to script")
args = parser.parse_args()
main(args.proteinfile)
