drawVectors.py
===========

**drawVectors** represents dipole and hydrophobicity vectors over a protein's structure. Opens
a ChimeraX GUI.

Usage: 
    $ chimera --script "drawVectors.py 1abc.pdb"

.. py:function:: drawVectors.main(proteinfile)

    Reads protein file to create structure.
    Visualizes protein structure, dipole vector and hydrophobic vector in
    Chimera X.

    :param proteinfile: PDB File
    :type proteinfile: str

.. py:function:: drawVectors.generateStructure(proteinfile)

    Creates protein structure from Builder from metrics module.

    :param proteinfile: File path to PDB file.
    :type proteinfile: Path
    :return: Protein structure

.. py:function:: drawVectors.createcoordinates(proteinfile)

    Creates coordinates for the center of mass, the dipole vector and th
    hydrophobic vector.

    :param proteinfile: File path to PDB file.
    :type proteinfile: Path
    :return: Coordinates of the center of mass, dipole vector and hydrophobic vector
    :rtype: tuple

.. py:function:: drawVectors.generateVectorFile(proteinfile)

    Creates vector files for dipole and hydrophobic vector.

    :param proteinfile: File path to PDB file.
    :type proteinfile: Path
    :return: Vectorfiles for dipole and hydrophobic vector.
    :rtype: str

.. py:function:: drawVectors.drawProtein(session, proteinfile)

    Displays protein structure in ChimeraX session.

    :param session: Global variable pointing to program instance, known to Chimera environment.
    :param proteinfile: File path to PDB file.

.. py:function:: drawVectors.drawVectors(session, *vectorfiles)

    Consume as many files as fed.

    :param session: Global variable pointing to program instance, known to Chimera environment.
    :param vectorfile: File to create vector in Chimera X.

.. py:function:: drawVectors.angles(session, proteinfile)

    Defines markers for the center of mass and the tip of the hydrophobic and
    dipole vector. Calculates the angle between the both vectors. Based on:
    https://rbvi.github.io/chimerax-recipes/mark_blobs/mark_blobs.html

    :param session: Global variable pointing to program instance, known to Chimera environment.
    :param proteinfile: File path to PDB file.
    :type proteinfile: Path

.. py:function:: drawVectors.drawArc(session, proteinfile)
    
    Draws arc between both vectors.

    :param session: Global variable pointing to program instance, known to Chimera environment.
    :param proteinfile: File path to PDB file.
    :type proteinfile: Path

.. py:function:: drawVectors.labels(session, proteinfile)

    Creates labels for hydrophobic vector, dipole vector and protein name.

    :param session: Global variable pointing to program instance, known to Chimera environment.
    :param proteinfile: File path to PDB file.
    :type proteinfile: Path