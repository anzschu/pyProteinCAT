import unittest
from pathlib import Path
from Bio.PDB import PDBParser

from src.metrics import Builder
from src.metrics.metrics import ModStructure, ModRes

class BuildingTests(unittest.TestCase):
    def setUp(self):
        self.pdb = Path("tests/pdb")
        self.parser = PDBParser(QUIET=1, structure_builder=Builder())

    def testBuilderWorks(self):
        for pdbfile in self.pdb.iterdir():
            s = self.parser.get_structure(pdbfile.stem, pdbfile)
            self.assertIsInstance(s, ModStructure)
            for residue in s.get_residues():
                self.assertIsInstance(residue, ModRes)
