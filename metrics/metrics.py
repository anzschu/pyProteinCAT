import warnings

from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.PDB.Residue import Residue, DisorderedResidue
from Bio.PDB.Polypeptide import three_to_one, standard_aa_names
from Bio.PDB.PDBExceptions import PDBConstructionWarning, PDBConstructionException

GLYXGLY_ASA = { # from Miller, Janin et al (1987)
    "A": 113,
    "R": 241,
    "N": 158,
    "D": 151,
    "C": 140,
    "Q": 189,
    "E": 183,
    "G":  85,
    "H": 194,
    "I": 182,
    "L": 180,
    "K": 211,
    "M": 204,
    "F": 218,
    "P": 143,
    "S": 122,
    "T": 146,
    "W": 259,
    "Y": 229,
    "V": 160,
}

class ModStructure(Structure):
    def calculate_sasa(self):
        """
        Perform SASA calculation with ShrakeRupley algorithm for the individual
        chains (monomers) and the whole structure (complex). Attach results to
        residue objects as Residue.monosasa and Residue.complexsasa.
        """
        calculator = ShrakeRupley()
        for chain in self.get_chains():
            calculator.compute(chain, level="R") # calculate on the monomers
            for res in chain.get_residues():
                res.monosasa = res.sasa.copy()

        calculator.compute(self, level="R") # calculate on the complex
        for res in self.get_residues():
            res.complexsasa = res.sasa.copy()

    def __str__(self):
        return f"ModStructure instance {self.id}"

    def __repr__(self):
        return self.__str__()

class ModRes(Residue):
    def __init__(self, id, resname, segid):
        Residue.__init__(self, id, resname, segid)
        self.resletter = three_to_one(resname) if resname in standard_aa_names else ""

    def is_buried(self, threshold=0.25):
        """
        Return boolean depending on ratio of residue SASA with respect to
        literature values in GLYXGLY_ASA and a threshold (default 0.25).
        """
        try:
            return self.monosasa / GLYXGLY_ASA[self.resletter] <= threshold
        except AttributeError:
            raise AttributeError("Attrib. monosasa not set: did you call structure.calculate_sasa()")

    def is_interface(self, thresh_buried=0.25, thresh_inter=0.05):
        """
        Return boolean depending on ratio of residue SASA change with respect to
        literature values in GLYXGLY_ASA and a threshold (default 0.05).
        """
        if self.is_buried(threshold = thresh_buried):
            return False
        else:
            return (self.monosasa - self.complexsasa) / GLYXGLY_ASA[self.resletter] > thresh_inter

    def __str__(self):
        return f"ModRes instance {self.resname}"

    def __repr__(self):
        return self.__str__()


class Builder(StructureBuilder):
    """
    Bespoke subclassing of Bio.PDB.StructureBuilder.StructureBuilder to
    incorporate ModRes and ModStructure into the Bio.PDB hierarchy.

    Usage:
    Pass a Builder instance to the PDBParser as a structure_builder parameter.

    Example:
    >>> parser = PDBParser(QUIET=1, structure_builder=Builder())
    >>> s = parser.get_structure("ab1", "5iy5_AB1.pdb")
    """
    def init_structure(self, structure_id):
        self.structure = ModStructure(structure_id)         # ---> bespoke class

    def init_residue(self, resname, field, resseq, icode):
        if field != " ":
            if field == "H":
                field = "H_" + resname

        res_id = (field, resseq, icode)
        self.residue = ModRes(res_id, resname, self.segid)  # ---> bespoke class
        self.chain.add(self.residue)


if __name__ =='__main__':
    parser = PDBParser(QUIET=1, structure_builder=Builder())
    s = parser.get_structure("ab1", "5iy5_AB1.pdb")
    s.calculate_sasa()
    for res in s.get_residues():
        print(res, res.is_interface())
