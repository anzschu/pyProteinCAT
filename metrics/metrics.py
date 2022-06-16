from ast import FormattedValue
from pydoc import apropos
import warnings
import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis
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

hydrophobicityscale = { # from Biophys J. 47:61-70(1985).
    'ALA': 0.100, 
    'ARG':  1.910,
    'ASN':  0.480,
    'ASP':  0.780,
    'CYS': -1.420,
    'GLN':  0.950,
    'GLU':  0.830,
    'GLY':  0.330,
    'HIS': -0.500,
    'ILE': -1.130,
    'LEU': -1.180,
    'LYS':  1.400,
    'MET': -1.590,
    'PHE': -2.120,
    'PRO':  0.730,
    'SER':  0.520,
    'THR':  0.070,
    'TRP': -0.510,
    'TYR': -0.210,
    'VAL': -1.270,
}

class ModStructure(Structure):
    
    def __init__(self, id): 
        Structure.__init__(self,id)
        self.sequence = self.sequencer()
        self.length = self.getlength()
        self.mw = self.molecularweight()
        self.aminopos = self.sequence.count('K') + self.sequence.count('R')
        self.aminoneg = self.sequence.count('D') + self.sequence.count('E')
        self.aminohyd = sum(self.sequence.count(value) for value in 'FLIV')

    def sequencer(self):
        ppb = PPBuilder()
        for pp in ppb.build_peptides(self):
            seq = pp.get_sequence()
        return seq        
    
    def getlength(self):
        '''
        Determine
        '''
        return len(self.sequence)
        
    def molecularweight(self):
        '''
        determine
        '''
        return ProteinAnalysis(str(self.sequence)).molecular_weight()
    
    def aminocount(self):
        '''
        Count
        '''
        aminocount = ProteinAnalysis(str(self.sequence)).count_amino_acids()
        aminopos = (aminocount['K'] + aminocount['R'])/self.length
        aminoneg = (aminocount['D'] + aminocount['E'])/self.length
        aminohyd = (aminocount['F'] + aminocount['L'] + aminocount['I'] + aminocount['V'])/self.length

        return aminopos, aminoneg, aminohyd

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

    def netcharge(self):
        '''
        Calculate
        '''
        return self.apos - self.aneg

    def netchargedensity(self):
        '''
        Calculate 
        '''
        netchargedensity = self.netcharge()/self.get_residues().monosasa
        return netchargedensity

    def dipolemoment(self):
        '''
        Calculate dipole moment
        '''
        dipolpos = np.zeros((1,3))
        dipolneg = np.zeros((1,3))
        dipolpos = sum(residue.center_of_mass for residue in self.get_residues() if residue.get_resname() in ['LYS', 'ARG'])
        dipolneg = sum(residue.center_of_mass for residue in self.get_residues() if residue.get_resname() in [ 'ASP', 'GLU'])
        dipolevector = dipolpos - dipolneg
        dipolemoment = 4.803 * np.linalg.norm(dipolevector)
        return dipolemoment

    def truehydrophobicity(self):
        '''
        Calculate
        '''
        truehydrophobicity = 0
        for residue in self.get_residues():
            if residue.get_resname() in GLYXGLY_ASA:
                sasaratio =  residue.sasa / GLYXGLY_ASA[residue.get_resname()]
                if sasaratio > 0.25:
                    truehydrophobicity += hydrophobicityscale[residue.get_resname()]
        return truehydrophobicity
    
    def hydrophobicmoment(self):
        '''
        Calculate first order hydrophobic moment source
        '''
        hydrophobicmoment = 0
        for residue in self.get_residues():
            if residue.get_resname() in GLYXGLY_ASA:
                hydrophobicmoment += (hydrophobicityscale[residue.get_resname()]* residue.sasa* (  residue.center_of_mass()- self.center_of_mass()))
    
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
    s = parser.get_structure("S6", "1ris.pdb")
    s.calculate_sasa()
    for res in s.get_residues():
        print(res, res.is_interface())
