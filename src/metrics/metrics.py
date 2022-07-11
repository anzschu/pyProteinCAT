import numpy as np
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa
from Bio.PDB.SASA import ShrakeRupley
from Bio.PDB.Structure import Structure
from Bio.PDB.StructureBuilder import StructureBuilder
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.PDB.Residue import Residue
from Bio.PDB.Polypeptide import three_to_one, standard_aa_names

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

hydrophobicityscale = { # from Guy, H.R., Biophys J. 47:61-70(1985)
    'A': 0.100, 
    'R':  1.910,
    'N':  0.480,
    'D':  0.780,
    'C': -1.420,
    'Q':  0.950,
    'E':  0.830,
    'G':  0.330,
    'H': -0.500,
    'I': -1.130,
    'L': -1.180,
    'K':  1.400,
    'M': -1.590,
    'F': -2.120,
    'P':  0.730,
    'S':  0.520,
    'T':  0.070,
    'W': -0.510,
    'Y': -0.210,
    'V': -1.270,
}

class ModStructure(Structure):
    
    def __init__(self, id): 
        Structure.__init__(self,id)
    
    def measure(self):
        '''
        Attach measurements as attributes to structure.
        '''
        self.sequence = self.sequencer()
        self.length = self.getlength()
        self.MWkDa = self.molecularweight()/1000
        self.fPos= self.aminocount()[0]/self.length
        self.fNeg = self.aminocount()[1]/self.length
        self.fFatty = self.aminocount()[2]/self.length
        self.sasa = self.calculate_sasa()
        self.nc = self.netcharge()
        self.ncd = self.netchargedensity()
        self.dpm = self.dipolemoment()
        self.guy = self.truehydrophobicity()
        self.hpm = self.hydrophobicmoment()
    
    def serializer(self):
        '''
        Return dictionary of all measurements from attributes.
        '''
        measurements = {
            'id': self.id,
            'length': self.length, 
            'mwkda': self.MWkDa,
            'fPos': self.fPos,
            'fNeg': self.fNeg,
            'fFatty': self.fFatty,
            'sasa': self.sasa,
            'nc': self.nc,
            'ncd': self.ncd,
            'dpm': self.dpm,
            'guy': self.guy,
            'hpm': self.hpm
        }
        return measurements

    def sequencer(self):
        '''
        Return polypeptide sequence from peptides.
        '''
        sequence = ""
        for residue in self.get_residues():
            sequence += residue.resletter
        return sequence
    
    def getlength(self):
        '''
        Return length of the protein sequence.
        '''
        return len(self.sequence)
        
    def molecularweight(self):
        '''
        Return molecular weight of the protein sequence.
        '''
        return ProteinAnalysis(str(self.sequence)).molecular_weight()
    
    def aminocount(self):
        '''
        Return the amount of positive, negative and hydrophobic residues.
        '''
        aminopos = self.sequence.count('K') + self.sequence.count('R')
        aminoneg = self.sequence.count('D') + self.sequence.count('E')
        aminohyd = sum(self.sequence.count(value) for value in 'FLIV')

        return aminopos, aminoneg, aminohyd

    def calculate_sasa(self):
        """
        Perform SASA calculation with ShrakeRupley algorithm for the individual
        chains (monomers). Attach results to residue objects as Residue.monosasa.
        """
        calculator = ShrakeRupley()
        for chain in self.get_chains():
            calculator.compute(chain, level="R") # calculate on the monomers
            for residue in chain.get_residues():
                residue.monosasa = residue.sasa.copy()
    
        return sum(residue.monosasa for residue in self.get_residues())

    def netcharge(self):
        '''
        Return the net charge of the protein from the amount of positive and negative amino acids.
        '''
        return self.aminocount()[0] - self.aminocount()[1]

    def netchargedensity(self):
        '''
        Return net charge density from net charge and total SASA of the
        residues.
        '''
        netchargedensity = self.netcharge()/self.sasa
        return netchargedensity

    def dipolemoment(self):
        '''
        Return dipole moment calculated from the dipole moments of the positive and negative residues.
        Source for dipolemoment equation: Felder, Prilusky, Silman, Sussman Nucleic Acids Research 2007
        '''
        dipolpos = np.zeros((1,3))
        dipolneg = np.zeros((1,3))
        dipolpos = sum(residue.center_of_mass() for residue in self.get_residues() if residue.resletter in 'KR')
        dipolneg = sum(residue.center_of_mass() for residue in self.get_residues() if residue.resletter in 'DE')
        dipolevector = dipolpos - dipolneg
        dipolemoment = 4.803 * np.linalg.norm(dipolevector)
        return dipolemoment

    def truehydrophobicity(self):
        '''
        Calculate total hydrophobicity for residues that are more than 25% exposed to the surface.
        '''
        truehydrophobicity = 0
        for residue in self.get_residues():
            #if residue.get_resname() in GLYXGLY_ASA:
            sasaratio =  residue.sasa / GLYXGLY_ASA[residue.resletter]
            if sasaratio > 0.25:
                truehydrophobicity += hydrophobicityscale[residue.resletter]
        return truehydrophobicity
    
    def hydrophobicmoment(self):
        '''
        Calculate first order hydrophobic moment. 
        Source for hydrophobicmoment equation: Silverman PNAS 2001, eq. 13
        '''
        hydrophobicmoment = 0
        for residue in self.get_residues():
            #if residue.get_resname() in GLYXGLY_ASA:
            hydrophobicmoment += (hydrophobicityscale[residue.resletter]* residue.sasa* (  residue.center_of_mass()- self.center_of_mass()))
        return np.linalg.norm(hydrophobicmoment)
    
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
        if not is_aa(resname, standard = True):
            return
        if field != " ":
            if field == "H":
                field = "H_" + resname

        residue_id = (field, resseq, icode)
        self.residue = ModRes(residue_id, resname, self.segid)  # ---> bespoke class
        self.chain.add(self.residue) 


if __name__ =='__main__':
    parser = PDBParser(QUIET=1, structure_builder=Builder())
    s = parser.get_structure("1ris", "data/data1/1ris.pdb")
    #s.calculate_sasa()
    s.measure()
    print(s.serializer())
