class Residue:
    def __init__(self, name, pos, chain=""):
        self.pos = pos
        self.name = name
        self.chain = chain

    def addChain(self, chain):
        self.chain = chain


    def is_upstream(self, other: 'Residue'):
        if self < other:
            return True
        elif self == other:
            raise RuntimeError("Different residues in same positon?")
        else:
            return False

    def __lt__(self, other):
        return self.pos < other.pos

    def __eq__(self, other):
        return self.pos == other.pos

    def __str__(self):
        return f"{self.name} @ {self.pos} in {self.chain}"

    def __repr__(self):
        return self.__str__()





if __name__ == '__main__':
    sequence = []
    with open("sequence.csv", "r") as csvfile:
        for line in csvfile.readlines():
            arguments = line.split(",")
            r = Residue(*arguments)
            sequence.append(r)

    r0, r1, _ = sequence

    print(r0.is_upstream(r1))
