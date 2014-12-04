#!/usr/bin/env python

class SeqPair(object):
    """Class to model DNA/mRNA sequence pairs"""

    def __init__(self, mseq, gseq, name):
        self.name = name
        self.mseq = mseq
        self.gseq = gseq
        self.gnuc_index = 0  #initialize counter for index in gen sequence
        self.mnuc_index = 0  #initialize counter for index in mRNA sequence
        self.codon_pos = 1 #keep track of position in codon
        self.aa_dict = {
            'F':['TTT','TTC'],
            'L':['TTA','TTG','CTT','CTC','CTA','CTG'],
            'I':['ATT','ATC','ATA'],
            'M':['ATG'],
            'V':['GTT','GTC','GTA','GTG'],
            'S':['TCT','TCC','TCA','TCG','AGT','AGC'],
            'P':['CCT','CCC','CCA','CCG'],
            'T':['ACT','ACC','ACA','ACG'],
            'A':['GCT','GCC','GCA','GCG'],
            'Y':['TAT','TAC'],
            'H':['CAT','CAC'],
            'Q':['CAA','CAG'],
            'N':['AAT','AAC'],
            'K':['AAA','AAG'],
            'D':['GAT','GAC'],
            'E':['GAA','GAG'],
            'C':['TGT','TGC'],
            'W':['TGG'],
            'R':['CGT','CGC','CGA','CGG','AGA','AGG'],
            'G':['GGT','GGC','GGA','GGG'],
            'STOP':['TAA','TAG','TGA'],
            }

    def index_nuc(self):
        return self.gnuc_index

    def index_mrna(self):
        return self.mnuc_index

    def index_position(self):
        return self.codon_pos

    def incr_nuc(self):
        self.gnuc_index += 1

    def incr_mrna(self):
        self.mnuc_index += 1

    def incr_pos(self):
        """maintains position within a codon as 1, 2, or 3"""
        if self.codon_pos < 3:
            self.codon_pos += 1
        else:
            self.codon_pos = 1

    def incr_all(self):
        """advances counters through DNA only"""
        self.incr_pos()
        self.incr_nuc()

    def lookup_gnuc(self):
        return self.gseq[self.gnuc_index]

    def lookup_mnuc(self):
        return self.mseq[self.mnuc_index]

    def lookup_gcodon(self):
        """returns DNA index as a codon"""
        if self.codon_pos == 1:  #uses codon position to determine slice
            return self.gseq[self.index_nuc():(self.index_nuc()+3)]
        elif self.codon_pos == 2:
            return self.gseq[(self.index_nuc()-1):(self.index_nuc()+2)]
        else:
            return self.gseq[(self.index_nuc()-2):self.index_nuc()+1]

    def lookup_mcodon(self):
        """returns mRNA index as a codon"""
        if self.codon_pos == 1:
            return self.mseq[self.index_mrna():(self.index_mrna()+3)]
        elif self.codon_pos == 2:
            return self.mseq[(self.index_mrna()-1):(self.index_mrna()+2)]
        else:
            return self.mseq[(self.index_mrna()-2):self.index_mrna()+1]

    def lookup_gaa(self):
        """returns aa specified by DNA codon"""
        codon = self.lookup_gcodon()
        #print codon
        for k1 in self.aa_dict.keys():
            for e in self.aa_dict.get(k1):
                if e == codon:
                    return k1

    def lookup_maa(self):
        """returns aa specified by mRNA codon"""
        codon = self.lookup_mcodon()
        #print codon
        for k1 in self.aa_dict.keys():
            for e in self.aa_dict.get(k1):
                if e == codon:
                    return k1


