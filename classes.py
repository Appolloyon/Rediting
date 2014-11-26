#!/usr/bin/env python

class Seq(object):
    """Class to model nucleotide/protein sequence pairs"""

    def __init__(self, mseq, gseq, name):
        self.name = name
        self.mseq = mseq
        self.gseq = gseq
        self.nuc_index = 0  #initialize counter for index in nuc sequence
        self.codon_index = 0  #separate counter for codons
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
            }

    def index_nuc(self):
        return self.nuc_index

    def index_codon(self):
        return self.codon_index

    def index_position(self):
        return self.codon_pos

    def incr_nuc(self):
        self.nuc_index += 1

    def incr_codon(self):
        if self.nuc_index % 3 == 0:
            self.codon_index += 1
        else:
            pass

    def incr_pos(self):
        if self.codon_pos < 3:
            self.codon_pos += 1
        else:
            self.codon_pos == 1

    def incr_all(self):
        self.incr_pos()
        self.incr_codon()
        self.incr_nuc()

    def lookup_gcodon(self):
        """returns nucleotide index as a codon"""
        if self.codon_pos == 1:
            return self.gseq[self.index_nuc():(self.index_nuc()+3)]
        elif self.codon_pos == 2:
            return self.gseq[(self.index_nuc()-1):(self.index_nuc()+1)]
        else:
            return self.gseq[(self.index_nuc()-2):self.index_nuc()]

    def lookup_mcodon(self):
        """returns nucleotide index as a codon"""
        if self.codon_pos == 1:
            return self.mseq[self.index_nuc():(self.index_nuc()+3)]
        elif self.codon_pos == 2:
            return self.mseq[(self.index_nuc()-1):(self.index_nuc()+1)]
        else:
            return self.mseq[(self.index_nuc()-2):self.index_nuc()]

    def lookup_gaa(self):
        """returns aa specified by codon"""
        codon = self.lookup_gcodon()
        #print codon
        for k1 in self.aa_dict.keys():
            for e in self.aa_dict.get(k1):
                if e == codon:
                    return k1

    def lookup_maa(self):
        """returns aa specified by codon"""
        codon = self.lookup_mcodon()
        #print codon
        for k1 in self.aa_dict.keys():
            for e in self.aa_dict.get(k1):
                if e == codon:
                    return k1


