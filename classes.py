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
        #print self.gnuc_index

    def incr_mrna(self):
        self.mnuc_index += 1
        #print self.mnuc_index

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


class CodonPair(object):
    """Class to model DNA/mRNA sequence pairs"""

    def __init__(self, mseq, gseq, name):
        self.name = name
        self.mseq = mseq
        self.gseq = gseq
        self.gnuc_index = 0  #initialize counter for index in gen sequence
        self.mnuc_index = 0  #initialize counter for index in mRNA sequence
        self.codon_pos = 1 #keep track of position in codon
        self.gnuc_aa_dict = {
            'F':{'TTT':0,'TTC':0},
            'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
            'I':{'ATT':0,'ATC':0,'ATA':0},
            'M':{'ATG':0},
            'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
            'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
            'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
            'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
            'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
            'Y':{'TAT':0,'TAC':0},
            'H':{'CAT':0,'CAC':0},
            'Q':{'CAA':0,'CAG':0},
            'N':{'AAT':0,'AAC':0},
            'K':{'AAA':0,'AAG':0},
            'D':{'GAT':0,'GAC':0},
            'E':{'GAA':0,'GAG':0},
            'C':{'TGT':0,'TGC':0},
            'W':{'TGG':0},
            'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
            'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0},
            'STOP':{'TAA':0,'TAG':0,'TGA':0}
            }
        self.mnuc_aa_dict = {
            'F':{'TTT':0,'TTC':0},
            'L':{'TTA':0,'TTG':0,'CTT':0,'CTC':0,'CTA':0,'CTG':0},
            'I':{'ATT':0,'ATC':0,'ATA':0},
            'M':{'ATG':0},
            'V':{'GTT':0,'GTC':0,'GTA':0,'GTG':0},
            'S':{'TCT':0,'TCC':0,'TCA':0,'TCG':0,'AGT':0,'AGC':0},
            'P':{'CCT':0,'CCC':0,'CCA':0,'CCG':0},
            'T':{'ACT':0,'ACC':0,'ACA':0,'ACG':0},
            'A':{'GCT':0,'GCC':0,'GCA':0,'GCG':0},
            'Y':{'TAT':0,'TAC':0},
            'H':{'CAT':0,'CAC':0},
            'Q':{'CAA':0,'CAG':0},
            'N':{'AAT':0,'AAC':0},
            'K':{'AAA':0,'AAG':0},
            'D':{'GAT':0,'GAC':0},
            'E':{'GAA':0,'GAG':0},
            'C':{'TGT':0,'TGC':0},
            'W':{'TGG':0},
            'R':{'CGT':0,'CGC':0,'CGA':0,'CGG':0,'AGA':0,'AGG':0},
            'G':{'GGT':0,'GGC':0,'GGA':0,'GGG':0},
            'STOP':{'TAA':0,'TAG':0,'TGA':0}
            }


    def index_nuc(self):
        return self.gnuc_index

    def index_mrna(self):
        return self.mnuc_index

    def index_position(self):
        return self.codon_pos

    def incr_nuc(self):
        self.gnuc_index += 1
        #print self.gnuc_index

    def incr_mrna(self):
        self.mnuc_index += 1
        #print self.mnuc_index

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

    def update_gcodons(self):
        """updates dict based on codon retrieved"""
        gcodon = self.lookup_gcodon()
        #print codon
        for k1 in self.gnuc_aa_dict.keys():
            for k2 in self.gnuc_aa_dict[k1].keys():
                if k2 == gcodon:
                    self.gnuc_aa_dict[k1][k2] += 1

    def update_mcodons(self):
        """updates dict based on codon retrieved"""
        mcodon = self.lookup_mcodon()
        #print codon
        for k1 in self.mnuc_aa_dict.keys():
            for k2 in self.mnuc_aa_dict[k1].keys():
                if k2 == mcodon:
                    self.mnuc_aa_dict[k1][k2] += 1

