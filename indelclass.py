#!/usr/bin/env python

class IndelSeq(object):
    """Class to keep track of index in a reference sequence"""

    def __init__(self, rseq, name):
        self.name = name
        self.rseq = rseq
        #self.qseq = qseq
        self.rseq_index = 0  #initialize counter for index in ref sequence
        #self.qseq_index = 0  #initialize counter for index in query sequence
        #self.codon_pos = 1 #keep track of position in codon
        #self.insert_dict = {}
        #self.del_dict = {}

    def index_ref(self):
        return self.rseq_index

    #def index_query(self):
        #return self.qseq_index

    #def index_position(self):
        #return self.codon_pos

    def incr_ref(self):
        self.rseq_index += 1
        #print self.gnuc_index

    #def incr_query(self):
        #self.qseq_index += 1
        #print self.mnuc_index

    #def incr_pos(self):
        #"""maintains position within a codon as 1, 2, or 3"""
        #if self.codon_pos < 3:
            #self.codon_pos += 1
        #else:
            #self.codon_pos = 1

    #def incr_all(self):
        #"""advances counters through both sequences"""
        #self.incr_ref()
        #self.incr_query()

    def lookup_rres(self):
        return self.rseq[self.rseq_index]

    def lookup_qres(self):
        return self.qseq[self.qseq_index]

    #def lookup_gcodon(self):
        #"""returns DNA index as a codon"""
        #if self.codon_pos == 1:  #uses codon position to determine slice
            #return self.gseq[self.index_nuc():(self.index_nuc()+3)]
        #elif self.codon_pos == 2:
            #return self.gseq[(self.index_nuc()-1):(self.index_nuc()+2)]
        #else:
            #return self.gseq[(self.index_nuc()-2):self.index_nuc()+1]

    #def lookup_mcodon(self):
        #"""returns mRNA index as a codon"""
        #if self.codon_pos == 1:
            #return self.mseq[self.index_mrna():(self.index_mrna()+3)]
        #elif self.codon_pos == 2:
            #return self.mseq[(self.index_mrna()-1):(self.index_mrna()+2)]
        #else:
            #return self.mseq[(self.index_mrna()-2):self.index_mrna()+1]

    #def lookup_gaa(self):
        #"""returns aa specified by DNA codon"""
        #codon = self.lookup_gcodon()
        #print codon
        #for k1 in self.aa_dict.keys():
            #for e in self.aa_dict.get(k1):
                #if e == codon:
                    #return k1

    #def lookup_maa(self):
        #"""returns aa specified by mRNA codon"""
        #codon = self.lookup_mcodon()
        #print codon
        #for k1 in self.aa_dict.keys():
            #for e in self.aa_dict.get(k1):
                #if e == codon:
                    #return k1


