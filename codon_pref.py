#!/usr/bin/env python

#import re
import sys

class Seq(object):
    """Class to model nucleotide/protein sequence pairs"""

    def __init__(self, seq1, seq2, name):
        self.name = name
        self.nuc_seq = seq1
        self.aa_seq = seq2
        self.aa_index = 0  #initialize counter for index in aa sequence
        self.nuc_index = 0  #separate counter for nuc sequence
        self.aa_dict = {
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
            }

    def aa_length(self):
        return len(self.aa_seq)

    def index_aa(self):
        return self.aa_index

    def incr_aa(self):
        self.aa_index += 1

    def index_nuc(self):
        return self.nuc_index

    def incr_nuc(self):
        self.nuc_index += 3

    def lookup_aa(self):
        return self.aa_seq[self.index_aa()]

    def lookup_nuc(self):
        """returns nucleotide index as a codon"""
        return self.nuc_seq[self.index_nuc():(self.index_nuc()+3)]

    def update_codons(self):
        """updates dict based on codon retrieved"""
        codon = self.lookup_nuc()
        #print codon
        for k1 in self.aa_dict.keys():
            for k2 in self.aa_dict[k1].keys():
                if k2 == codon:
                    self.aa_dict[k1][k2] += 1


def all_equal(l):
    """each time it is called, compares all amino acids in a given column of
    the alignment, returns True if all are equal"""
    prev = ''
    counter = 0  #need to know if object is first looked at
    for obj in l:
        cur = obj.lookup_aa()  #get the current amino acid
        if counter == 0:
            prev = cur
            counter += 1
        else:
            if cur == prev:
                pass
            else:
                return False  #only True if all amino acids are equal
    #print prev
    return True


nuc_file = sys.argv[1]
prot_file = sys.argv[2]
outfile = "howe_project_test_v3.txt"

with open(nuc_file, 'U') as f1, open(prot_file, 'U') as f2:
    nuc_dict = {}
    prot_dict = {}
    for curline in f1:
        if curline.startswith(">"):
            curline = curline.strip(">").strip("\n")
            ID = curline
            nuc_dict[ID] = ''
        else:
            curline = curline.strip("\n")
            nuc_dict[ID] += curline

    for curline in f2:
        if curline.startswith(">"):
            curline = curline.strip(">").strip("\n")
            ID = curline
            prot_dict[ID] = ''
        else:
            curline = curline.strip("\n")
            prot_dict[ID] += curline

    obj_list = []
    for a,b in zip(nuc_dict.keys(), prot_dict.keys()):
        if a == b:
            #creates an object only when same organism
            obj_list.append(Seq(nuc_dict.get(a), prot_dict.get(b), a))
        else:
            pass

    num = obj_list[0].aa_length()  #determine number of aa's in alignment
    #print num
    reps = 0
    while reps < num:  #need to go through each position
        if all_equal(obj_list):
            for obj in obj_list:
                obj.update_codons()
                obj.incr_aa()
                obj.incr_nuc()
        else:
            for obj in obj_list:
                if obj.lookup_aa() == '-':
                    obj.incr_aa()  #only increase nuc index if aa present
                else:
                    obj.incr_aa()
                    obj.incr_nuc()
        reps += 1

with open(outfile, 'w') as o:
    for obj in obj_list:
        o.write(obj.name + "\n")
        for k1 in obj.aa_dict:
            o.write(k1 + "\n")
            for k2 in obj.aa_dict[k1].keys():  #write aa
                o.write("\t" + k2 + ":" + str(obj.aa_dict[k1][k2]) + "\n")  #write codon usage
        o.write("\n")

"""
read in aligned nucleotide and protein seq files

parse each as per other programs

write a loop to initialize a new object for each sequence with the first full
word of the description line as the var name and the sequence as one of the
object's attributes

now the actual program part:
for each object, get the current amino acid for the index value
check each one against each other, if all of them match, get the associated
	codon used, and then add it to the codon usage variable
advance the index for the nucleotide sequence three

if the amino acids do not match, advance the index for the nucleotide sequence
	only if the current value is not '-'

no matter what though, advance the amino acid sequence index forward one

continue on in this manner until the end of each sequence is reached

write output to a file which has the sequence ID and the associated list of
	codons and the number of times they are used
"""
