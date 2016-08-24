#!/usr/bin/env python

import copy

class SeqPair(object):
    """Class to model DNA/mRNA sequence pairs"""

    def __init__(self, mseq, gseq, name):
        self.name = name
        self.mseq = mseq
        self.gseq = gseq
        # Next two lines initialize index counters
        self.gnuc_index = 0
        self.mnuc_index = 0
        # Could keep track of codon position in each but we take
        # the genomic sequence as correct, in case of indels, etc.
        self.codon_pos = 1
        # These dicts need to be separate as they are populated separately
        # This could in theory be done more succintly with a deep copy
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
        self.mnuc_aa_dict = copy.deepcopy(self.gnuc_aa_dict)
        # Unlike the last two dict which keep track of codons
        # This last dict is to track all possible base changes
        self.transition_dict = {
            'a_t':0, 'a_g':0, 'a_c':0,
            't_a':0, 't_g':0, 't_c':0,
            'g_a':0, 'g_t':0, 'g_c':0,
            'c_a':0, 'c_t':0, 'c_g':0
            }
        self.first_pos_transdict = copy.deepcopy(self.transition_dict)
        self.second_pos_transdict = copy.deepcopy(self.transition_dict)
        self.third_pos_transdict = copy.deepcopy(self.transition_dict)


    def index_nuc(self):
        """simple access method"""
        return self.gnuc_index

    def index_mrna(self):
        """simple access method"""
        return self.mnuc_index

    def index_position(self):
        """simple access method"""
        return self.codon_pos

    def incr_nuc(self):
        """increments DNA index only"""
        self.gnuc_index += 1

    def incr_mrna(self):
        """increments mRNA index only"""
        self.mnuc_index += 1

    def incr_pos(self):
        """maintains position within a codon as 1, 2, or 3"""
        if self.codon_pos < 3:
            # If less than three we simply increment
            self.codon_pos += 1
        else:
            # If the codon position hits three than we
            # Cannot increment and need to reset it to 1
            self.codon_pos = 1

    def incr_all(self):
        """advances counters through DNA only"""
        self.incr_pos()
        self.incr_nuc()

    def lookup_gnuc(self):
        """looks up current DNA nucleotide"""
        return self.gseq[self.gnuc_index]

    def lookup_mnuc(self):
        """looks up current mRNA nucleotide"""
        return self.mseq[self.mnuc_index]

    def lookup_gcodon(self):
        """returns DNA index as a codon"""
        # Uses codon position to determine slice
        if self.codon_pos == 1:
            # Gets current and next two
            return self.gseq[self.index_nuc():(self.index_nuc()+3)]
        elif self.codon_pos == 2:
            # Gets one behind and one in front
            return self.gseq[(self.index_nuc()-1):(self.index_nuc()+2)]
        else:
            # Gets two behind and current
            return self.gseq[(self.index_nuc()-2):self.index_nuc()+1]

    def lookup_mcodon(self):
        """returns mRNA index as a codon"""
        # Exact same mechanic as for gcodons
        if self.codon_pos == 1:
            return self.mseq[self.index_mrna():(self.index_mrna()+3)]
        elif self.codon_pos == 2:
            return self.mseq[(self.index_mrna()-1):(self.index_mrna()+2)]
        else:
            return self.mseq[(self.index_mrna()-2):self.index_mrna()+1]

    def lookup_gaa(self):
        """returns aa specified by DNA codon"""
        # Uses aa dict but does not populate it
        codon = self.lookup_gcodon()
        for k1 in self.gnuc_aa_dict.keys():
            for e in self.gnuc_aa_dict.get(k1):
                if e == codon:
                    return k1
        # Contingency in case of non-canonical bases
        return 'X'

    def lookup_maa(self):
        """returns aa specified by mRNA codon"""
        # Uses aa dict but does not populate it
        codon = self.lookup_mcodon()
        for k1 in self.mnuc_aa_dict.keys():
            for e in self.mnuc_aa_dict.get(k1):
                if e == codon:
                    return k1
        # Contingency in case of non-canonical bases
        return 'X'

    def update_gcodons(self):
        """updates DNA dict based on codon retrieved"""
        # Actually updates the dict
        # If the codon is not recognized it will not be added
        gcodon = self.lookup_gcodon()
        for k1 in self.gnuc_aa_dict.keys():
            for k2 in self.gnuc_aa_dict[k1].keys():
                if k2 == gcodon:
                    self.gnuc_aa_dict[k1][k2] += 1

    def update_mcodons(self):
        """updates mRNA dict based on codon retrieved"""
        # Same mechanic
        mcodon = self.lookup_mcodon()
        for k1 in self.mnuc_aa_dict.keys():
            for k2 in self.mnuc_aa_dict[k1].keys():
                if k2 == mcodon:
                    self.mnuc_aa_dict[k1][k2] += 1

    def update_base_pref(self,tdict):
        """general update function"""
        gnuc = self.lookup_gnuc()
        mnuc = self.lookup_mnuc()

        if gnuc == "A" and mnuc == "T":
            tdict['a_t'] += 1
        elif gnuc == "A" and mnuc == "G":
            tdict['a_g'] += 1
        elif gnuc == "A" and mnuc == "C":
            tdict['a_c'] += 1
        elif gnuc == "T" and mnuc == "A":
            tdict['t_a'] += 1
        elif gnuc == "T" and mnuc == "G":
            tdict['t_g'] += 1
        elif gnuc == "T" and mnuc == "C":
            tdict['t_c'] += 1
        elif gnuc == "G" and mnuc == "A":
            tdict['g_a'] += 1
        elif gnuc == "G" and mnuc == "T":
            tdict['g_t'] += 1
        elif gnuc == "G" and mnuc == "C":
            tdict['g_c'] += 1
        elif gnuc == "C" and mnuc == "A":
            tdict['c_a'] += 1
        elif gnuc == "C" and mnuc == "T":
            tdict['c_t'] += 1
        elif gnuc == "C" and mnuc == "G":
            tdict['c_g'] += 1
        else:
            # Won't record information regarding
            # non-canonical/uncertain bases
            pass

    def update_transdict(self):
        """updates transition_dict based on observed bases"""
        self.update_base_pref(self.transition_dict)

    def update_first_pos_transdict(self):
        """updates dict only if first codon position"""
        self.update_base_pref(self.first_pos_transdict)

    def update_second_pos_transdict(self):
        """updates dict only if second codon position"""
        self.update_base_pref(self.second_pos_transdict)

    def update_third_pos_transdict(self):
        """updates dict only if third codon position"""
        self.update_base_pref(self.third_pos_transdict)


class RefPair(object):
    """Leighter weight class to model DNA/ref sequence pairs.
    This is similar to the SeqPair class but has slightly different
    mechanics and is less computationally intensive"""

    def __init__(self, rseq, gseq, name):
        self.name = name
        self.rseq = rseq
        self.gseq = gseq
        # Next two lines initialize counters
        self.rnuc_index = 0
        self.gnuc_index = 0
        # Next two lines initialize codon positions
        self.rcodon_pos = 1
        self.gcodon_pos = 1

    def index_nuc(self):
        """simple access method"""
        return self.gnuc_index

    def index_ref(self):
        """simple access method"""
        return self.rnuc_index

    def index_rposition(self):
        """simple access method"""
        return self.rcodon_pos

    def index_gposition(self):
        """simple access method"""
        return self.gcodon_pos

    def incr_ref(self):
        """increments ref index only"""
        self.rnuc_index += 1

    def incr_gen(self):
        """increments gen index only"""
        self.gnuc_index += 1

    def incr_rpos(self):
        """maintains position within a ref codon as 1, 2, or 3"""
        if self.rcodon_pos < 3:
            self.rcodon_pos += 1
        else:
            self.rcodon_pos = 1

    def incr_gpos(self):
        """maintatins position within gen codon as 1, 2, or 3"""
        if self.gcodon_pos < 3:
            self.gcodon_pos += 1
        else:
            self.gcodon_pos = 1

    def incr_all_ref(self):
        """advances counters through ref only"""
        self.incr_rpos()
        self.incr_ref()

    def incr_all_gen(self):
        """advances counters through gen only"""
        self.incr_gpos()
        self.incr_gen()


class Simulation(object):
    """Class to track cumulative information from multiple runs of the
    weighted sequence mutation program DNA/mRNA sequence pairs"""

    def __init__(self):
        self.transition_dict = {
            'a_t':0, 'a_g':0, 'a_c':0,
            't_a':0, 't_g':0, 't_c':0,
            'g_a':0, 'g_t':0, 'g_c':0,
            'c_a':0, 'c_t':0, 'c_g':0
            }
        self.codon_pos_dict = {
                1:0, 2:0, 3:0
            }

    def update_transdict(self,start_base,end_base):
        """updates transition_dict based on observed bases"""
        if start_base == "A" and end_base == "T":
            self.transition_dict['a_t'] += 1
        elif start_base == "A" and end_base == "G":
            self.transition_dict['a_g'] += 1
        elif start_base == "A" and end_base == "C":
            self.transition_dict['a_c'] += 1
        elif start_base == "T" and end_base == "A":
            self.transition_dict['t_a'] += 1
        elif start_base == "T" and end_base == "G":
            self.transition_dict['t_g'] += 1
        elif start_base == "T" and end_base == "C":
            self.transition_dict['t_c'] += 1
        elif start_base == "G" and end_base == "A":
            self.transition_dict['g_a'] += 1
        elif start_base == "G" and end_base == "T":
            self.transition_dict['g_t'] += 1
        elif start_base == "G" and end_base == "C":
            self.transition_dict['g_c'] += 1
        elif start_base == "C" and end_base == "A":
            self.transition_dict['c_a'] += 1
        elif start_base == "C" and end_base == "T":
            self.transition_dict['c_t'] += 1
        elif start_base == "C" and end_base == "G":
            self.transition_dict['c_g'] += 1
        else:
            # Won't record information regarding
            # non-canonical/uncertain bases
            pass

    def update_codon_dict(self,codon_pos):
        """Updates codon_pos_dict based on observed position"""
        self.codon_pos_dict[codon_pos] += 1

    def get_value(self,tdict):
        """Finds the total for all values in a dictionary"""
        total = 0
        for val in tdict.values():
            total += val
        return total

    def get_codon_percent(self):
        """Prints information about editing codon positions"""
        total = self.get_value(self.codon_pos_dict)
        for k,v in self.codon_pos_dict.items():
            print "%s pos %.2f" % (k,(100 * v/float(total)))

    def get_base_conversion(self):
        """Prints information about base conversions"""
        total = self.get_value(self.transition_dict)
        for k,v in self.transition_dict.items():
            print "%s num %.2f" % (k,(100 * v/float(total)))

    def write_sim_information(self,name,file_obj):
        """Writes information to a file"""
        codon_total = float(self.get_value(self.codon_pos_dict))
        base_total = float(self.get_value(self.transition_dict))
        file_obj.write("%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,"
            "%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f" % (name,
            (100 * self.codon_pos_dict.get(1)/codon_total),
            (100 * self.codon_pos_dict.get(2)/codon_total),
            (100 * self.codon_pos_dict.get(3)/codon_total),
            (100 * self.transition_dict.get('a_t')/base_total),
            (100 * self.transition_dict.get('a_g')/base_total),
            (100 * self.transition_dict.get('a_c')/base_total),
            (100 * self.transition_dict.get('t_a')/base_total),
            (100 * self.transition_dict.get('t_g')/base_total),
            (100 * self.transition_dict.get('t_c')/base_total),
            (100 * self.transition_dict.get('g_a')/base_total),
            (100 * self.transition_dict.get('g_t')/base_total),
            (100 * self.transition_dict.get('g_c')/base_total),
            (100 * self.transition_dict.get('c_a')/base_total),
            (100 * self.transition_dict.get('c_t')/base_total),
            (100 * self.transition_dict.get('c_g')/base_total)))
        file_obj.write("\n")
