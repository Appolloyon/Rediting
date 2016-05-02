#!/usr/bin/env python

import sys
import subprocess
import scipy.stats as st
import matplotlib.pyplot as plt

from util import sequence

infiles = sys.argv[3:]
outfile = sys.argv[2]
num_gens = int(sys.argv[1])


# Global variable for possible bases in generated sequences
bases = 'AGTC'

# This bit will get the relevant information
gene_dict = {}

for infile in infiles:
    name = infile.split('.')[0]
    subprocess.call(["/Users/cklinger/git/Rediting/rediting/detect_edit_clusters.py", "-in",
        infile, "-out", "test.csv", "-n", name, "-r", "rna", "-gen", "gen", "--simulation"])
    with open('tempfile.csv','U') as i:
        for line in i:
            llist = line.strip('\n').split(',')
        gene_dict[llist[0]] = []
        gene_dict[llist[0]].extend([llist[1],llist[2],llist[3]])
        edit_list = []
        for val in llist[4:len(llist)]:
            edit_list.append(float(val))
        gene_dict[llist[0]].append(edit_list)

with open(outfile,'w') as o:
    for k,v in sorted(gene_dict.items()):
        p_values = []
        length = int(v[0])
        num_edits = int(v[1])
        exp_dist = v[3]
        for i in range(num_gens):
            seq1 = sequence.generate_start_sequence(length,bases)
            # Add mutations to simulate edits
            seq2 = sequence.mutate_sequence(seq1,num_edits,bases)
            sim_dist = sequence.calc_cluster_score(seq1,seq2)
            # Check for appropriate t test
            p_val = st.levene(exp_dist,sim_dist)
            #if rmath.equal_vars(exp_dist,sim_dist):
                #p_val = st.ttest_ind(exp_dist,sim_dist)
            #else:
                #p_val = st.ttest_ind(exp_dist,sim_dist,equal_var=False)
            # 1st value is the test value itself, second is the p-value
            p_values.append(p_val[1])

            #n, bins, patches = plt.hist(sim_dist,50,facecolor='green')
            #n, bins, patches = plt.hist(exp_dist,50,facecolor='blue')
            #plt.show()

        sig_pvals = 0
        for p_val in p_values:
            if p_val < 0.05:
                sig_pvals += 1
        percent_sig = (float(sig_pvals)/len(p_values))
        o.write("%s,%s,%s,%.2f" % (k,length,num_edits,percent_sig))
        o.write("\n")

        #print k
        #n, bins, patches = plt.hist(p_values,50)
        #plt.xlabel('p value')
        #plt.ylabel('frequency')
        #plt.show()


