from Bio import Phylo, SeqIO, AlignIO
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from collections import defaultdict
import numpy as np
import json
from Bio.Seq import translate

if __name__ == '__main__':
    import argparse
    parser = parser = argparse.ArgumentParser(description='assemble a script for pymol',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--reference', type=str, help="reference genome of EVD68")
    parser.add_argument('--tree', type=str, help="input tree (refined)")
    parser.add_argument('--mutations', help="mutation JSON file")
    parser.add_argument('--output', help="name of pymol file")
    args = parser.parse_args()

    ref_seq = SeqIO.read(args.reference, 'genbank')
    T = Phylo.read(args.tree, 'newick')
    proteins = ref_seq.features[3:-1]

    with open(args.mutations, 'r') as fh:
        tmp = json.load(fh)
        muts = tmp["nodes"]

    codon_changes_syn = defaultdict(list)
    codon_changes_nonsyn = defaultdict(list)
    exclude_terminal = False

    for n in T.find_clades():
        for protein in proteins:
            pname = protein.qualifiers['gene'][0]
            for pos in range(protein.location.start, protein.location.end,3):
                codon_pos = (pos-protein.location.start)//3
                anc_codon = muts[n.name]['sequence'][pos:pos+3]
                try:
                    anc_aa = translate(anc_codon)
                except:
                    continue
                for c in n:
                    if exclude_terminal and c.is_terminal():
                        continue
                    der_codon = muts[c.name]['sequence'][pos:pos+3]
                    if 'N' in der_codon or '-' in der_codon:
                        continue
                    if anc_codon!=der_codon:
                        try:
                            der_aa = translate(der_codon)
                        except:
                            continue
                        if anc_aa==der_aa:
                            codon_changes_syn[(pname, codon_pos)].append((anc_codon, der_codon))
                        else:
                            codon_changes_nonsyn[(pname, codon_pos)].append((anc_codon, der_codon))


    fs = 15
    plt.figure()
    plt.hist([len(x) for x in codon_changes_syn.values()], bins=range(25), density=True, alpha=0.6, label='synonymous')
    plt.hist([len(x) for x in codon_changes_nonsyn.values()], bins=range(25), density=True, alpha=0.6, label='non-synonymous')
    # plt.yscale('log')
    plt.xlabel('number of mutations per codon',fontsize=fs)
    plt.ylabel('distribution',fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.legend(fontsize=fs)
    plt.savefig('syn_vs_nonsyn.pdf')

    codon_changes_nonsyn_sorted = [x for x in sorted(codon_changes_nonsyn.items(), key=lambda x:len(x[1])) #/(len(x[1]) + len(codon_changes_syn[x[0]])))
                                   if x[0][0]=='VP1']

    plt.figure()
    gene = 'VP1'
    length = 309
    sorted_vals = np.sort([len(codon_changes_nonsyn[(gene, pos)]) for pos in range(length)])
    sorted_indices = np.argsort([len(codon_changes_nonsyn[(gene, pos)]) for pos in range(length)])
    sorted_indices_rev = {v:i for i,v in enumerate(sorted_indices)}
    plt.plot(np.arange(length,0,-1),sorted_vals, label=f'all positions in {gene}')

    BC_positions = [89,91,95,96,97,102]
    plt.scatter([length-sorted_indices_rev[p] for p in BC_positions], [len(codon_changes_nonsyn[(gene, pos)]) for pos in BC_positions], c='C1', label='BC-loop')
    DE_positions = [139,140,141,142,143,144,145,148]
    plt.scatter( [length-sorted_indices_rev[p] for p in DE_positions], [len(codon_changes_nonsyn[(gene, pos)]) for pos in DE_positions], c='C2', label='DE-loop')
    # C_term_positions = range(10)
    # plt.scatter([len(codon_changes_nonsyn[(gene, pos)]) for pos in C_term_positions], [length-sorted_indices_rev[p] for p in C_term_positions], c='C3', label='C-terminus')
    plt.legend(fontsize=fs)
    plt.ylabel('number of non-synonymous mutations',fontsize=fs)
    plt.xlabel(f'variability rank within {gene}',fontsize=fs)
    plt.tick_params(labelsize=0.8*fs)
    plt.savefig('vp1_nonsyn_stats.pdf')

    plt.figure()
    plt.hist([len(x) for n,x in codon_changes_nonsyn.items() if n[0]==gene], bins=range(0,24,3))
    plt.hist([len(codon_changes_nonsyn[(gene, p)]) for p in BC_positions+DE_positions], bins=range(0,24,3))
    plt.hist([len(codon_changes_nonsyn[(gene, p)]) for p in DE_positions], bins=range(0,24,3))

    print(f"pos\tnon-syn\tsyn")
    for p in range(89,103):
        print(f"{p+1}\t{len(codon_changes_nonsyn[('VP1',p)])}\t{len(codon_changes_syn[('VP1',p)])}")

    print(f"pos\tnon-syn\tsyn")
    for p in range(139,150):
        print(f"{p+1}\t{len(codon_changes_nonsyn[('VP1',p)])}\t{len(codon_changes_syn[('VP1',p)])}")
