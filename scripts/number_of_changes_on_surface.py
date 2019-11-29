from Bio import Phylo, SeqIO, AlignIO
import matplotlib.pyplot as plt
import numpy as np
import json


if __name__ == '__main__':
    import argparse
    parser = parser = argparse.ArgumentParser(description='assemble a script for pymol',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tree', type=str, help="input tree (refined)")
    parser.add_argument('--mutations', help="mutation JSON file")
    parser.add_argument('--output', help="name of pymol file")
    args = parser.parse_args()


    proteins = {"VP1":[-12, 297], "VP2":[0,248], "VP3":[0,234]}
    T = Phylo.read(args.tree, 'newick')
    for n in T.find_clades(order='postorder'):
        if n.is_terminal():
            n.ntips = 1
        else:
            n.ntips = np.sum([c.ntips for c in n])

    with open(args.mutations, 'r') as fh:
        tmp = json.load(fh)
        aa_muts = tmp["nodes"]
        annotation =tmp["annotations"]

    surface_exposed_manually={
        "VP1":[(-13,0), (70,73), (74,75), (76,89),(91,92),(93,94),(126,138), (149,160), (203,208), (210,216), (226,232), (257,260), (267,296)],
        "VP2":[(63,64),(65,78), (86,89),(132,157),(215,225),(243,247)],
        "VP3":[(58,80), (87,92), (137,142),(143,145), (179,182),(204,209), (231,234), (234,247)],
        "VP4":[]
    }

    # fix the VP1 numbering. Shouldn't be necessary after nextstrain uses updated annotation
    surface_exposed_manually['VP1'] = [(b+13, e+13) for b,e in surface_exposed_manually['VP1']]

    number_of_changes = {gene: np.zeros((1+annotation[gene]["end"]-annotation[gene]["start"])//3, dtype=int)
                         for gene in annotation}
    for n in T.find_clades():
        if n.name in aa_muts:
            for gene in aa_muts[n.name]["aa_muts"]:
                for m in aa_muts[n.name]["aa_muts"][gene]:
                    if m[0]!='X' and m[-1]!='X' and n.ntips>1:
                        number_of_changes[gene][int(m[1:-1])-1] += 1


    gene_length = {"VP1":309, "VP2": 248, "VP3": 234,"VP4":69}

    fs=16

    fig, axs=plt.subplots(1,1, figsize=(11,2.5))
    offset = 0
    for pi, p in enumerate(['VP4', 'VP2', 'VP3', 'VP1']):
        axs.plot(np.arange(len(number_of_changes[p]))+offset, number_of_changes[p])

        for b,e in surface_exposed_manually[p]:
            axs.plot([b+offset, e+offset], [-.5,-.5], lw=8, c="C%d"%pi)

        offset += gene_length[p]

    plt.ylabel('variability', fontsize=fs)
    plt.xlabel('position along polyprotein', fontsize=fs)
    plt.text(0,-.9, 'outer surface:', fontsize=fs*0.7)
    plt.tick_params(labelsize=fs*0.8)
    plt.tight_layout()
    plt.savefig(args.output)

