from Bio import Phylo, SeqIO, AlignIO
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from collections import defaultdict
import numpy as np
import json

from treetime.treeanc import mutations


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

    gene_length = {"VP1":309, "VP2": 248, "VP3": 234,"VP4":69}
    length_surface = {y: sum([x[1]-x[0] for x in surface_exposed_manually[y]]) for y in surface_exposed_manually}
    mutations_surface = defaultdict(int)
    mutations_other = defaultdict(int)
    colors = {"VP1": 'purple', "VP2":"orange", "VP3":"green", "VP4":"blue"}

    # fix the VP1 numbering. Shouldn't be necessary after nextstrain uses updated annotation
    surface_exposed_manually['VP1'] = [(b+13, e+13) for b,e in surface_exposed_manually['VP1']]

    number_of_changes = {gene: np.zeros((1+annotation[gene]["end"]-annotation[gene]["start"])//3, dtype=int)
                         for gene in annotation}
    for n in T.find_clades():
        if n.name in aa_muts:
            for gene in aa_muts[n.name]["aa_muts"]:
                for m in aa_muts[n.name]["aa_muts"][gene]:
                    if m[0]!='X' and m[-1]!='X' and n.ntips>1:
                        pos = int(m[1:-1])
                        number_of_changes[gene][pos-1] += 1
                        if gene in surface_exposed_manually and any([(x[0]<=pos and pos<=x[1]) for x in surface_exposed_manually[gene]]):
                            mutations_surface[gene] += 1
                        elif gene in surface_exposed_manually:
                            mutations_other[gene] += 1



    fs=16

    fig, axs=plt.subplots(1,1, figsize=(11,2.5))
    offset = 0
    for pi, p in enumerate(['VP4', 'VP2', 'VP3', 'VP1']):
        if p=="VP1":
            axs.add_patch(Rectangle((89+offset,0), 14, 10, ec='silver', fc='silver', linewidth=2))
            axs.add_patch(Rectangle((140+offset,0), 13, 10, ec='silver', fc='silver', linewidth=2))
            axs.add_patch(Rectangle((279+offset,0), 30, 10, ec='silver', fc='silver', linewidth=2))
        axs.plot(np.arange(len(number_of_changes[p]))+offset, number_of_changes[p], c=colors[p])
        plt.text(0.5*(2*offset+len(number_of_changes[p])) + (50 if p=='VP1' else 0),7, p,
                fontsize=fs,horizontalalignment='center')

        for b,e in surface_exposed_manually[p]:
            axs.plot([b+offset, e+offset], [-.5,-.5], lw=8, c=colors[p])#"C%d"%pi)

        offset += gene_length[p]

    plt.ylabel('variability', fontsize=fs)
    plt.xlabel('position along polyprotein', fontsize=fs)
    plt.text(0,-.9, 'outer surface:', fontsize=fs*0.7)
    plt.tick_params(labelsize=fs*0.8)
    plt.tight_layout()
    plt.savefig(args.output)

    from scipy.stats import fisher_exact
    for y in mutations_surface:
        print(f"{y}: mutation per site surface: {mutations_surface[y]/length_surface[y]:1.2f}")
        print(f"{y}: mutation per site other: {mutations_other[y]/(gene_length[y] - length_surface[y]):1.2f}")

        OR, pval = fisher_exact([[mutations_surface[y], length_surface[y]], [mutations_other[y], gene_length[y] - length_surface[y]]])
        print(f"{y}: OR={OR:1.3f}, p={pval:1.3e}")

    all_muts_surface = sum(list(mutations_surface.values()))
    all_surface = sum(list(length_surface.values()))
    all_muts_other = sum(list(mutations_other.values()))
    all_length_other = sum([gene_length[y] - length_surface[y] for y in gene_length])
    print(f"VP1-VP4: mutation per site surface: {all_muts_surface/all_surface:1.2f}")
    print(f"VP1-VP4: mutation per site other: {all_muts_other/all_length_other:1.2f}")

    OR, pval = fisher_exact([[all_muts_surface, all_surface], [all_muts_other, all_length_other]])
    print(f"VP1-VP4: OR={OR:1.3f}, p={pval:1.3e}")
