'''
Calculates the number of changes that happen along the tree to tips from
different clades and plots this as a function of time for different clades
'''
from Bio import Phylo
import matplotlib.pyplot as plt
import json
import numpy as np



regions = {
    "BC":("VP1", 89,103),
    "DE":("VP1", 139,148),
    "cTerm": ("VP1", 279,309),
    "VP2_central":("VP2", 134,156),
}

colors = {
    "B3":"#DC2F24", "B2":"#E67932", "A2":"#3E58CF", "A1":"#67A0CB"
    }

def count_mutations(T,muts, gene, start, stop):
    T.root.mut_count = 0
    for n in T.get_nonterminals(order='preorder'):
        for c in n:
            nmuts = 0
            if c.name in muts and gene in muts[c.name]['aa_muts']:
                for m in muts[c.name]['aa_muts'][gene]:
                    a,p,d = m[0], int(m[1:-1])-1, m[-1]
                    if a!='X' and d!='X' and p>=start and p<stop:
                        nmuts+=1
            c.mut_count = n.mut_count + nmuts


if __name__ == '__main__':
    import argparse
    parser = parser = argparse.ArgumentParser(description='assemble a script for pymol',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tree', type=str, help="input tree (refined)")
    parser.add_argument('--mutations', help="json with mutations")
    parser.add_argument('--clades', help="JSON file with clades")
    parser.add_argument('--branch-lengths', help="JSON file with branch-length")
    parser.add_argument('--output', help="file name to write the alignemt to")
    args = parser.parse_args()

    T = Phylo.read(args.tree, 'newick')

    with open(args.clades, 'r') as fh:
      clades = json.load(fh)["nodes"]

    with open(args.mutations, 'r') as fh:
      muts = json.load(fh)["nodes"]

    with open(args.branch_lengths, 'r') as fh:
      numdates = json.load(fh)["nodes"]

    clades_to_label = ['B3', 'B2', 'A1', 'A2']
    time_points = np.linspace(1990,2019,101)
    traj = {}
    for region in regions:
        count_mutations(T, muts, *regions[region])
        traj[region] = {x:[[] for t in time_points] for x in clades_to_label}
        for n in T.get_terminals():
            clade = clades[n.name]['clade_membership']
            if clade in clades_to_label:
                path_to_root = [T.root] + T.get_path(target=n)
                mut_counts = np.array([[0,0]] + [[numdates[x.name]['numdate'], x.mut_count] for x in path_to_root])
                tmp_trajs = []
                for ti,t in enumerate(time_points):
                    if t>mut_counts[-1][0]:
                        break

                    nmuts = mut_counts[mut_counts[:,0].searchsorted(t),1]
                    traj[region][clade][ti].append(nmuts)


    fs=16
    fig, axs = plt.subplots(3,1, figsize=(7.5,11)) #, sharex=True)
    plt.tight_layout()
    [axs[0].axvline(xs, color="gainsboro", linestyle="-") for xs in [1990,1995,2000,2005,2010,2015,2020]]
    [axs[0].axhline(xs, color="silver", linestyle=":") for xs in [0,1,2,3,4,5,6]]
    [axs[1].axvline(xs, color="gainsboro", linestyle="-") for xs in [1990,1995,2000,2005,2010,2015,2020]]
    [axs[1].axhline(xs, color="silver", linestyle=":") for xs in [0,1,2,3]]
    [axs[2].axvline(xs, color="gainsboro", linestyle="-") for xs in [1990,1995,2000,2005,2010,2015,2020]]
    [axs[2].axhline(xs, color="silver", linestyle=":") for xs in [0,2,4,6,8]]

    for ci,clade in enumerate(clades_to_label):
        nmuts = 0.08*(ci-1.5) + np.array([np.mean(x) for x in traj['BC'][clade]])
        axs[0].plot(time_points, nmuts, label=clade, lw=3, c=colors[clade])

    axs[0].set_ylabel('Cumulative mutations\nin the BC loop', fontsize=fs)
    for ci,clade in enumerate(clades_to_label):
        nmuts = 0.08*(ci-1.5) + np.array([np.mean(x) for x in traj['DE'][clade]])
        axs[1].plot(time_points, nmuts, label=clade, lw=3, c=colors[clade])

    axs[1].set_yticks([0,1,2,3])
    axs[1].set_ylabel('Cumulative mutations\nin the DE loop', fontsize=fs)
    for ci,clade in enumerate(clades_to_label):
        nmuts = 0.08*(ci-1.5) + np.array([np.mean(x) + np.mean(y) for x,y in zip(traj['BC'][clade],traj['DE'][clade])])
        axs[2].plot(time_points, nmuts, label=clade, lw=3, c=colors[clade])

    axs[2].set_ylabel('Sum of BC and DE loop', fontsize=fs)
    axs[2].set_xlabel('Time')
    axs[0].legend(fontsize=fs)
    for ax in axs:
        ax.tick_params(labelsize=0.8*fs)
    plt.tight_layout()
    plt.savefig(args.output)