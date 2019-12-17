import argparse
from Bio import Phylo, SeqIO, AlignIO
from matplotlib import cm
import numpy as np
import json


if __name__ == '__main__':
    parser = parser = argparse.ArgumentParser(description='assemble a script for pymol',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--reference', type=str, help="reference genome of VP1")
    parser.add_argument('--tree', type=str, help="input tree (refined)")
    parser.add_argument('--clades', help="clades JSON file")
    parser.add_argument('--mutations', help="mutation JSON file")
    parser.add_argument('--output', help="name of pymol file")
    args = parser.parse_args()

    ref_seq = SeqIO.read(args.reference, 'genbank')

    T = Phylo.read(args.tree, 'newick')
    # calculate number of tips for every node
    for n in T.find_clades(order='postorder'):
        if n.is_terminal():
            n.ntips = 1
        else:
            n.ntips = np.sum([c.ntips for c in n])

    with open(args.mutations, 'r') as fh:
        tmp = json.load(fh)
        aa_muts = tmp["nodes"]
        annotation =tmp["annotations"]

    with open(args.clades, 'r') as fh:
        clades = json.load(fh)["nodes"]

    # map clade names to the internal node that defines the clade
    clade_to_node = {v['clade_annotation']:k for k,v in clades.items() if 'clade_annotation' in v}

    # calculate the number of amino acid mutations at each position in each gene.
    # conditions on mutations that don't happen on terminal nodes
    number_of_changes = {gene: np.zeros((1+annotation[gene]["end"]-annotation[gene]["start"])//3, dtype=int) for gene in annotation}

    for n in T.find_clades():
        if n.name in aa_muts:
            for gene in aa_muts[n.name]["aa_muts"]:
                for m in aa_muts[n.name]["aa_muts"][gene]:
                    if m[0]!='X' and m[-1]!='X' and n.ntips>1:
                        number_of_changes[gene][int(m[1:-1])-1] += 1

    cmap = cm.inferno

    # start to assemble pymol string of coloring of by number of differences
    diff_color_str = ""
    chain_start = {"VP1":12, "VP2":0, "VP3":0, "VP4":0}
    chain_names = {"VP1":"4wm8_0003 and chain A", "VP2":"4wm8_0003 and chain B",
                   "VP3":"4wm8_0003 and chain C", "VP4":"4wm8_0003 and chain D"}
    background_colors= {"VP1": "lightpink", "VP2":"paleyellow", "VP3":"palegreen", "VP4":"gray50"}

    # assign colors to chains and highlight changing positions
    for gene, cs in chain_start.items():
        for pos in range(cs-1,len(number_of_changes[gene])):
            pseq = pos-cs+1
            chain =chain_names[gene]
            if number_of_changes[gene][pos]>2:
                col = "firebrick"
            elif number_of_changes[gene][pos]:
                col = "ruby"
            else:
                col = background_colors[gene]

            diff_color_str += f"sele res_{gene}_{pseq}, {chain} and resi {pseq}\n"
            # diff_color_str += f"show surface, res_{gene}_{pseq}\n"
            diff_color_str += f"color {col}, res_{gene}_{pseq}\n"

    # handle changes at the beginning of VP1 separately, since they are actually in
    # VP3 -- this won't be necessary in the future after the annotation used in nextstrain is update
    for pos in range(0,13):
        gene='VP1'
        chain = chain_names['VP3']
        pseq = pos + 235
        if number_of_changes[gene][pos]>1:
            col = "firebrick"
        elif number_of_changes[gene][pos]:
            col = "ruby"
        else:
            col = background_colors['VP3']

        diff_color_str += f"sele res_{gene}_{pseq}, {chain} and resi {pseq}\n"
        # diff_color_str += f"show surface, res_{gene}_{pseq}\n"
        diff_color_str += f"color {col}, res_{gene}_{pseq}\n"

    # delete all monomers that aren't used
    delete_str = "\n".join([f"delete 4wm8_{i:04d}" for i in range(6,61)])

    chain_start = 12
    BC_loop = f"{90-chain_start}-{103-chain_start}"
    DE_loop = f"{140-chain_start}-{153-chain_start}"
    EVd68_epi = f"{282-chain_start}-{308-chain_start}"


    pymol_script = \
f"""
fetch 4wm8, type=pdb1, multiplex=1
split_state 4wm8
{delete_str}

orient

hide everything

sele vp1_1, 4wm8_0001 and chain A
sele vp2_1, 4wm8_0001 and chain B
sele vp3_1, 4wm8_0001 and chain C
sele vp4_1, 4wm8_0001 and chain D

sele vp1_2, 4wm8_0002 and chain A
sele vp2_2, 4wm8_0002 and chain B
sele vp3_2, 4wm8_0002 and chain C
sele vp4_2, 4wm8_0002 and chain D

sele vp1_3, 4wm8_0003 and chain A
sele vp2_3, 4wm8_0003 and chain B
sele vp3_3, 4wm8_0003 and chain C
sele vp4_3, 4wm8_0003 and chain D

sele vp1_4, 4wm8_0004 and chain A
sele vp2_4, 4wm8_0004 and chain B
sele vp3_4, 4wm8_0004 and chain C
sele vp4_4, 4wm8_0004 and chain D

sele vp1_5, 4wm8_0005 and chain A
sele vp2_5, 4wm8_0005 and chain B
sele vp3_5, 4wm8_0005 and chain C
sele vp4_5, 4wm8_0005 and chain D

sele BC_loop, 4wm8_0005 and chain A and resi {BC_loop}
sele DE_loop, 4wm8_0005 and chain A and resi {DE_loop}
sele EVd68_epi, 4wm8_0005 and chain A and resi {EVd68_epi}

show surface, 4wm8_0001
show surface, 4wm8_0002
show surface, 4wm8_0003
show surface, 4wm8_0004
show surface, 4wm8_0005

color gray40, 4wm8_0002
color gray40, 4wm8_0004
color gray80, 4wm8_0005

color teal, BC_loop
color teal, DE_loop
color teal, EVd68_epi

color deeppurple, vp1_1
color sand, vp2_1
color forest, vp3_1
color gray, vp4_1

{diff_color_str}

"""

    with open(args.output, 'w') as fh:
        fh.write(pymol_script)
