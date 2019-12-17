from Bio import Phylo, SeqIO, AlignIO
from matplotlib import cm
import numpy as np
import json



if __name__ == '__main__':
    import argparse
    parser = parser = argparse.ArgumentParser(description='assemble a script for pymol',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tree', type=str, help="input tree (refined)")
    parser.add_argument('--alignment', help="alignment")
    parser.add_argument('--clades', help="JSON file with clades")
    parser.add_argument('--output', help="file name to write the alignemt to")
    args = parser.parse_args()

    T = Phylo.read(args.tree, 'newick')
    aln = {s.name:s for s in AlignIO.read(args.alignment, 'fasta')}
    aa_ref_seq = aln[T.root.name]

    with open(args.clades, 'r') as fh:
      clades = json.load(fh)["nodes"]

    clade_to_node = {v['clade_annotation']:k for k,v in clades.items() if 'clade_annotation' in v}

    A1_rep = '692-OsakaCity-JPN-2013'
    A2_rep = 'EVD68_NLD_007_181008_NFLG'
    A2_early = '2123800027'
    B1_rep = 'EV-D68/Homosapiens/USA/SSENT38/2014'
    B3_US_rep = 'EVD68_BEL_006_181010_NFLG'
    B3_ESP_rep = 'EVD68_ESP_007_180925_NFLG'

    # print sequence difference in the epitopes
    vp2_epi = [134,156]
    vp2_root = np.array(list(aa_ref_seq[vp2_epi[0]:vp2_epi[1]]))

    with open(args.output, 'w') as fh:
        outstr = f"root\t{aa_ref_seq.seq[vp2_epi[0]:vp2_epi[1]]}"
        fh.write(outstr+'\n')
        for n, seq in [('A2', aln[A2_rep]), ('A1', aln[A1_rep]),
                       ('A2early', aln[A2_early]),
                       ('B-base', aln[clade_to_node['B']]),
                       ('B1', aln[B1_rep]),
                       ('B3-US', aln[B3_US_rep]), ('B3-EU', aln[B3_ESP_rep])]:
            vp2_seq = np.array(list(seq.seq[vp2_epi[0]:vp2_epi[1]]))
            vp2_seq[vp2_seq==vp2_root] = '.'
            outstr = f"{n}\t{''.join(vp2_seq)}"
            print(outstr)
            fh.write(outstr + '\n')

