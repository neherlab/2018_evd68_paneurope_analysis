
import os, sys
import re
import time
import numpy as np
from Bio import Phylo
from collections import defaultdict
from augur.utils import read_metadata, read_node_data, write_json, read_config, read_lat_longs, read_colors, json_to_tree
from augur.export_v2 import parse_node_data_and_metadata
import pandas as pd

if __name__ == '__main__':
    import argparse

    parser = parser = argparse.ArgumentParser(description='track lineages for figure plot',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tree', help="input tree (refined)")           
    parser.add_argument('--branch-lengths', help="branch length JSON")
    parser.add_argument('--meta', help="metadata TSV")
    parser.add_argument('--output', help="output figure PNG")
    args = parser.parse_args()

    treefile = args.tree
    branchfile = args.branch_lengths
    metadatafile = args.meta

    sampSets = {'2014-5':{'startTime':2014, 'endTime':2015, 'time': [], 'lineages': []},
                '2016-7':{'startTime':2016, 'endTime':2017, 'time': [], 'lineages': []},
                '2018-9':{'startTime':2018, 'endTime':2019, 'time': [], 'lineages': []},
                '2011-9':{'startTime':2011, 'endTime':2019, 'time': [], 'lineages': []}}

    origKeys = list(sampSets.keys())

    for key in origKeys:
        realKey = key

        params = sampSets[realKey]

        startTime = params['startTime']
        endTime = params['endTime']
        time = params['time']
        num_lineages = params['lineages']

        T = Phylo.read(treefile, 'newick')
        node_data = read_node_data([branchfile])
        node_data, node_attrs, node_data_names, metadata_names = parse_node_data_and_metadata(T, [branchfile], metadatafile)
        rate = node_data['clock']['rate']

        for node in T.find_clades(order='postorder'):
            data = node_data['nodes'][node.name]
            node.date = data['date']
            node.num_date = data['numdate']
            raw_data = node_attrs[node.name]
            node.region = raw_data['region'] if 'region' in raw_data else ''
            node.branch_length = data['branch_length']/rate

        #set parents to avoid excess tree-traversal
        for node in T.find_clades(order='preorder'):
            for child in node:
                child.parent = node

        #including or not 'future_tips' doesn't make much difference.
        all_tips = []
        tips_current = {}
        future_tips = {}
        for node in T.find_clades(terminal=True, order='postorder'):
            all_tips.append(node.name)
            #Remove this node as it is 2016-01-01 (estimated) and from B1 clade,
            #so adds a whole other deep coalescence to the graph!!
            if node.name == "US/MO/14-18949":
                continue

            if node.num_date >= startTime and node.num_date < endTime: 
                tips_current[node.name] = node
            elif node.num_date >= endTime:
                future_tips[node.name] = node

        removed = []

        tips_to_remove = [ti for ti in all_tips if ti not in removed + list(tips_current.keys())]

        #remove tips not being inspected right now
        for tip in tips_to_remove:
            T.prune(tip)

        #set parents again to avoid excess tree-traversal
        for node in T.find_clades(order='preorder'):
            for child in node:
                child.parent = node

        #New approach - set first 'count' time to be latest tip (date nearest to present)
        latestDate = sorted([v.num_date for v in tips_current.values()], reverse=True)[0]
        time = [latestDate]
        num_lineages = [len(tips_current)]

        Phylo.write(T, "results/{}_vp1_tree.nwk".format(realKey), "newick")