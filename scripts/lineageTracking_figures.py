## Lineage tracking - Exploration of scripts
# FOR WHOLE GENOME RUN!

## For more code, including diff ways to plot, look at original
#code (and earlier versions) in entero/scripts/finalLineageScripts

import os, sys
import re
import time
import numpy as np
from Bio import Phylo
from collections import defaultdict
from augur.utils import read_metadata, read_node_data, write_json, read_config, read_lat_longs, read_colors, json_to_tree
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as patches
from textwrap import wrap
from augur.export_v2 import parse_node_data_and_metadata
import pandas as pd
from scipy.stats import iqr
import matplotlib.ticker as ticker

if __name__ == '__main__':
    import argparse

    parser = parser = argparse.ArgumentParser(description='track lineages for figure plot',
                formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tree', help="input tree (refined)")           
    parser.add_argument('--branch-lengths', help="branch length JSON")
    parser.add_argument('--clades', help="clades JSON file")
    parser.add_argument('--meta', help="metadata TSV")
    parser.add_argument('--output', help="output figure PNG")
    args = parser.parse_args()

    #treefile = "results/tree_genome.nwk"
    #branchfile = "results/branch_lengths_genome.json"
    #cladefile = "results/clades_genome.json"
    #metadatafile = "results/metadata-ages.tsv"

    treefile = args.tree
    branchfile = args.branch_lengths
    cladefile = args.clades
    metadatafile = args.meta

    print("treefile", treefile)
    print("branchfile", branchfile)
    print("cladefile", cladefile)
    print("metafile", metadatafile)

    #T = Phylo.read(treefile, 'newick')
    #node_data = read_node_data([branchfile, cladefile])
    node_data = read_node_data([cladefile])


    sampSets = {'2014-5':{'startTime':2014, 'endTime':2015, 'time': [], 'lineages': []},
                '2016-7':{'startTime':2016, 'endTime':2017, 'time': [], 'lineages': []},
                '2018-9':{'startTime':2018, 'endTime':2019, 'time': [], 'lineages': []},
                '2011-9':{'startTime':2011, 'endTime':2019, 'time': [], 'lineages': []}}

    EUUStips = {}
    allRegionTips = {}
    all2014Tips = {}
    all2016Tips = {}
    all2018Tips = {}

    origKeys = list(sampSets.keys())

    for key in origKeys:

        for i in range(0,2):
            if i==0:
                realKey = key
            else:
                realKey = key
                realKey = key+'_west'
                tmpParam = sampSets[key]
                sampSets[realKey] = {'startTime':tmpParam['startTime'], 
                                'endTime':tmpParam['endTime'], 'time': [], 'lineages': []}

            params = sampSets[realKey]

            startTime = params['startTime']
            endTime = params['endTime']
            time = params['time']
            num_lineages = params['lineages']

            T = Phylo.read(treefile, 'newick')
            node_data = read_node_data([branchfile, cladefile])
            #raw_strain_info = collect_strain_info(node_data, metadatafile)
            node_data, node_attrs, node_data_names, metadata_names = parse_node_data_and_metadata(T, [branchfile, cladefile], metadatafile)
            rate = node_data['clock']['rate']

            for node in T.find_clades(order='postorder'):
                data = node_data['nodes'][node.name]
                node.clade_membership = data['clade_membership']
                node.date = data['date']
                node.num_date = data['numdate']
                #raw_data = raw_strain_info[node.name]
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
                if i==0:
                    if node.num_date >= startTime and node.num_date < endTime: 
                        tips_current[node.name] = node
                    elif node.num_date >= endTime:
                        future_tips[node.name] = node
                else:
                    if node.num_date >= startTime and node.num_date < endTime and node.region in ['north_america', 'europe']: 
                        tips_current[node.name] = node
                    elif node.num_date >= endTime and node.region in ['north_america', 'europe']:
                        future_tips[node.name] = node

            removed = []
            # This counts lineages that extend beyond the window as a 'tip'
            # we don't want to do this!
        #    for name, node in future_tips.items():
        #        collapse_until = node
        #        while collapse_until.parent.num_date > endTime:
        #            collapse_until = collapse_until.parent

        #        if collapse_until == node: #coalesce is inside or before timezone, trim and leave.
        #            node.num_date = endTime
        #            tips_current[node.name] = node
        #        else: #coalesce is after timezone, 'trim' coalesce and add
        #            collapse_until.collapse_all()
        #            while collapse_until.count_terminals() > 1:
        #                removed.append(collapse_until[0].name)
        #                T.prune(collapse_until[0].name)
        #            collapse_until[0].num_date = endTime
        #            tips_current[collapse_until[0].name] = collapse_until[0]

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

            Phylo.write(T, "results/{}_genome_tree.nwk".format(realKey), "newick")

            ### LOOP START
            if realKey=='2011-9_west':
                EUUStips = tips_current.copy()
            if realKey=='2011-9':
                allRegionTips = tips_current.copy()
            if realKey=='2014-5':
                all2014Tips = tips_current.copy()
            if realKey=='2016-7':
                all2016Tips = tips_current.copy()
            if realKey=='2018-9':
                all2018Tips = tips_current.copy()

            while len(tips_current) > 1:
                parents = defaultdict(list)
                for name, tip in tips_current.items():
                    parents[tip.parent.num_date].append(tip.parent)
                    parents[tip.parent.num_date] = list(set(parents[tip.parent.num_date]))

                earliestParentTime = sorted(parents.keys(), reverse=True)[0]

                earliestParents = parents[earliestParentTime]
                if key == "2016-7" and "NODE_0001300" in tips_current:
                    print("earliestPar is {} and node still present".format(earliestParents))

                brkrly = False
                for par in earliestParents:
                    for child in par:
                        if child.name not in tips_current:
                            print("Child before: {}, parent: {}".format(child.num_date, par.num_date))
                            child.num_date = par.num_date + 0.000000000001
                            print("Child after: {}, parent: {}".format(child.num_date, par.num_date))
                            brkrly = True
                            print("WARNING! Negative branch length for {}? ({})".format(child.name, child.branch_length))

                if brkrly:
                    continue
                
                for par in earliestParents:
                    for child in par:
                        del tips_current[child.name]
                    tips_current[par.name] = par

                time.append(earliestParentTime)
                num_lineages.append(len(tips_current))

            
            sampSets[realKey]['time'] = time
            sampSets[realKey]['lineages'] = num_lineages

    for key,val in sampSets.items():
        #print("For years {} to {}:".format(val['startTime'], val['endTime']))
        print("For data {}:".format(key))
        print("\tNum data points: {}".format(len(val['lineages'])))
        print("\tMax lineages: {}".format(max(val['lineages'])))



    allDats = [val.num_date for key,val in EUUStips.items()]
    #allDats = [val.num_date for key,val in allRegionTips.items()]
    othDats = [val.num_date for key,val in allRegionTips.items() if val.region not in ['north_america', 'europe']]
    worldDates = [val.num_date for key,val in allRegionTips.items()]

    tipDat2014 = [val.num_date-2015 for key,val in all2014Tips.items()]
    tipDat2016 = [val.num_date-2017 for key,val in all2016Tips.items()]
    tipDat2018 = [val.num_date-2019 for key,val in all2018Tips.items()]


    linestyle = {"2014-5": '-',
                "2014-5_west": '--',
                "2016-7": '-',
                "2016-7_west": '--',
                "2018-9": '-',
                "2018-9_west": '--',
                "2011-9": '-',
                "2011-9_west": '--'}

    clrs = {"2014-5": 'C0',
                "2014-5_west": 'C0',
                "2016-7": 'C1',
                "2016-7_west": 'C1',
                "2018-9": 'C2',
                "2018-9_west": 'C2',
                "2011-9": 'C3',
                "2011-9_west": 'C3'}

    for epis in ['2014-5', '2016-7', '2018-9']:
        #epis = '2014-5'
        endT = sampSets[epis]['endTime']
        yrsB2 = [abs(x-(endT-2)) for x in sampSets[epis]['time']]
        actYr2 = sampSets[epis]['time'][yrsB2.index(min(yrsB2))]
        ln2Yrs = sampSets[epis]['lineages'][yrsB2.index(min(yrsB2))]
        yrsB4 = [abs(x-(endT-4)) for x in sampSets[epis]['time']]
        actYr4 = sampSets[epis]['time'][yrsB4.index(min(yrsB4))]
        ln4Yrs = sampSets[epis]['lineages'][yrsB4.index(min(yrsB4))]

        rlt = "{}:\n\t 2 yrs before:\t{}\t{}\n\t 4 yrs before:\t{}\t{}".format(epis, actYr2,ln2Yrs,actYr4,ln4Yrs)
        print(rlt)


    #import ipdb; ipdb.set_trace()

    sliding_window = True
    window_size = 2 #in months - 1 month or 2 months. If sliding window, will slide by 1/2 month and 1 month, respectively

    plt.figure(1, figsize=(7,5.5))
    ax1 = plt.subplot(211)
    [ax1.axvline(xs, color="silver") for xs in list(range(-10,1))]
    [ax1.axhline(xs, color="silver", linestyle=":") for xs in [1, 5,10,50,100]]

    #add the IQR of sample dates for all 3 epidemics combined
    tipDates = np.concatenate((tipDat2014, tipDat2016, tipDat2018))
    iqrVal = iqr(tipDates)
    lef = np.mean(tipDates)-(iqrVal/2)
    rig = np.mean(tipDates)+(iqrVal/2)
    rightmostTime = 0
    ax1.axvline(lef-rightmostTime, color="black", linestyle="--")
    ax1.axvline(rig-rightmostTime, color="black", linestyle="--")

    ax4 = plt.subplot(212)
    [ax4.axvline(xs, color="silver") for xs in list(range(2013,2020))]
    [ax4.axvline(xs, color="silver") for xs in list(range(-10,1))]
    [ax4.axhline(xs, color="silver", linestyle=":") for xs in [0,5,10,15]]

    #plt.subplots_adjust(top=0.99, bottom=0.05, hspace=0.3)
    plt.subplots_adjust(top=0.95, bottom=0.05, hspace=0.32, left=0.135)
    #plt.tight_layout()

    for key, val in sampSets.items():
        if "_west" not in key:
            yrs = val['time'] 
            num_lin = val['lineages']
            maxyr = max(yrs)
            maxln = max(num_lin)

            sTime = val['startTime']
            eTime = val['endTime']

            lab = key

            if key != "2011-9":

                ax1.step([x-eTime for x in yrs], [x for x in num_lin],
                    c=clrs[key], linestyle=linestyle[key],
                    where='post', label=lab)

            if key != "2011-9":
                # make evenly spaced bins for all years 2013-2019 (12 per year or 6 per year, depending)
                # if sliding window, sliding window bins will be 1/2 month (for 1 month) or 1 month off (for 2 months)
                num_bins = 72
                slide_size = 0.0415
                if window_size == 2:
                    num_bins = 36
                    slide_size = 0.0833

                histdat = np.histogram(yrs, bins=num_bins, range=(eTime-6,eTime))
                bins2 = [ij-slide_size for ij in histdat[1]]

                #put the yrs data (x axis of 2nd panel) into the bins
                ret = pd.cut(x=yrs, bins=histdat[1])
                #put the sample dates into the bins
                smps = pd.cut(x=worldDates, bins=histdat[1])

                ret2 = pd.cut(x=yrs, bins=bins2)
                smps2 = pd.cut(x=worldDates, bins=bins2)

                cats = ret.categories.copy()
                allcats = cats
                if sliding_window:
                    allcats = cats.union(ret2.categories).sort_values()

                # for each bin, find out how many samples it has,
                # and find out the slope (change in num lineages / change in time)
                slp = {}
                for intr in ret.categories:
                    indx = np.where(ret==intr)[0]
                    rr = 0
                    old_rr = 0
                    rise = 0
                    run = 0
                    old_run = 0
                    if len(indx) == 1:
                        rise = 1
                        run = intr.left-intr.right
                        rr = rise/run

                    if len(indx) != 0 and len(indx) != 1:
                        ys = [yrs[j] for j in indx]
                        lins = [num_lin[j] for j in indx]
                        old_run = min(ys)-max(ys)
                        run = intr.left-intr.right
                        rise = min(lins)-max(lins)
                        rr = rise/run
                        old_rr = rise/old_run

                    slp[intr] = {'slope': rr, 'old_slope': old_rr, 'samps': len(np.where(smps==intr)[0]),
                        'lin_change':rise, 'time_change': run, 'old_time_change':old_run}

                if sliding_window:
                    for intr in ret2.categories:
                        indx = np.where(ret2==intr)[0]
                        rr = 0
                        old_rr = 0
                        rise = 0
                        run = 0
                        old_run = 0
                        if len(indx) != 0 and len(indx) != 1:
                            ys = [yrs[j] for j in indx]
                            lins = [num_lin[j] for j in indx]

                            old_run = min(ys)-max(ys)
                            run = intr.left-intr.right
                            rise = min(lins)-max(lins)
                            rr = rise/run
                            old_rr = rise/old_run
                        slp[intr] = {'slope': rr, 'old_slope': old_rr, 'samps': len(np.where(smps2==intr)[0]),
                            'lin_change':rise, 'time_change':run, 'old_time_change':old_run}

                #put these values into arrays so easy to plot
                xvals = []
                xvalspts = []
                slope_y = []
                samps_y = []
                othsamps_y = []
                linchange_y = []
                #for key, val in slp.items():
                for keyk in allcats:
                    val = slp[keyk]
                    xvals.append(np.mean([keyk.left, keyk.right]))
                    xvalspts.append((keyk.left, keyk.right))
                    slope_y.append(val['slope'])
                    samps_y.append(val['samps'])
                # othsamps_y.append(val['oth_samps'])
                    linchange_y.append(val['lin_change'])

                #plot
                xvals2 = [x for x in xvals if x>=(maxyr-2)]
                linchange_y2 = [y for y,x in zip(linchange_y,xvals) if x>=(maxyr-2)]
                maxlin1 = max(num_lin)
                ax4.plot( [x for x in xvals2], [abs(k)/maxlin1/2*100 for k in linchange_y2], #c=clrs[key], 
                    c="black", label=lab) 
                if lab == "2018-9":
                    slp_2018 = slp.copy()


    startAx1 = ax1.get_xlim()
    fixedAx3 = (1992.7989979662855, 2020.2476667635106)

    # panel 1 info
    ax1.set_yscale('log')
    ax1.set_ylabel('\n'.join(wrap("Number of Lineages Remaining (log)",20)))
    ax1.legend(loc="upper left")
    ax1.set_xlim((-6,startAx1[1]-1))
    ax1.set_xlabel("Years before end of epidemic year (2019, 2017, 2015)")
    ax1a = ax1.twinx()
    ax1a.set_yscale('log')
    ax1a.set_ylim(ax1.get_ylim())

    ax1.text(-7.1,200, "B", fontsize=50)
    ax4.text(2011.9, 14, "C", fontsize=50)

    #panel 4 info
    ws = " per Month"
    #ax4.set_ylabel('\n'.join(wrap("Percent of Mean Total Lineages Coalescing{}".format(ws),25)))
    ax4.set_ylabel("Percent of Mean\nTotal Lineages\nCoalescing per Month")
    ax4.set_xlim((fixedAx3[1]-(startAx1[1]-(-6)),fixedAx3[1]-1))
    ax4a = ax4.twinx()
    ax4a.plot( xvals, samps_y, c="C4", label="Samples"), #linestyle=(0,(.8,.8)), )
    ax4a.set_ylabel('\n'.join(wrap("Mean Number of Samples per Month".format(ws),20)), color="C4")
    ax4a.tick_params(axis='y', labelcolor="C4")

    #plt.show()
    #plt.savefig("../figures/trialplot.png")
    plt.savefig(args.output)
    plt.close()






