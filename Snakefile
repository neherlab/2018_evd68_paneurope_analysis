# This should include rules for making the figures for the
# Europe EV-D68 2018 paper

rule clade_age_fig:
    input:
        auspice_tree = "../enterovirus_d68/vp1/results/metadata_subgenotype_2018y.tsv"
    output:
        fig = "figures/age_clade_plot-fig3.pdf",
        suppfig_bootstrap = "figures/supp-BootstrapAgeDist.pdf",
        suppfig_perc_age = "figures/supp-age_by_clade.pdf"
    shell:
        """
        Rscript scripts/age_clade_figure.R {input.auspice_tree} {output.fig} \
            {output.suppfig_bootstrap} {output.suppfig_perc_age}
        """

# This also creates the genome trees, with only tips inside epidemic periods.
rule lineage_tracking_fig:
    input:
        treefile = "../enterovirus_d68/genome/results/tree_2018y.nwk",
        branchfile = "../enterovirus_d68/genome/results/branch_lengths_2018y.json",
        cladefile = "../enterovirus_d68/genome/results/clades_2018y.json",
        metafile = "../enterovirus_d68/genome/results/metadata-ages.tsv"
    output:
        fig = "figures/lineage_panels-fig2.pdf"
    shell:
        """
        python scripts/lineageTracking_figures.py \
            --tree {input.treefile} --branch-lengths {input.branchfile} \
            --clades {input.cladefile} --meta {input.metafile} \
            --output {output.fig}
        """

rule supp_sample_date_fig: #INPUT FILE MAY NEED UPDATING IN FUTURE!
    input:
        swedish_meta = "../enterovirus_d68/data/20190611_Karolinska-region.csv"
    output:
        fig = "figures/supp-sampleHist.pdf"
    shell:
        """
        Rscript scripts/2018_sampleDate_hist.r {input.swedish_meta} {output.fig} 
        """

rule pymol_script:
    input:
        treefile = "../enterovirus_d68/genome/results/tree_2018y.nwk",
        cladefile = "../enterovirus_d68/genome/results/clades_2018y.json",
        mutationfile = "../enterovirus_d68/genome/results/aa_muts_2018y.json",
        ref = "../enterovirus_d68/genome/config/ev_d68_reference_genome.gb"
    output:
        pymol_script = "results/pymol_script_vp1-vp4.pml"
    shell:
        """
        python scripts/make_pymol_file.py --tree {input.treefile} \
            --clades {input.cladefile} --mutations {input.mutationfile} \
            --reference {input.ref} --output {output.pymol_script}
        """

rule number_of_mutation_fig:
    input:
        treefile = "../enterovirus_d68/genome/results/tree_2018y.nwk",
        mutationfile = "../enterovirus_d68/genome/results/aa_muts_2018y.json",
    output:
        fig_name = "figures/number_of_changes_on_surface.pdf"
    shell:
        """
        python scripts/number_of_changes_on_surface.py --tree {input.treefile} \
            --mutations {input.mutationfile} --output {output.fig_name}
        """

rule epitope_changes:
    input:
        treefile = "../enterovirus_d68/vp1/results/tree_2018y.nwk",
        alignment = "../enterovirus_d68/vp1/results/aa_alignment_2018y_VP1.fasta",
        cladefile = "../enterovirus_d68/vp1/results/clades_2018y.json",
    output:
        epitope_changes = "results/epitope_changes.txt"
    shell:
        """
        python scripts/epitope_changes.py --tree {input.treefile} \
                --alignment {input.alignment} --clades {input.cladefile} \
                --output {output.epitope_changes}
        """


rule skyline:
     input:
         tree = "results/2011-9_genome_tree.nwk",
         aln = "../enterovirus_d68/genome/results/aligned_2018y.fasta",
         dates = "../enterovirus_d68/genome/results/metadata.tsv"
     output:
         skyline = "results/skyline/skyline.tsv"
     params:
         outdir = "results/skyline",
         npoints = 150
     shell:
         """
         treetime --tree {input.tree} --aln {input.aln} --dates {input.dates} --name-column strain \
 		 --coalescent skyline --n-skyline {params.npoints} --outdir {params.outdir}
         """

# make a new metadata which only includes 5 countries + 'rest of europe' + 'rest of world'
# for the traits migration estimation
# Also makes reduced regions of just 'china' 'europe' 'north_america' and 'rest of world'
rule meta_reduce_countries:
    input:
        meta = "../enterovirus_d68/vp1/results/metadata-ages.tsv"
    output:
        meta = "results/reduced_meta.tsv"
    shell:
        """
        Rscript scripts/reduce_metadata_countries.R {input.meta} {output.meta} 
        """

# Makes the reduced VP1 trees that include only tips from within the epidemic periods.
rule reduced_vp1_trees:
    input:
        treefile = "../enterovirus_d68/vp1/results/tree_2018y.nwk",
        branchfile = "../enterovirus_d68/vp1/results/branch_lengths_2018y.json",
        metafile = "../enterovirus_d68/vp1/results/metadata-ages.tsv"
    output:
        trees = ["results/2011-9_vp1_tree.nwk", "results/2018-9_vp1_tree.nwk",
                "results/2014-5_vp1_tree.nwk", "results/2016-7_vp1_tree.nwk"]
    shell:
        """
        python scripts/epidemic_only_trees.py \
            --tree {input.treefile} --branch-lengths {input.branchfile} \
            --meta {input.metafile} 
        """

# Note the trees used for this (produced by lineage_tracking_fig rule) have the tree branches
# in YEARS !!
# This uses VP1 trees and data with REDUCED COUNTRIES AND REGIONS
rule mugration:
    input:
        #ree = "results/{dset}_genome_tree.nwk",
        tree = rules.reduced_vp1_trees.output.trees,
        meta = rules.meta_reduce_countries.output.meta
    output:
        model = "results/mugration_{attr}_{dset}/GTR.txt"
    params:
        outdir = "results/mugration_{attr}_{dset}"
    shell:
        """
        treetime mugration --tree {input.tree} --states {input.meta} --attribute {wildcards.attr} --outdir {params.outdir} --missing-data nan
        """
