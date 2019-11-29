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
