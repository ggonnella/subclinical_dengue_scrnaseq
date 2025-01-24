#!/bin/bash

# Collect the figures for the manuscript from the results of the steps

THIS_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $THIS_SCRIPT_DIR

STEPSDIR=../steps

function fig1 {
  mkdir -p Fig1/elements
  cp $STEPSDIR/005.0.V.T.splitted_CD8EM/results/plots/umap.by_cond.singleR_fine.mod.svg \
    Fig1/elements/Fig1B.svg
  cp $STEPSDIR/004.0.V.T.integrated/results/plots/singleR_proportion_test.svg \
    Fig1/elements/Fig1C.svg
  cp $STEPSDIR/005.0.V.T.splitted_CD8EM/results/plots/umap.by_cond.Terminal_effector_CD4_T_cells.svg \
    Fig1/elements/Fig1D.svg
  cp $STEPSDIR/005.0.V.T.splitted_CD8EM/results/plots/umap.by_cond.MAIT_cells.svg \
    Fig1/elements/Fig1E.svg
  cp $STEPSDIR/005.0.V.T.splitted_CD8EM/results/plots/umap.by_cond.Gamma_delta_T_cells.svg \
    Fig1/elements/Fig1F.svg
}

function fig2 {
  mkdir -p Fig2/elements
  cp $STEPSDIR/009.0.V.T.metabolism_plots/results/plots/scMetabolism.Dotplot.T.svg \
    Fig2/elements/Fig2A.svg
  cp $STEPSDIR/005.0.V.T.splitted_CD8EM/results/plots/phase_marker.MKI67.svg \
    Fig2/elements/Fig2B.svg
  cp $STEPSDIR/005.0.V.T.splitted_CD8EM/results/plots/umap.by_cond.phase.svg \
    Fig2/elements/Fig2C.svg
  cp $STEPSDIR/012.0.V.T.isolated_CD8EM/results/plots/isolated_CD8EM.umap_plot_by_subcluster.svg \
    Fig2/elements/Fig2D.svg
  cp $STEPSDIR/012.0.V.T.isolated_CD8EM/results/plots/isolated_CD8EM.celltype_proportions.svg \
    Fig2/elements/Fig2E.svg
  cp $STEPSDIR/012.0.V.T.isolated_CD8EM/results/plots/isolated_CD8EM.feature_plot_MAP3K8.svg \
    Fig2/elements/Fig2F.svg
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/selected_genes_heatmap/heatmap_prolif_CD8TEM.pdf \
    Fig2/elements/Fig2G.pdf
}

function fig4 {
  mkdir -p Fig4/elements
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/dotplots/T.de_genes_dotplot.Naive_CD4_T_cells.svg \
    Fig4/elements/Fig4A.1.svg
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/dotplots/T.de_genes_dotplot.T_regulatory_cells.svg \
    Fig4/elements/Fig4A.2.svg
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/dotplots/T.de_genes_dotplot.Th1_cells.svg \
    Fig4/elements/Fig4A.3.svg
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/dotplots/T.de_genes_dotplot.Th2_cells.svg \
    Fig4/elements/Fig4A.4.svg
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/dotplots/T.de_genes_dotplot.Th1_Th17_cells.svg \
    Fig4/elements/Fig4A.5.svg
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/dotplots/T.de_genes_dotplot.Th17_cells.svg \
    Fig4/elements/Fig4A.6.svg
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/dotplots/T.de_genes_dotplot.Follicular_helper_T_cells.svg \
    Fig4/elements/Fig4A.7.svg
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/dotplots/T.de_genes_dotplot.Naive_CD8_T_cells.svg \
    Fig4/elements/Fig4A.8.svg
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/dotplots/T.de_genes.dotplot.Effector_memory_CD8_T_cells.svg \
    Fig4/elements/Fig4A.9.svg
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/dotplots/T.de_genes.dotplot.Central_memory_CD8_T_cells.svg \
    Fig4/elements/Fig4A.10.svg
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/dotplots/T.de_genes_dotplot.Gamma-delta_T_cells.svg \
    Fig4/elements/Fig4A.11.svg
  cp $STEPSDIR/010.0.V.T.de_results/results/plots/selected_genes_heatmap/heatmap_T_subtypes.pdf \
    Fig4/elements/Fig4B.pdf
}

function fig5 {
  mkdir -p Fig5/elements
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.clonotype_size_umapharm.svg \
    Fig5/elements/Fig5A.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.CD4__TEMRA_clonotype_size.svg \
    Fig5/elements/Fig5B.1.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.Prolif__CD8__EM_clonotype_size.svg \
    Fig5/elements/Fig5B.2.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.celltypes_clonotypes_barplots.svg \
    Fig5/elements/Fig5C.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.alpha_beta_pairings.CD4__TEMRA.svg \
    Fig5/elements/Fig5D.1.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.alpha_beta_pairings.CD8_EM.svg \
    Fig5/elements/Fig5D.2.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_clonotypes.VDJ.by_sample.CD4__TEMRA.svg \
    Fig5/elements/Fig5E.1.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_clonotypes.VDJ.by_sample.CD8__EM.svg \
    Fig5/elements/Fig5E.2.svg
}

function fig6 {
  mkdir -p Fig6/elements
  cp $STEPSDIR/107.0.V.B.integrated/results/plots/umap.celltypes.svg \
    Fig6/elements/Fig6A.svg
  cp $STEPSDIR/108.0.V.B.splitted_plasmablasts/results/plots/singleR_proportion_test.svg \
    Fig6/elements/Fig6B.svg
  cp $STEPSDIR/108.0.V.B.splitted_plasmablasts/results/plots/umap.by_cond.clustering_res_0.4.svg \
    Fig6/elements/Fig6C.svg
  cp $STEPSDIR/108.0.V.B.splitted_plasmablasts/results/plots/phase_marker.MKI67.svg \
    Fig6/elements/Fig6D.svg
  cp $STEPSDIR/108.0.V.B.splitted_plasmablasts/results/plots/umap.by_cond.phase.svg \
    Fig6/elements/Fig6E.svg
  cp $STEPSDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Naive_B_cells.svg \
    Fig6/elements/Fig6F.1.svg
  cp $STEPSDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Non-switched_memory_B_cells.svg \
    Fig6/elements/Fig6F.2.svg
  cp $STEPSDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Switched_memory_B_cells.svg \
    Fig6/elements/Fig6F.3.svg
  cp $STEPSDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Plasmablasts.svg \
    Fig6/elements/Fig6F.4.svg
  cp $STEPSDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Proliferating_plasmablasts.svg \
    Fig6/elements/Fig6F.5.svg
  cp $STEPSDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Non-proliferating_plasmablasts.svg \
    Fig6/elements/Fig6F.6.svg
  cp $STEPSDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Exhausted_B_cells.svg \
    Fig6/elements/Fig6F.7.svg
}

function fig7 {
  mkdir -p Fig7/elements
  # TODO: vertical order of the dots, can it be done automatically?
  cp $STEPSDIR/113.0.V.B.BCR_results/results/plots/B.clonotype_size_umapharm.svg \
    Fig7/elements/Fig7A.svg
  cp $STEPSDIR/113.0.V.B.BCR_results/results/plots/B.celltypes_clonotypes_barplots.svg \
    Fig7/elements/Fig7B.svg
  cp $STEPSDIR/113.0.V.B.BCR_results/results/plots/isotype_umap.svg \
    Fig7/elements/Fig7C.svg
  # TODO: S left H right
  cp $STEPSDIR/113.0.V.B.BCR_results/results/plots/heavy_chain_Cgene_group_Switched_memory_B_cells.svg \
    Fig7/elements/Fig7D.1.svg
  # TODO: S left H right
  cp $STEPSDIR/113.0.V.B.BCR_results/results/plots/heavy_chain_Cgene_group_Plasmablasts.svg \
    Fig7/elements/Fig7D.2.svg
  # TODO: S left H right
  cp $STEPSDIR/113.0.V.B.BCR_results/results/plots/heavy_chain_Cgene_group_Exhausted_B_cells.svg \
    Fig7/elements/Fig7D.3.svg
  cp $STEPSDIR/114.0.V.B.enclone_plots/results/plots/Non_proliferating_plasmablasts.Hospitalized.svg \
    Fig7/elements/Fig7E.1.svg
  cp $STEPSDIR/114.0.V.B.enclone_plots/results/plots/Proliferating_plasmablasts.Hospitalized.svg \
    Fig7/elements/Fig7E.2.svg
  cp $STEPSDIR/113.0.V.B.BCR_results/results/plots/SHM_boxplot.H.Plasmablasts_pr_vs_nonpr.svg \
    Fig7/elements/Fig7F.svg
  # TODO: check differences in signif
  cp $STEPSDIR/113.0.V.B.BCR_results/results/plots/B.VH_VL_pairings.HPb.svg \
    Fig7/elements/Fig7G.svg
  # TODO: PDF: fix
  cp $STEPSDIR/113.0.V.B.BCR_results/results/plots/most_freq_CDR3s_aa.vdj.pdf \
    Fig7/elements/Fig7H.svg
}

function figS4 {
  mkdir -p FigS4/elements
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.Non_prolif__CD8__EM_clonotype_size.svg \
    FigS4/elements/FigS4A.1.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.CD8__TEMRA_clonotype_size.svg \
    FigS4/elements/FigS4A.2.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.GD_clonotype_size.svg \
    FigS4/elements/FigS4A.3.svg
  cd $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.alpha_beta_pairings.CD8__TEMRA.svg \
    FigS4/elements/FigS4B.svg
}

function figS5 {
  mkdir -p FigS5/elements
  # TODO: check diff in .1.
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_CDR3aa.VDJ.by_sample.CD4__Naive.svg \
    FigS5/elements/FigS5.1.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_CDR3aa.VDJ.by_sample.Tfh.svg \
    FigS5/elements/FigS5.2.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_CDR3aa.VDJ.by_sample.Th1.svg \
    FigS5/elements/FigS5.3.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_CDR3aa.VDJ.by_sample.Th1_Th17.svg \
    FigS5/elements/FigS5.4.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_CDR3aa.VDJ.by_sample.Th17.svg \
    FigS5/elements/FigS5.5.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_CDR3aa.VDJ.by_sample.Th2.svg \
    FigS5/elements/FigS5.6.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_CDR3aa.VDJ.by_sample.Tregs.svg \
    FigS5/elements/FigS5.7.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_CDR3aa.VDJ.by_sample.CD8__Naive.svg \
    FigS5/elements/FigS5.8.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_CDR3aa.VDJ.by_sample.CD8__CM.svg \
    FigS5/elements/FigS5.9.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_CDR3aa.VDJ.by_sample.CD8__TEMRA.svg \
    FigS5/elements/FigS5.10.svg
}

function figS6 {
  mkdir -p FigS6/elements
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/VDJ_cdr3s_aa.overlap.svg \
    FigS6/elements/FigS6.A.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/VJ_cdr3s_aa.overlap.svg \
    FigS6/elements/FigS6.B.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/VDJ_VJ_cdr3s_aa.overlap.svg \
    FigS6/elements/FigS6.C.svg
}

function figS7 {
  mkdir -p FigS7/elements
  # source step: 113
  cp $STEPSDIR/113.0.V.B.BCR_results/results/plots/SHM_boxplot.Non-switched_memory_B_cells.svg \
    FigS7/elements/FigS7A.1.svg
  cp $STEPSDIR/113.0.V.B.BCR_results/results/plots/SHM_boxplot.Exhausted_B_cells.svg \
    FigS7/elements/FigS7A.2.svg
  cp $STEPSDIR/113.0.V.B.BCR_results/results/plots/SHM_boxplot.Switched_memory_B_cells.svg \
    FigS7/elements/FigS7A.3.svg
  # TODO: missing
    FigS7/elements/FigS7B.svg
}

fig1
fig2
fig4
fig5
fig6
fig7
figS4
figS5
figS6
figS7

