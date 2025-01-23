#!/bin/bash

# Collect the figures for the manuscript from the results of the steps

STEPSDIR=../steps

function fig1 {
  mkdir -p Fig1/elements
  cp $STEPSDIR/008.0.V.T.splitted_CD8EM/results/plots/umap.by_cond.singleR_fine.mod.svg \
    Fig1/elements/Fig1B.svg
  cp ../steps/... Fig1/elements/Fig1C.svg # MUST BE ADDED TO INTEGRATED
  cp ../steps/... Fig1/elements/Fig1D.svg # MUST BE OUTPUT AS PLOT FROM SPLITTED
  cp ../steps/... Fig1/elements/Fig1E.svg # MUST BE OUTPUT AS PLOT FROM SPLITTED
  cp ../steps/... Fig1/elements/Fig1F.svg # MUST BE OUTPUT AS PLOT FROM SPLITTED
}

function fig2 {
  mkdir -p Fig2/elements
  cp $STEPSDIR/009.0.V.T.metabolism_plots/results/plots/scMetabolism.Dotplot.T.svg \
    Fig2/elements/Fig2A.svg
  cp $STEPSDIR/008.0.V.T.splitted_CD8EM/results/plots/phase_marker.MKI67.svg \
    Fig2/elements/Fig2B.svg
  cp $STEPSDIR/008.0.V.T.splitted_CD8EM/results/plots/umap.by_cond.phase.svg \
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
  # TODO: rename to T.... (step 013)
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/B.clonotype_size_umapharm.svg \
    Fig5/elements/Fig5A.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.CD4__TEMRA_clonotype_size.svg \
    Fig5/elements/Fig5B.1.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.Prolif__CD8__EM_clonotype_size.svg \
    Fig5/elements/Fig5B.2.svg
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.celltypes_clonotypes_barplots.svg \
    Fig5/elements/Fig5C.svg
  cp $STEPSDIR/ \
    Fig5/elements/Fig5D.1.svg
  cp $STEPSDIR/ \
    Fig5/elements/Fig5D.2.svg
  # TODO: fix colors (step 013)
  cp $STEPSDIR/013.0.V.T.TCR_results/results/plots/T.most_freq_clonotypes.VDJ.by_sample.CD4__TEMRA.svg \
    Fig5/elements/Fig5E.1.svg
  # TODO: save plot for unsplitted EM and fix colors (step 013)
  cp $STEPSDIR/ \
    Fig5/elements/Fig5E.2.svg
}

function fig6 {
  mkdir -p Fig6/elements
  # TODO: save plot from step 107
  cp $STEPDIR/ \
    Fig6/elements/Fig6A.svg
  # TODO: save plot from step 108
  cp $STEPDIR/ \
    Fig6/elements/Fig6B.svg
  # TODO: save plot from step 108
  cp $STEPDIR/ \
    Fig6/elements/Fig6C.svg
  # TODO: save plot from step 108
  cp $STEPDIR/ \
    Fig6/elements/Fig6D.svg
  # TODO: save plot from step 108
  cp $STEPDIR/ \
    Fig6/elements/Fig6E.svg
  cp $STEPDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Naive_B_cells.svg \
    Fig6/elements/Fig6F.1.svg
  cp $STEPDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Non-switched_memory_B_cells.svg \
    Fig6/elements/Fig6F.2.svg
  cp $STEPDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Switched_memory_B_cells.svg \
    Fig6/elements/Fig6F.3.svg
  cp $STEPDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Plasmablasts.svg \
    Fig6/elements/Fig6F.4.svg
  # TODO: missing!
  cp $STEPDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Proliferating_plasmablasts.svg \
    Fig6/elements/Fig6F.5.svg
  cp $STEPDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Non-proliferating_plasmablasts.svg \
    Fig6/elements/Fig6F.6.svg
  cp $STEPDIR/110.0.V.B.de_results/results/plots/dotplots/B.de_genes_dotplot.Exhausted_B_cells.svg \
    Fig6/elements/Fig6F.7.svg
}

function fig7 {
  mkdir -p Fig7/elements
  # TODO: vertical order of the dots, can it be done automatically?
  cp $STEPDIR/113.0.V.B.BCR_results/results/plots/B.clonotype_size_umapharm.svg \
    Fig7/elements/Fig7A.svg
  cp $STEPDIR/113.0.V.B.BCR_results/results/plots/B.celltypes_clonotypes_barplots.svg \
    Fig7/elements/Fig7B.svg
  # TODO: save from step 113
  cp $STEPDIR/ \
    Fig7/elements/Fig7C.svg
  # TODO: save from step 113
  cp $STEPDIR \
    Fig7/elements/Fig7D.1.svg
  # TODO: save from step 113
  cp $STEPDIR \
    Fig7/elements/Fig7D.2.svg
  # TODO: save from step 113
  cp $STEPDIR \
    Fig7/elements/Fig7D.3.svg
  cp $STEPDIR \
    Fig7/elements/Fig7E.1.svg
  cp $STEPDIR \
    Fig7/elements/Fig7E.2.svg
  cp $STEPDIR \
    Fig7/elements/Fig7F.1.svg
  cp $STEPDIR \
    Fig7/elements/Fig7F.2.svg
  cp $STEPDIR \
    Fig7/elements/Fig7G.svg
  cp $STEPDIR \
    Fig7/elements/Fig7H.1.svg
  cp $STEPDIR \
    Fig7/elements/Fig7H.2.svg
}

function figS4 {
  mkdir -p FigS4/elements
  # source step: 113
    FigS4/elements/FigS4A.1.svg
    FigS4/elements/FigS4A.2.svg
    FigS4/elements/FigS4A.3.svg
    FigS4/elements/FigS4B.svg
}

function figS5 {
  mkdir -p FigS5/elements
  # source step: 113
    FigS5/elements/FigS5.1.svg
    FigS5/elements/FigS5.2.svg
    FigS5/elements/FigS5.3.svg
    FigS5/elements/FigS5.4.svg
    FigS5/elements/FigS5.5.svg
    FigS5/elements/FigS5.6.svg
    FigS5/elements/FigS5.7.svg
    FigS5/elements/FigS5.8.svg
    FigS5/elements/FigS5.9.svg
    FigS5/elements/FigS5.10.svg
}

function figS6 {
  mkdir -p FigS6/elements
  # source step: 113
    FigS6/elements/FigS6.1.svg
    FigS6/elements/FigS6.2.svg
    FigS6/elements/FigS6.3.svg
}

function figS7 {
  mkdir -p FigS7/elements
  # source step: 113
    FigS7/elements/FigS7A.1.svg
    FigS7/elements/FigS7A.2.svg
    FigS7/elements/FigS7A.3.svg
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

