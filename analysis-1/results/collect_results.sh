#!/bin/bash

# Collect the figures for the manuscript from the results of the steps

THIS_SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
cd $THIS_SCRIPT_DIR
STEPSDIR=../steps

function getresult { local STEP=$1 SUB=$2 SRC=$3 FIG=$4 DST=$5
  SRCFN=$STEPSDIR/$STEP/results/$SUB/$SRC
  if [ -f $SRCFN ]; then
    mkdir -p $FIG/elements
    cp $SRCFN $FIG/elements/$DST
    echo "  $DST collected"
    return 0
  else
    echo "  $DST missing: $SRC not found in results of step $STEP"
    return 1
  fi
}

function getplot { local STEP=$1 SRC=$2; FIG=$3; DST=$4
  getresult $STEP plots $SRC $FIG $DST
  return $?
}

function fig_report { fig=$1 n_missing=$2
  if [ $n_missing -eq 0 ]; then
    echo -e "$fig completed"
  else
    echo -e "$fig missing $n_missing subplots"
  fi
  echo
}

function fig1 {
  fig=Fig1
  echo "=== $fig ==="
  n_missing=0
  getplot 005.0.V.T.splitted_CD8EM umap.by_cond.singleR_fine.mod.svg \
    $fig Fig1B.svg
  n_missing=$[n_missing+$?]
  getplot 004.0.V.T.integrated singleR_proportion_test.svg \
    $fig Fig1C.svg
  n_missing=$[n_missing+$?]
  getplot 005.0.V.T.splitted_CD8EM umap.by_cond.Terminal_effector_CD4_T_cells.svg \
    $fig Fig1D.svg
  n_missing=$[n_missing+$?]
  getplot 005.0.V.T.splitted_CD8EM umap.by_cond.MAIT_cells.svg \
    $fig Fig1E.svg
  n_missing=$[n_missing+$?]
  getplot 005.0.V.T.splitted_CD8EM umap.by_cond.Gamma_delta_T_cells.svg \
    $fig Fig1F.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function fig2 {
  fig=Fig2
  echo "=== $fig ==="
  n_missing=0
  getplot 006.0.V.T.metabolism_plots scMetabolism.Dotplot.T.svg \
    $fig Fig2A.svg
  n_missing=$[n_missing+$?]
  getplot 005.0.V.T.splitted_CD8EM phase_marker.MKI67.svg \
    $fig Fig2B.svg
  n_missing=$[n_missing+$?]
  getplot 005.0.V.T.splitted_CD8EM umap.by_cond.phase.svg \
    $fig Fig2C.svg
  n_missing=$[n_missing+$?]
  getplot 009.0.V.T.isolated_CD8EM isolated_CD8EM.umap_plot_by_subcluster.svg \
    $fig Fig2D.svg
  n_missing=$[n_missing+$?]
  getplot 009.0.V.T.isolated_CD8EM isolated_CD8EM.celltype_proportions.svg \
    $fig Fig2E.svg
  n_missing=$[n_missing+$?]
  getplot 009.0.V.T.isolated_CD8EM isolated_CD8EM.feature_plot_MAP3K8.svg \
    $fig Fig2F.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results selected_genes_heatmap/heatmap_prolif_CD8TEM.pdf \
    $fig Fig2G.pdf
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function fig4 {
  fig=Fig4
  echo "=== $fig ==="
  n_missing=0
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Naive_CD4_T_cells.svg \
    Fig4 Fig4A.1.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.T_regulatory_cells.svg \
    Fig4 Fig4A.2.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Th1_cells.svg \
    Fig4 Fig4A.3.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Th2_cells.svg \
    Fig4 Fig4A.4.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Th1_Th17_cells.svg \
    Fig4 Fig4A.5.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Th17_cells.svg \
    Fig4 Fig4A.6.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Follicular_helper_T_cells.svg \
    Fig4 Fig4A.7.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Naive_CD8_T_cells.svg \
    Fig4 Fig4A.8.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Effector_memory_CD8_T_cells.svg \
    Fig4 Fig4A.9.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Central_memory_CD8_T_cells.svg \
    Fig4 Fig4A.10.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Gamma-delta_T_cells.svg \
    Fig4 Fig4A.11.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results selected_genes_heatmap/heatmap_T_subtypes.pdf \
    Fig4 Fig4B.pdf
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function fig5 {
  fig=Fig5
  echo "=== $fig ==="
  n_missing=0
  getplot 010.0.V.T.TCR_results T.clonotype_size_umapharm.svg \
    Fig5 Fig5A.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.CD4__TEMRA_clonotype_size.svg \
    Fig5 Fig5B.1.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.Prolif__CD8__EM_clonotype_size.svg \
    Fig5 Fig5B.2.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.celltypes_clonotypes_barplots.svg \
    Fig5 Fig5C.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.alpha_beta_pairings.CD4__TEMRA.svg \
    Fig5 Fig5D.1.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.alpha_beta_pairings.CD8_EM.svg \
    Fig5 Fig5D.2.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.CD4__TEMRA.svg \
    Fig5 Fig5E.1.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.CD8__EM.svg \
    Fig5 Fig5E.2.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function fig6 {
  fig=Fig6
  echo "=== $fig ==="
  n_missing=0
  getplot 104.0.V.B.integrated umap.celltypes.svg \
    Fig6 Fig6A.svg
  n_missing=$[n_missing+$?]
  getplot 105.0.V.B.splitted_plasmablasts singleR_proportion_test.svg \
    Fig6 Fig6B.svg
  n_missing=$[n_missing+$?]
  getplot 105.0.V.B.splitted_plasmablasts umap.by_cond.clustering_res_0.4.svg \
    Fig6 Fig6C.svg
  n_missing=$[n_missing+$?]
  getplot 105.0.V.B.splitted_plasmablasts phase_marker.MKI67.svg \
    Fig6 Fig6D.svg
  n_missing=$[n_missing+$?]
  getplot 105.0.V.B.splitted_plasmablasts umap.by_cond.phase.svg \
    Fig6 Fig6E.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Naive_B_cells.svg \
    Fig6 Fig6F.1.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Non-switched_memory_B_cells.svg \
    Fig6 Fig6F.2.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Switched_memory_B_cells.svg \
    Fig6 Fig6F.3.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Plasmablasts.svg \
    Fig6 Fig6F.4.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Proliferating_plasmablasts.svg \
    Fig6 Fig6F.5.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Non-proliferating_plasmablasts.svg \
    Fig6 Fig6F.6.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Exhausted_B_cells.svg \
    Fig6 Fig6F.7.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function fig7 {
  fig=Fig7
  echo "=== $fig ==="
  n_missing=0
  # TODO: vertical order of the dots, can it be done automatically?
  getplot 110.0.V.B.BCR_results B.clonotype_size_umapharm.svg \
    Fig7 Fig7A.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results B.celltypes_clonotypes_barplots.svg \
    Fig7 Fig7B.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results isotype_umap.svg \
    Fig7 Fig7C.svg
  n_missing=$[n_missing+$?]
  # TODO: S left H right
  getplot 110.0.V.B.BCR_results heavy_chain_Cgene_group_Switched_memory_B_cells.svg \
    Fig7 Fig7D.1.svg
  n_missing=$[n_missing+$?]
  # TODO: S left H right
  getplot 110.0.V.B.BCR_results heavy_chain_Cgene_group_Plasmablasts.svg \
    Fig7 Fig7D.2.svg
  n_missing=$[n_missing+$?]
  # TODO: S left H right
  getplot 110.0.V.B.BCR_results heavy_chain_Cgene_group_Exhausted_B_cells.svg \
    Fig7 Fig7D.3.svg
  n_missing=$[n_missing+$?]
  getplot 111.0.V.B.enclone_plots Non_proliferating_plasmablasts.Hospitalized.svg \
    Fig7 Fig7E.1.svg
  n_missing=$[n_missing+$?]
  getplot 111.0.V.B.enclone_plots Proliferating_plasmablasts.Hospitalized.svg \
    Fig7 Fig7E.2.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results SHM_boxplot.H.Plasmablasts_pr_vs_nonpr.svg \
    Fig7 Fig7F.svg
  n_missing=$[n_missing+$?]
  # TODO: check differences in signif
  getplot 110.0.V.B.BCR_results B.VH_VL_pairings.HPb.svg \
    Fig7 Fig7G.svg
  n_missing=$[n_missing+$?]
  # TODO: PDF: fix
  getplot 110.0.V.B.BCR_results most_freq_CDR3s_aa.vdj.pdf \
    Fig7 Fig7H.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function figS4 {
  fig=FigS4
  echo "=== $fig ==="
  n_missing=0
  getplot 010.0.V.T.TCR_results T.Non_prolif__CD8__EM_clonotype_size.svg \
    FigS4 FigS4A.1.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.CD8__TEMRA_clonotype_size.svg \
    FigS4 FigS4A.2.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.GD_clonotype_size.svg \
    FigS4 FigS4A.3.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.alpha_beta_pairings.CD8__TEMRA.svg \
    FigS4 FigS4B.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function figS5 {
  fig=FigS5
  echo "=== $fig ==="
  n_missing=0
  # TODO: check diff in .1.
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.CD4__Naive.svg \
    FigS5 FigS5.1.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.Tfh.svg \
    FigS5 FigS5.2.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.Th1.svg \
    FigS5 FigS5.3.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.Th1_Th17.svg \
    FigS5 FigS5.4.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.Th17.svg \
    FigS5 FigS5.5.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.Th2.svg \
    FigS5 FigS5.6.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.Tregs.svg \
    FigS5 FigS5.7.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.CD8__Naive.svg \
    FigS5 FigS5.8.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.CD8__CM.svg \
    FigS5 FigS5.9.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.CD8__TEMRA.svg \
    FigS5 FigS5.10.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function figS6 {
  fig=FigS6
  echo "=== $fig ==="
  n_missing=0
  getplot 010.0.V.T.TCR_results VDJ_cdr3s_aa.overlap.svg \
    FigS6 FigS6.A.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results VJ_cdr3s_aa.overlap.svg \
    FigS6 FigS6.B.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results VDJ_VJ_cdr3s_aa.overlap.svg \
    FigS6 FigS6.C.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function figS7 {
  fig=FigS7
  echo "=== $fig ==="
  n_missing=0
  getplot 110.0.V.B.BCR_results SHM_boxplot.Non-switched_memory_B_cells.svg \
    FigS7 FigS7A.1.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results SHM_boxplot.Exhausted_B_cells.svg \
    FigS7 FigS7A.2.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results SHM_boxplot.Switched_memory_B_cells.svg \
    FigS7 FigS7A.3.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results B.VH_VL_pairings.HPb.svg \
    Fig7 Fig7G.svg
  n_missing=$[n_missing+$?]
  # TODO: missing
  #  FigS7 FigS7B.svg
  n_missing=$[n_missing+1]
  fig_report $fig $n_missing
}

echo "#### T cells analysis main figures ####"
echo
fig1
fig2
fig4
fig5
echo

echo "#### B cells analysis main figures ####"
echo
fig6
fig7
echo

echo "#### T cells analysis supplementary figures ####"
echo
figS4
figS5
figS6
echo

echo "#### B cells analysis supplementary figures ####"
figS7

