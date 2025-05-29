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
    $fig ${fig}B.svg
  n_missing=$[n_missing+$?]
  getplot 004.0.V.T.integrated singleR_proportion_test.svg \
    $fig ${fig}C.svg
  n_missing=$[n_missing+$?]
  getplot 005.0.V.T.splitted_CD8EM umap.by_cond.Terminal_effector_CD4_T_cells.svg \
    $fig ${fig}D.svg
  n_missing=$[n_missing+$?]
  getplot 005.0.V.T.splitted_CD8EM umap.by_cond.MAIT_cells.svg \
    $fig ${fig}E.svg
  n_missing=$[n_missing+$?]
  getplot 005.0.V.T.splitted_CD8EM umap.by_cond.Gamma_delta_T_cells.svg \
    $fig ${fig}F.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function fig2 {
  fig=Fig2
  echo "=== $fig ==="
  n_missing=0
  getplot 006.0.V.T.metabolism_plots scMetabolism.Dotplot.T.svg \
    $fig ${fig}A.svg
  n_missing=$[n_missing+$?]
  getplot 005.0.V.T.splitted_CD8EM phase_marker.MKI67.svg \
    $fig ${fig}B.svg
  n_missing=$[n_missing+$?]
  getplot 005.0.V.T.splitted_CD8EM umap.by_cond.phase.svg \
    $fig ${fig}C.svg
  n_missing=$[n_missing+$?]
  getplot 009.0.V.T.isolated_CD8EM isolated_CD8EM.umap_plot_by_subcluster.svg \
    $fig ${fig}D.svg
  n_missing=$[n_missing+$?]
  getplot 009.0.V.T.isolated_CD8EM isolated_CD8EM.celltype_proportions.svg \
    $fig ${fig}E.svg
  n_missing=$[n_missing+$?]
  getplot 009.0.V.T.isolated_CD8EM isolated_CD8EM.feature_plot_MAP3K8.svg \
    $fig ${fig}F.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function fig4 {
  fig=Fig4
  echo "=== $fig ==="
  n_missing=0
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Naive_CD4_T_cells.svg \
    ${fig} ${fig}A.1.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.T_regulatory_cells.svg \
    ${fig} ${fig}A.2.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Th1_cells.svg \
    ${fig} ${fig}A.3.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Th2_cells.svg \
    ${fig} ${fig}A.4.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Th1_Th17_cells.svg \
    ${fig} ${fig}A.5.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Th17_cells.svg \
    ${fig} ${fig}A.6.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Follicular_helper_T_cells.svg \
    ${fig} ${fig}A.7.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Naive_CD8_T_cells.svg \
    ${fig} ${fig}A.8.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Effector_memory_CD8_T_cells.svg \
    ${fig} ${fig}A.9.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Central_memory_CD8_T_cells.svg \
    ${fig} ${fig}A.10.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results dotplots/T.de_genes_dotplot.Gamma-delta_T_cells.svg \
    ${fig} ${fig}A.11.svg
  n_missing=$[n_missing+$?]
  getplot 007.0.V.T.de_results selected_genes_heatmap/heatmap_T_subtypes.pdf \
    ${fig} ${fig}B.pdf
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function fig5 {
  fig=Fig5
  echo "=== $fig ==="
  n_missing=0
  getplot 010.0.V.T.TCR_results T.clonotype_size_umapharm.svg \
    ${fig} ${fig}A.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.CD4__TEMRA_clonotype_size.svg \
    ${fig} ${fig}B.1.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.Prolif__CD8__EM_clonotype_size.svg \
    ${fig} ${fig}B.2.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.celltypes_clonotypes_barplots.svg \
    ${fig} ${fig}C.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.alpha_beta_pairings.CD4__TEMRA.svg \
    ${fig} ${fig}D.1.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.alpha_beta_pairings.CD8_EM.svg \
    ${fig} ${fig}D.2.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.CD4__TEMRA.svg \
    ${fig} ${fig}E.1.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.CD8__EM.svg \
    ${fig} ${fig}E.2.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function fig6 {
  fig=Fig6
  echo "=== $fig ==="
  n_missing=0
  getplot 104.0.V.B.integrated umap.celltypes.svg \
    ${fig} ${fig}A.svg
  n_missing=$[n_missing+$?]
  getplot 105.0.V.B.splitted_plasmablasts singleR_proportion_test.pdf \
    ${fig} ${fig}B.pdf
  n_missing=$[n_missing+$?]
  getplot 105.0.V.B.splitted_plasmablasts umap.by_cond.clustering_res_0.4.svg \
    ${fig} ${fig}C.svg
  n_missing=$[n_missing+$?]
  getplot 105.0.V.B.splitted_plasmablasts phase_marker.MKI67.svg \
    ${fig} ${fig}D.svg
  n_missing=$[n_missing+$?]
  getplot 105.0.V.B.splitted_plasmablasts umap.by_cond.phase.svg \
    ${fig} ${fig}E.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Naive_B_cells.svg \
    ${fig} ${fig}F.1.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Non-switched_memory_B_cells.svg \
    ${fig} ${fig}F.2.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Switched_memory_B_cells.svg \
    ${fig} ${fig}F.3.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Proliferating_plasmablasts.svg \
    ${fig} ${fig}F.4.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Non-proliferating_plasmablasts.svg \
    ${fig} ${fig}F.5.svg
  n_missing=$[n_missing+$?]
  getplot 107.0.V.B.de_results dotplots/B.de_genes_dotplot.Exhausted_B_cells.svg \
    ${fig} ${fig}F.6.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function fig7 {
  fig=Fig7
  echo "=== $fig ==="
  n_missing=0
  getplot 110.0.V.B.BCR_results B.clonotype_size_umapharm.svg \
    ${fig} ${fig}A.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results B.celltypes_clonotypes_barplots.svg \
    ${fig} ${fig}B.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results isotype_umap.svg \
    ${fig} ${fig}C.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results heavy_chain_Cgene_group_Switched_memory_B_cells.svg \
    ${fig} ${fig}D.1.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results heavy_chain_Cgene_group_Plasmablasts.svg \
    ${fig} ${fig}D.2.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results heavy_chain_Cgene_group_Exhausted_B_cells.svg \
    ${fig} ${fig}D.3.svg
  n_missing=$[n_missing+$?]
  getplot 111.0.V.B.enclone_plots Non_proliferating_plasmablasts.Hospitalized.svg \
    ${fig} ${fig}E.1.svg
  n_missing=$[n_missing+$?]
  getplot 111.0.V.B.enclone_plots Proliferating_plasmablasts.Hospitalized.svg \
    ${fig} ${fig}E.2.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results SHM_boxplot.H.Plasmablasts_pr_vs_nonpr.svg \
    ${fig} ${fig}F.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results B.VH_VL_pairings.HPb.svg \
    ${fig} ${fig}G.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results most_freq_CDR3s_aa.vdj.svg \
    ${fig} ${fig}H.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function figS5 {
  fig=FigS5
  echo "=== $fig ==="
  n_missing=0
  getplot 005.0.V.T.splitted_CD8EM umap.by_sample.singleR_fine.mod.svg \
    $fig ${fig}.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function figS6 {
  fig=FigS6
  echo "=== $fig ==="
  n_missing=0
  getplot 010.0.V.T.TCR_results T.Non_prolif__CD8__EM_clonotype_size.svg \
    ${fig} ${fig}A.1.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.CD8__TEMRA_clonotype_size.svg \
    ${fig} ${fig}A.2.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.GD_clonotype_size.svg \
    ${fig} ${fig}A.3.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.alpha_beta_pairings.CD8__TEMRA.svg \
    ${fig} ${fig}B.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function figS7 {
  fig=FigS7
  echo "=== $fig ==="
  n_missing=0
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.CD4__Naive.svg \
    ${fig} ${fig}.1.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.Tfh.svg \
    ${fig} ${fig}.2.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.Th1.svg \
    ${fig} ${fig}.3.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.Th1_Th17.svg \
    ${fig} ${fig}.4.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.Th17.svg \
    ${fig} ${fig}.5.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.Th2.svg \
    ${fig} ${fig}.6.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.Tregs.svg \
    ${fig} ${fig}.7.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.CD8__Naive.svg \
    ${fig} ${fig}.8.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.CD8__CM.svg \
    ${fig} ${fig}.9.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results T.most_freq_CDR3aa.VDJ.by_sample.CD8__TEMRA.svg \
    ${fig} ${fig}.10.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function figS8 {
  fig=FigS8
  echo "=== $fig ==="
  n_missing=0
  getplot 010.0.V.T.TCR_results VDJ_cdr3s_aa.overlap.svg \
    ${fig} ${fig}.A.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results VJ_cdr3s_aa.overlap.svg \
    ${fig} ${fig}.B.svg
  n_missing=$[n_missing+$?]
  getplot 010.0.V.T.TCR_results VDJ_VJ_cdr3s_aa.overlap.svg \
    ${fig} ${fig}.C.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function figS9 {
  fig=FigS9
  echo "=== $fig ==="
  n_missing=0
  getplot 105.0.V.B.splitted_plasmablasts umap.by_sample.singleR_fine.mod.svg \
    $fig ${fig}.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

function figS10 {
  fig=FigS10
  echo "=== $fig ==="
  n_missing=0
  getplot 110.0.V.B.BCR_results SHM_boxplot.Non-switched_memory_B_cells.svg \
    $fig ${fig}A.1.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results SHM_boxplot.Exhausted_B_cells.svg \
    $fig ${fig}A.2.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results SHM_boxplot.Switched_memory_B_cells.svg \
    $fig ${fig}A.3.svg
  n_missing=$[n_missing+$?]
  getplot 110.0.V.B.BCR_results B.VH_VL_pairings.Non_switched_memory_B_cells.svg \
    $fig ${fig}B.svg
  n_missing=$[n_missing+$?]
  fig_report $fig $n_missing
}

if [ "$1" != "B" ]; then
  echo "#### T cells analysis main figures ####"
  echo
  fig1
  fig2
  fig4
  fig5
  echo
fi

if [ "$1" != "T" ]; then
  echo "#### B cells analysis main figures ####"
  echo
  fig6
  fig7
  echo
fi

if [ "$1" != "B" ]; then
  echo "#### T cells analysis supplementary figures ####"
  echo
  figS5
  figS6
  figS7
  figS8
  echo
fi

if [ "$1" != "T" ]; then
  echo "#### B cells analysis supplementary figures ####"
  figS9
  figS10
fi
