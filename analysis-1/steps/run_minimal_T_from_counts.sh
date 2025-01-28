#!/bin/bash
cd 001.0.P.T.create_seurat_obj/rundir && ./run_all_samples.sh && cd -
cd 002.0.P.T.qc_filter/rundir/ && ./run_all_samples.sh && cd -
knit2html 003.0.P.T.merging.seurat/rundir/merging.Rmd
knit2html 003.0.V.T.merged/rundir/post_merging.Rmd
knit2html 004.0.P.T.integration.harmony/rundir/harmony.Rmd
knit2html 004.0.P.T.integration.harmony/rundir/harmony_post_processing.Rmd
knit2html 004.0.V.T.integrated/rundir/post_harmony_analysis.Rmd
knit2html 005.0.P.T.split_CD8EM/rundir/split.Rmd
knit2html 005.0.V.T.splitted_CD8EM/rundir/post_split_analysis.Rmd
knit2html 007.0.P.T.compute_de/rundir/compute_de.Rmd
knit2html 007.0.V.T.de_results/rundir/de_dotplots.Rmd
knit2html 007.0.V.T.de_results/rundir/deseq_tables_and_vulcanos.Rmd
knit2html 007.0.V.T.de_results/rundir/selected_genes_heatmap.Rmd
knit2html 009.0.P.T.isolate_CD8EM/rundir/isolate_effmemCD8.Rmd
knit2html 009.0.V.T.isolated_CD8EM/rundir/analysis_isolated_CD8EM.Rmd &
knit2html 010.0.P.T.TCR_analysis.platypus/rundir/platypus_create.Rmd
knit2html 010.0.P.T.TCR_analysis.platypus/rundir/post_processing.Rmd
knit2html 010.0.P.T.TCR_analysis.platypus/rundir/further_post_processing.Rmd
knit2html 010.0.V.T.TCR_results/rundir/gene_usage.Rmd &
knit2html 010.0.V.T.TCR_results/rundir/cdr3_region.Rmd &
knit2html 010.0.V.T.TCR_results/rundir/clonotypes.Rmd &
wait
