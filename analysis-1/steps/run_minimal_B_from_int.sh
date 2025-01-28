#!/bin/bash
knit2html 104.0.P.B.integration.harmony/rundir/harmony_post_processing.Rmd input_integrated_rdata=TRUE
knit2html 104.0.V.B.integrated/rundir/post_harmony_analysis.Rmd
knit2html 105.0.P.B.split_plasmablasts/rundir/split.Rmd
knit2html 105.0.V.B.splitted_plasmablasts/rundir/post_split_analysis.Rmd &
knit2html 107.0.P.B.de_analysis/rundir/compute_de.Rmd
knit2html 107.0.V.B.de_results/rundir/de_dotplots.Rmd &
knit2html 107.0.V.B.de_results/rundir/deseq_tables_and_vulcanos.Rmd &
knit2html 109.0.P.B.isolate_plasmablasts/rundir/isolate_plasmablasts.Rmd
knit2html 109.0.V.B.isolated_plasmablasts/rundir/analysis_isolated_plasmablasts.Rmd &
knit2html 110.0.P.B.BCR_analysis.platypus/rundir/platypus_create.Rmd
knit2html 110.0.P.B.BCR_analysis.platypus/rundir/post_processing.Rmd
knit2html 110.0.P.B.BCR_analysis.platypus/rundir/further_post_processing.Rmd
knit2html 110.0.V.B.BCR_results/rundir/gene_usage.Rmd &
knit2html 110.0.V.B.BCR_results/rundir/cdr3_region.Rmd &
knit2html 110.0.V.B.BCR_results/rundir/clonotypes.Rmd &
knit2html 110.0.V.B.BCR_results/rundir/shm_analysis.Rmd &
knit2html 111.0.P.B.enclone_preparation/rundir/output_barcode_table.Rmd
knit2html 111.0.V.B.enclone_plots/rundir/create_enclone_plots.Rmd
wait
