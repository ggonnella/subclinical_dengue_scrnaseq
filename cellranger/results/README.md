This directory shall contain the outs directory of Cellranger
(version 6 or 7) renamed to for T cells ``sample_1_T``, ``sample_2_T``, ...
and for B cells ``sample_1_B``, ``sample_2_B``, ...

The analysis pipeline will use the following files/subdirectories where
``$SAMPLE`` is the sample name (e.g. ``sample_1_T``) and
``$VDJ`` is either ``vdj_b`` or ``vdj_t``:
1. sample h5 counts file (used to create Seurat objects) named
   ``$SAMPLE/per_sample_outs/$SAMPLE/count/sample_feature_bc_matrix.h5`` (Cellranger 6) or
   ``$SAMPLE/per_sample_outs/$SAMPLE/count/sample_filtered_feature_bc_matrix.h5`` (Cellranger 7)
2. filtered contig annotations CSV file (used to filter by VDJ availability) named
   ``$SAMPLE/per_sample_outs/$SAMPLE/$VDJ/filtered_contig_annotations.csv``
3. all contig annotations JSON file (used to prepare honeycomb plots with enclone) named
   ``$SAMPLE/outs/multi/$VDJ/all_contigs_annotations.json``
4. vdj directory (used for creating Platypus vgm objects) named
   ``$SAMPLE/per_samples_outs/$SAMPLE/$VDJ``
