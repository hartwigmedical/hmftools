#!/bin/bash

base_dir=/hpc/cog_bioinf/cuppen/project_data/Luan_projects/datasets/HMF_DR010_DR047/
R_SCRIPT=$base_dir/scripts/extract_sigs/merge_sig_vectors.R

guixr load-profile ~/.guix-profile --<<EOF
Rscript $R_SCRIPT
EOF