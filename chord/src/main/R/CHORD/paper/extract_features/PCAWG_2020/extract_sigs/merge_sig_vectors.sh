#!/bin/bash

guixr load-profile ~/.guix-profile --<<EOF
Rscript /hpc/cuppen/projects/P0013_WGS_patterns_Diagn/datasets/processed/PCAWG_2020/scripts/extract_sigs/merge_sig_vectors.R
EOF
