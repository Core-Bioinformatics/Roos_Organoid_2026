#!/bin/bash

# this calls a custom-made python script that follows the RNA velocity workflow
conda activate rna_velocity
python3 rna_velocity_dynamical.py \
    "output/data_for_python" \
    "output/rna_velocity" \
    "cellranger/looms" \
    50 \
    30 \
    10

