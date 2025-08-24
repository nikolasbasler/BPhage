#!/bin/bash -l

mkdir -p $VSC_SCRATCH/BPhage/
cd $VSC_SCRATCH/BPhage/

wget https://zenodo.org/records/16937256/files/BPhage_test_dataset.tar.gz
tar -xzvf BPhage_test_dataset.tar.gz