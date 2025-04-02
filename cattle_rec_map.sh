
#!/bin/bash
# Process cattle recombination rate maps for stdpopsim

# Prepare by creating and entering cattle_rec_map directory
# git clone
# cd stdpopsim_cattle_rec_map

if [[ $(basename $(pwd)) != "cattle_rec_map" ]]; then
    echo "We must be in the cattle_rec_map directory"
    exit 1
fi

# Ma et al. (2015) Cattle Sex-Specific Recombination and Genetic Control from a Large Pedigree Analysis
# https://doi.org/10.1371/journal.pgen.1005387
mkdir -p Ma_et_al_2015
cd Ma_et_al_2015
# Tried wget or curl to download files from https://datadryad.org/dataset/doi:10.5061/dryad.q2q84
# but couldn't make it work:(
# Manually save cattle_rmap.txt, crossover_female.data, and crossover_male.data into this folder
# See cattle_rmap_to_RateMap_initial.R for exploration on how to remap this map from modifed
# UMD3.1 assembly to ARS-UCD1.2 assembly
# This script below actually does the remapping and avoids shrinking the map
Rscript ./cattle_rmap_to_RateMap.R
cd ..

# Shen et al. (2018) Characterization of recombination features and the genetic basis in multiple cattle breeds
# https://doi.org/10.1186/s12864-018-4705-y
mkdir -p Shen_et_al_2018
cd Shen_et_al_2018
wget https://static-content.springer.com/esm/art%3A10.1186%2Fs12864-018-4705-y/MediaObjects/12864_2018_4705_MOESM1_ESM.rmap
cd ..

# Schnabel (2018) ARS-UCD1.2 Cow Genome Assembly: Mapping of all existing variants
# https://www.animalgenome.org/repository/cattle/UMC_bovine_coordinates
# Ma et al. (2015) and Shen et al (2018) used a modification of the UMD3.1 assembly.
# Schnabel (2018) provides a map of all of the variants from existing bovine genotyping assays
# onto the ARS-UCD1.2 cow genome assembly
mkdir -p Schnabel_2018
cd Schnabel_2018
wget https://www.animalgenome.org/repository/cattle/UMC_bovine_coordinates/UMC_marker_names_180910.zip
wget https://www.animalgenome.org/repository/cattle/UMC_bovine_coordinates/UMC_marker_names_180910.zip.md5
md5sum UMC_marker_names_180910.zip
# b6d554793e9447fed41ed52ff0d1f529 UMC_marker_names_180910.zip
cat UMC_marker_names_180910.zip.md5
# b6d554793e9447fed41ed52ff0d1f529 UMC_marker_names_180910.zip
# check that the md5sum matches!
unzip UMC_marker_names_180910.zip


