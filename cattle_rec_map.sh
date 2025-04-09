
#!/bin/bash
# Process cattle recombination rate maps for stdpopsim

# This script is intended to make it easier to download and process
# the recombination rate maps, but in almost all!? cases there are
# some manual steps involved, so this script is mostly a convenience
# documentation of the process!

# Prepare and enter the repo
# git clone https://github.com/gregorgorjanc/stdpopsim_cattle_rec_map.git
# cd stdpopsim_cattle_rec_map

if [[ $(basename $(pwd)) != "stdpopsim_cattle_rec_map" ]]; then
    echo "We must be in the stdpopsim_cattle_rec_map directory"
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

# Qanbari and Wittenburg (2020) Male recombination map of the autosomal genome in German Holstein
# https://doi.org/10.1186/s12711-020-00593-z
mkdir -p Qanbari_and_Wittenburg_2020
cd Qanbari_and_Wittenburg_2020
wget https://static-content.springer.com/esm/art%3A10.1186%2Fs12711-020-00593-z/MediaObjects/12711_2020_593_MOESM2_ESM.xlsx
# Manually convert 12711_2020_593_MOESM2_ESM.xlsx to 12711_2020_593_MOESM2_ESM.csv
# Hmm, this file has some rounded-up physical position.
# I have contacted Doerte Wittenburg, who sent me physical_map.txt file and
# the updated map (with some markers removed due to likely assembly errors) GeneticMap_Holstein.xlsx
# (I converted this to csv). This last file is from the work published in
# Melzer et al. (2023) CLARITY: a Shiny app for interactive visualisation of the bovine physical-genetic map
# https://doi.org/10.3389/fgene.2023.1082782
Rscript ./map_to_RateMap.R
cd ..

# Brekke et al. (2023) Variation and genetic control of individual recombination rates in Norwegian Red dairy cattle
# https://doi.org/10.3168/jds.2022-22368
mkdir -p Brekke_et_al_2023
cd Brekke_et_al_2023
# Tried wget or curl to download files from https://figshare.com/articles/dataset/Autosomal_linkage_map_NRF/20976067/1?file=37270660
# but couldn't make it work:(
# Manually save all_chr_linkage_map_JDS.xls into this folder
# Manually convert all_chr_linkage_map_JDS.xls to all_chr_linkage_map_JDS.csv
# Contacted Cathrine who send a more detailed NRF_LinkageMap.txt
cd ..
