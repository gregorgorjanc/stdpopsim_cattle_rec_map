# stdpopsim_cattle_rec_map
Converting published cattle recombination rate (genetic) maps for [stdpopsim](https://github.com/popsim-consortium/stdpopsim).

[Here](https://tskit.dev/msprime/docs/stable/api.html#msprime.RateMap.read_hapmap)
is the description of the HapMap file that we are aiming to produce and
[here](https://popsim-consortium.github.io/stdpopsim-docs/stable/development.html#adding-a-recombination-genetic-map-or-annotation)
is the description of the process for adding such a file to stdpopsim.

[Here](https://github.com/popsim-consortium/stdpopsim/issues/602)
is the GitHub issue about adding genetic map(s) for cattle to stdpopsim

Several publications estimated genetic maps for cattle (let us know of the publications we missed!):
  * [Ma et al. (2015) Cattle Sex-Specific Recombination and Genetic Control from a Large Pedigree Analysis](https://doi.org/10.1371/journal.pgen.1005387) (for Holstein, 186,927 three-generation families, 59,309 autosomal SNPs, 8.5+ million maternal and paternal recombination events, giving 25.5 Morgans for males and 23.2 Morgans for females; without the X chromosome) - the challenge here is that this one is based on the modified UMD3.1 assembly from 2009 (see Note1 below) - considered in this repo, but assembly changes made porting this map hard/impossible --> ABANDONED (for now, until the author reruns estimation against the new assembly)
  * [Shen et al. (2018) Characterization of recombination features and the genetic basis in multiple cattle breeds](https://doi.org/10.1186/s12864-018-4705-y) (for Holstein - taking Ma et al. (2015) result, Jersey, Brown Swiss, and Ayrshire, 161,309 three-generation families, 58,982 autosomal SNPs) - the challenge here is that this one is based on the modified UMD3.1 assembly from 2009 (see Note1 below) - considered in this repo, but assembly changes made porting this map hard/impossible --> ABANDONED (for now, until the author reruns estimation against the new assembly)
  * [Qanbari and Wittenburg (2020) Male recombination map of the autosomal genome in German Holstein](https://doi.org/10.1186/s12711-020-00593-z) (for male Holstein, 44,696 autosomal SNP, 876 half-sib families with 39-4236 progeny (in total over 367K progeny), 8.9 million paternal recombinations, giving 24.43 Morgans for 2486 Mbp and a refined estimate of 25.35 Morgans) - improved in [Melzer et al. (2023) CLARITY: a Shiny app for interactive visualisation of the bovine physical-genetic map](https://doi.org/10.3389/fgene.2023.1082782) - considered in this repo, did some additional processing, but this is fundamentally a male genetic map --> ABANDONED (in favour of Wittenburg et al. (202?))
  * [Zhang et al. (2020) Genetic architecture of individual variation in recombination rate on the X chromosome in cattle](https://doi.org/10.1038/s41437-020-0341-9) --> TODO https://github.com/gregorgorjanc/stdpopsim_cattle_rec_map/issues/3
  * [Brekke et al. (2023) Variation and genetic control of individual recombination rates in Norwegian Red dairy cattle](https://doi.org/10.3168/jds.2022-22368) (for Norwegian Red, 35,880 autosomal SNP for 110,555 individuals, giving 2,492.9 cM for males and 2,308.9 cM for females) --> DONE
  * Wittenburg et al. (202?) A Resource of Bovine Genetic Maps --> TODO https://github.com/gregorgorjanc/stdpopsim_cattle_rec_map/issues/1

The script `cattle_rec_map.sh` tries to download various files (some have to be downloaded manually and for some I got updates from the authors!) into publication-specific folders, which also contain conversion scripts.

Note1: Publication of ARS-UCD1.2 assembly is [Rosen et al. (2020) De novo assembly of the cattle reference genome with single-molecule sequencing](https://academic.oup.com/gigascience/article/9/3/giaa021/5810242) and [Schnabel (2018) ARS-UCD1.2 Cow Genome Assembly: Mapping of all existing variants](https://www.animalgenome.org/repository/cattle/UMC_bovine_coordinates) provides mapping of various SNP array genotypes based on the UMD3.1 assembly to this newer assembly. I tried this with the Ma et al. (2015) map (see `Ma_et_al_2015/cattle_rmap_to_RateMap.R`) following discussion at [here](https://github.com/popsim-consortium/stdpopsim/issues/602#issuecomment-2767373686) and [here](https://github.com/popsim-consortium/stdpopsim/issues/602#issuecomment-2767826026), but got stuck on what to do with some larger blocks of chromosomes. Here is the GitHub [issue](https://github.com/popsim-consortium/stdpopsim/issues/571) and [PR](https://github.com/popsim-consortium/stdpopsim/pull/574) on lifting over recombination maps within stdpopsim.







