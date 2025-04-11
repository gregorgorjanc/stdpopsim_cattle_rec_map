#!/usr/bin/env Rscript

# setwd(dir = "~/Storages/GitBox/popsim/stdpopsim/stdpopsim_cattle_rec_map/Brekke_et_al_2023/")

# This script imports two files and merges them to create a HapMap object and
# plots the rate and map along each chromosome.

map <- read.table(
  file = "NRF_LinkageMap.txt",
  header = TRUE
)
head(map)
nrow(map) # 35880

# Ensure appropriate order!
map <- map[order(map$chr, map$bp), ]
head(map)

summary(map$bp)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 53134  21830515  43360008  48632905  70180764 157892202
# --> all good here
summary(map$cM)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00   25.04   45.01   46.05   63.95  126.08
# --> all good here

chrs <- unique(map$chr)
tot <- c(0, 0, 0)
for (chr in chrs) {
  # chr <- 1
  sel <- map$chr == chr
  tot <- tot +
    c(max(map$male_cM[sel]), max(map$female_cM[sel]), max(map$cM[sel]))
}
tot
# 2493.052 2308.798 2400.925

map$bp_diff <- NA
chrs <- unique(map$chr)
for (chr in chrs) {
  sel <- map$chr == chr
  map$bp_diff[sel] <- diff(c(0, map$bp[sel]))
}
head(map)
summary(map$bp_diff)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#     23   35113   52419   69096   84798 2549027
# --> all good here

chrs <- unique(map$chr)
for (chr in chrs) {
  sel <- map$chr == chr
  print(head(map[sel, ]))
  print(tail(map[sel, ]))
}
# --> all good here

chrs <- unique(map$chr)
map$cM_diff <- NA
for (chr in chrs) {
  # chr <- 1
  sel <- map$chr == chr
  map$cM_diff[sel] <- diff(c(0, map$cM[sel]))
}
map$recRate <- (map$cM_diff) / (map$bp_diff / 10^6)
# / 10^6 to get cM / Mbp
head(map)
head(map, n = 20)

summary(map$recRate)
#     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#  0.0000    0.3722    0.8789    1.6167    1.8165 1347.8261
# --> some high rates, but visual inspection of the abovge plots suggests this
#    is largely OK, maybe that 1347 is an outlier ... but keeping it in to
#    follow the publication

# Change format from "distance between marker and it's PREVIOUS marker" to
# "distance between marker and it's NEXT marker" to satisfy HapMap/RateMap format
#
# Following https://tskit.dev/msprime/docs/stable/api.html#msprime.RateMap.read_hapmap
# we have this example:
#
# Chromosome  Position(bp)  Rate(cM/Mb)  Map(cM)
# chr10       48232         0.1614       0.002664
# chr10       48486         0.1589       0.002705
# chr10       50009         0.159        0.002947
#
# and these key texts:
#
# "The value in the rate column in a given line gives the constant rate between
# the physical position in THAT line (inclusive) and the physical position on the
# NEXT line (exclusive)."
#
# "... the first row has a nonzero genetic map position (last column, cM),
# implying a nonzero recombination rate before that position, that is assumed to
# extend to the start of the chromosome (at position 0 bp). However, if the first
# line has a nonzero bp position (second column) and a zero genetic map position
# (last column, cM), then the recombination rate before that position is unknown,
# producing missing data.
#
# Chromosome  Position(bp)  Rate(cM/Mb)  Map(cM)
# chr10       0             x            0
# chr10       48232         y=0.1614     0.002664
# chr10       48486         z=0.1589     0.002705
# chr10       50009         0.159        0.002947
# x = (0.002664 - 0) / ((48232 - 0) / 10^6) = 0.05523304
# y = (0.002705 - 0.002664) / ((48486 - 48232) / 10^6) = 0.1614
# z = (0.002947 - 0.002705) / ((50009 - 48486) / 10^6) = 0.1589

# So, now we have this:
# chr                           snpid      bp     cM male_cM female_cM bp_diff cM_diff   recRate
# 1   1              BovineHD0100000015  716721 0.0000   0.000     0.000  716721  0.0000 0.0000000
# 2   1          Hapmap43437-BTA-101873  776231 0.0185   0.037     0.000   59510  0.0185 0.3108721
# 3   1              ARS-BFGL-NGS-16466  907810 0.3125   0.293     0.332  131579  0.2940 2.2343991
# 4   1              BovineHD0100000096  987646 0.3750   0.305     0.445   79836  0.0625 0.7828549
# 5   1 Hapmap34944-BES1_Contig627_1906 1032564 0.3855   0.326     0.445   44918  0.0105 0.2337593
# 6   1              ARS-BFGL-NGS-98142 1110393 0.4560   0.386     0.526   77829  0.0705 0.9058320
# The first row above is from 0-716721, meaning we should change to
#   chr                           snpid      bp     cM male_cM female_cM bp_diff cM_diff   recRate
# 1   1                      Chr1_START       0 0.0000   0.000     0.000  716721  0.0000 0.0000000
# 2   1              BovineHD0100000015  716721 0.0185   0.037     0.000   59510  0.0185 0.3108721
# 3   1          Hapmap43437-BTA-101873  776231 0.3125   0.293     0.332  131579  0.2940 2.2343991
# 4   1              ARS-BFGL-NGS-16466  907810 0.3750   0.305     0.445   79836  0.0625 0.7828549
# 5   1              BovineHD0100000096  987646 0.3855   0.326     0.445   44918  0.0105 0.2337593
# 6   1 Hapmap34944-BES1_Contig627_1906 1032564 0.4560   0.386     0.526   77829  ...
# 7   1              ARS-BFGL-NGS-98142 1110393 ...

chrs <- unique(map$chr)
for (chr in chrs) {
  # chr <- 1
  cat("\n")
  print(paste0("Chromosome: ", chr))
  sel <- map$chr == chr
  tmp <- map[sel, ]
  tmp0 <- tmp[1, , drop = FALSE]
  tmp0$snpid <- paste0("Chr", chr, "_START")
  tmp0$bp <- 0
  tmp0$bp_diff <- 0
  tmp0$cM <- 0
  tmp0$cM_diff <- 0
  tmp0$male_cM <- 0
  tmp0$female_cM <- 0
  tmp0$recRate <- 0
  tmp <- rbind(tmp0, tmp)
  n <- nrow(tmp)
  cols <- c("bp_diff", "cM_diff", "recRate")
  tmp[1:(n - 1), cols] <- tmp[2:n, cols]
  print(head(tmp))
  # tail(tmp)
  tmp[n, "recRate"] <- 0
  print(tail(tmp))
  cat("\n")
  if (chr == 1) {
    map2 <- tmp
  } else {
    map2 <- rbind(map2, tmp)
  }
}

HapMap <- map2[, c("chr", "bp", "recRate", "cM")]
colnames(HapMap) <- c("Chromosome", "Position(bp)", "Rate(cM/Mb)", "Map(cM)")
chrs <- unique(HapMap$Chromosome)
head(HapMap)
for (chr in chrs) {
  # chr <- 1
  sel <- HapMap$Chromosome == chr
  tmp <- HapMap[sel, ]
  tmp <- tmp[order(tmp$`Position(bp)`), ]
  # head(tmp)
  write.table(
    file = paste0("chr_", chr, "_HapMap.txt"),
    x = tmp,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    sep = " "
  )
  plot(
    tmp$`Rate(cM/Mb)` ~ tmp$`Position(bp)`,
    type = "l",
    xlab = "Position (bp)",
    ylab = "Rate (cM/Mb)",
    main = paste0("Chromosome ", chr)
  )
  dev.copy(png, file = paste0("chr_", chr, "_rate.png"))
  dev.off()
  plot(
    tmp$`Map(cM)` ~ tmp$`Position(bp)`,
    type = "l",
    xlab = "Position (bp)",
    ylab = "Map (cM)",
    main = paste0("Chromosome ", chr)
  )
  dev.copy(png, file = paste0("chr_", chr, "_map.png"))
  dev.off()
}
