#!/usr/bin/env Rscript

# setwd(dir = "~/Storages/GitBox/popsim/stdpopsim/stdpopsim_cattle_rec_map/Qanbari_and_Wittenburg_2020/")

# This script imports two files and merges them to create a RateMap object and
# plots the rate and map along each chromosome.

map <- read.csv(
  file = "GeneticMap_Holstein.csv",
  header = TRUE,
  sep = ","
)
head(map)
nrow(map) # 44696 before, but now we have 44631
sum(map$candidate_misplacement) # 0 because the missaligned markers were removed

# Ensure appropriate order!
map <- map[order(map$Chr, map$bp_position), ]
head(map)

# Handle NAs
summary(map$bp_position)
#  Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 52724  21741094  43610191  48902424  70098936 157503718
# --> all good here
summary(map$cM_deterministic)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 0.00   20.95   41.39   43.00   61.22  129.43
# --> all good here
summary(map$cM_likelihood)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's
#    0.00   25.15   44.90   46.21   64.21  126.95     210
# --> OK, no data for some markers to estimate recombination rate using likelihood
#     approach. What should we do - we will be calculating differences in these
#     values bellow, which will propagate these NAs!
sel <- is.na(map$cM_likelihood)
map[sel, ]
which(sel)[1]
map[c(263:265), ]
# I will interpolate ...
n <- nrow(map)
for (row in which(sel)) {
  # row <- which(sel)[1]
  if (row > 1 & row < n) {
    tmp <- map[c(row + c(-1, 0, 1)), ]
    if (all(tmp$Chr == tmp$Chr[1])) {
      NATest <- is.na(tmp$cM_likelihood[1]) | is.na(tmp$cM_likelihood[3])
      if (!NATest) {
        diffCM <- map$cM_likelihood[row + 1] - map$cM_likelihood[row - 1]
        diffBP <- map$bp_position[row + 1] - map$bp_position[row - 1]
        diffBP2 <- map$bp_position[row] - map$bp_position[row - 1]
        map$cM_likelihood[row] <- map$cM_likelihood[row - 1] +
          diffCM / diffBP * diffBP2
      } else {
        print("Marker")
        print(map[row, ])
        print("Rows")
        print(tmp)
        warning("NA's values in the neighborhood too!")
      }
    } else {
      if (tmp$Chr[1] == tmp$Chr[2]) {
        # I am at the end of a chromosome, so taking previous cM value
        map$cM_likelihood[row] <- map$cM_likelihood[row - 1]
      } else {
        print("Marker")
        print(map[row, ])
        print("Rows")
        print(tmp)
        stop("I need three rows within the chromosome to interpolate!")
      }
    }
  } else {
    stop("I need three rows to interpolate!")
  }
}
summary(map$cM_likelihood)
# Still have 21 NA's ...
# I will manually fix these ...
sel <- is.na(map$cM_likelihood)
map[sel, ]
# 1&2 easy one - constant value
s <- which(sel)[1]
map[s + -2:2, ]
map[c(7725, 7726), "cM_likelihood"] <- 20.47567
map[s + -2:2, ]
# 3&4
s <- which(sel)[3]
map[s + -2:2, ]
diffCM <- diff(map[s + c(-1, 2), "cM_likelihood"])
diffBP <- diff(map[s + c(-1, 2), "bp_position"])
diffBP2 <- diff(map[s + c(-1, 0), "bp_position"])
diffBP3 <- diff(map[s + c(-1, 1), "bp_position"])
map$cM_likelihood[s] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP2
map$cM_likelihood[s + 1] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP3
map[s + -2:2, ]
# 5&6
s <- which(sel)[5]
map[s + -2:2, ]
diffCM <- diff(map[s + c(-1, 2), "cM_likelihood"])
diffBP <- diff(map[s + c(-1, 2), "bp_position"])
diffBP2 <- diff(map[s + c(-1, 0), "bp_position"])
diffBP3 <- diff(map[s + c(-1, 1), "bp_position"])
map$cM_likelihood[s] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP2
map$cM_likelihood[s + 1] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP3
map[s + -2:2, ]
# 7&8
s <- which(sel)[7]
map[s + -2:2, ]
diffCM <- diff(map[s + c(-1, 2), "cM_likelihood"])
diffBP <- diff(map[s + c(-1, 2), "bp_position"])
diffBP2 <- diff(map[s + c(-1, 0), "bp_position"])
diffBP3 <- diff(map[s + c(-1, 1), "bp_position"])
map$cM_likelihood[s] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP2
map$cM_likelihood[s + 1] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP3
map[s + -2:2, ]
# 9&10
s <- which(sel)[9]
map[s + -2:2, ]
diffCM <- diff(map[s + c(-1, 2), "cM_likelihood"])
diffBP <- diff(map[s + c(-1, 2), "bp_position"])
diffBP2 <- diff(map[s + c(-1, 0), "bp_position"])
diffBP3 <- diff(map[s + c(-1, 1), "bp_position"])
map$cM_likelihood[s] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP2
map$cM_likelihood[s + 1] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP3
map[s + -2:2, ]
# 11&12
s <- which(sel)[11]
map[s + -2:2, ]
diffCM <- diff(map[s + c(-1, 2), "cM_likelihood"])
diffBP <- diff(map[s + c(-1, 2), "bp_position"])
diffBP2 <- diff(map[s + c(-1, 0), "bp_position"])
diffBP3 <- diff(map[s + c(-1, 1), "bp_position"])
map$cM_likelihood[s] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP2
map$cM_likelihood[s + 1] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP3
map[s + -2:2, ]
# 13&14
s <- which(sel)[13]
map[s + -2:2, ]
diffCM <- diff(map[s + c(-1, 2), "cM_likelihood"])
diffBP <- diff(map[s + c(-1, 2), "bp_position"])
diffBP2 <- diff(map[s + c(-1, 0), "bp_position"])
diffBP3 <- diff(map[s + c(-1, 1), "bp_position"])
map$cM_likelihood[s] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP2
map$cM_likelihood[s + 1] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP3
map[s + -2:2, ]
# 15&16&17
s <- which(sel)[15]
map[s + -2:3, ]
diffCM <- diff(map[s + c(-1, 3), "cM_likelihood"])
diffBP <- diff(map[s + c(-1, 3), "bp_position"])
diffBP2 <- diff(map[s + c(-1, 0), "bp_position"])
diffBP3 <- diff(map[s + c(-1, 1), "bp_position"])
diffBP4 <- diff(map[s + c(-1, 2), "bp_position"])
map$cM_likelihood[s] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP2
map$cM_likelihood[s + 1] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP3
map$cM_likelihood[s + 2] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP4
map[s + -2:3, ]
# 18&19
s <- which(sel)[18]
map[s + -2:2, ]
diffCM <- diff(map[s + c(-1, 2), "cM_likelihood"])
diffBP <- diff(map[s + c(-1, 2), "bp_position"])
diffBP2 <- diff(map[s + c(-1, 0), "bp_position"])
diffBP3 <- diff(map[s + c(-1, 1), "bp_position"])
map$cM_likelihood[s] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP2
map$cM_likelihood[s + 1] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP3
map[s + -2:2, ]
# 20&21
s <- which(sel)[20]
map[s + -2:2, ]
diffCM <- diff(map[s + c(-1, 2), "cM_likelihood"])
diffBP <- diff(map[s + c(-1, 2), "bp_position"])
diffBP2 <- diff(map[s + c(-1, 0), "bp_position"])
diffBP3 <- diff(map[s + c(-1, 1), "bp_position"])
map$cM_likelihood[s] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP2
map$cM_likelihood[s + 1] <- map$cM_likelihood[s - 1] +
  diffCM / diffBP * diffBP3
map[s + -2:2, ]

# Doerte says: "Column recrate gives you the recombination rate between
# the marker named and its predecessor. (That’s why the first entry is zero.)
# “Mbp_inter_marker_distance” is the physical distance between the marker named
# and its predecessor. To achieve your target format, it should be fine if you
# take recrate_adjacent.... * 100 / Mbp_inter.... to get Rate (cM/Mb). Right?"

# recrate_adjacent.... * 100 is to convert Morgan (recombination rate) to cM

# Recalculating the intermarker distance to avoid loosing precision in using
# the rounded Mbp_inter_marker_distance

# Note that this data is in format of "rate & distance between marker and it's PREVIOUS marker" to
# So, now we have this:
# head(map)
#   Chr                            Name bp_position cM_deterministic    recrate
# 1   1              ARS-BFGL-NGS-16466      907810       0.00000000 0.00000000
# 2   1             ARS-BFGL-NGS-105096      993060       0.00091345 0.00000913
# 3   1 Hapmap34944-BES1_Contig627_1906     1032564       0.00318736 0.00002274
# Then correct "distance between marker and it's PREVIOUS marker" is:
# 907810-0=907810
# 993060-907810=85250
# 1032564-993060=39504

chrs <- unique(map$Chr)
for (chr in chrs) {
  sel <- map$Chr == chr
  map$inter_marker_distance[sel] <- diff(c(0, map$bp_position[sel]))
}
head(map)
summary(map$inter_marker_distance)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    0   28403   39999   55616   64142 2_750_827
# We will divide by inter_marker_distance so let's check where we have 0s
sel <- map$inter_marker_distance == 0
sum(sel) # 29
map[sel, ]
for (row in which(sel)) {
  # row <- which(sel)[1]
  cat("\n")
  print(paste0("Row: ", row))
  print(map[row + -2:2, ])
}
# TODO: Ask Doerte what to do here!?
#       This will impact the below calculation of:
#       * Rate_deterministic
#       * Rate_likelihood

# Should we use cM_deterministic or cM_likelihood?
# Doert says that likelihood-based approach gives better results in case of huge
# data (as it was for Holstein) but in general that’s not true for smaller data.
# Trying to calculate recrate_adjacent_likelihood
chrs <- unique(map$Chr)
map$recrate_adjacent_likelihood <- NA
for (chr in chrs) {
  # chr <- 1
  sel <- map$Chr == chr
  map$recrate_adjacent_likelihood[sel] <- diff(
    c(0, map$cM_likelihood[sel] / 100)
    # / 100 to convert to Morgans (to keep the below equations the same)
  )
  map$Rate_deterministic <- (map$recrate_adjacent_deterministic * 100) /
    (map$inter_marker_distance / 10^6)
  map$Rate_likelihood <- (map$recrate_adjacent_likelihood * 100) /
    (map$inter_marker_distance / 10^6)
  # * 100 to convert M to cM
  # / 10^6 to get cM / Mbp
}
head(map)
head(map, n = 20)

with(map, plot(recrate_adjacent_likelihood ~ recrate_adjacent_deterministic))
dev.copy(png, file = "recrate_adjacent_likelihood_vs_deterministic.png")
dev.off()
with(
  map,
  plot(recrate_adjacent_likelihood ~ recrate_adjacent_deterministic, log = "xy")
)
dev.copy(png, file = "recrate_adjacent_likelihood_vs_deterministic_log.png")
dev.off()
# TODO: this looks odd! - this should be effectively a straight line (more or less)

summary(map$recrate_adjacent_likelihood)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -1.700e-08  0.000e+00  0.000e+00  5.678e-04  7.834e-04  1.531e-02
# TODO: where are the negative values coming from?!
summary(map$recrate_adjacent_deterministic)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0000000 0.0003704 0.0004976 0.0005398 0.0006640 0.0027658
selL <- is.finite(map$recrate_adjacent_likelihood)
selD <- is.finite(map$recrate_adjacent_deterministic)
summary(map$recrate_adjacent_likelihood[selL])
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.
# -1.700e-08  0.000e+00  0.000e+00  5.678e-04  7.834e-04  1.531e-02
# TODO: where are the negative values coming from?!
summary(map$recrate_adjacent_deterministic[selD])
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0000000 0.0003704 0.0004976 0.0005398 0.0006640 0.0027658

summary(map$Rate_deterministic)
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#  0.0000  0.6983  1.1621     Inf  1.8207     Inf
# TODO: deal with Infs - by handling inter_marker_distance == 0
summary(map$Rate_likelihood)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.      NA's
# -0.000072  0.000000  0.000000       Inf  1.517572       Inf        18
# TODO: deal with Infs - by handling inter_marker_distance == 0
# TODO: NA's?
selL <- is.finite(map$Rate_deterministic)
selD <- is.finite(map$Rate_likelihood)
summary(map$Rate_deterministic[selD])
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max.
# 0.00     0.70     1.16     3.09     1.82 66768.00
summary(map$Rate_likelihood[selL])
#  Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.000     0.000     0.000     1.676     1.516 17822.717

with(map, plot(Rate_likelihood ~ Rate_deterministic))
dev.copy(png, file = "rate_likelihood_vs_deterministic.png")
dev.off()
with(map, plot(Rate_likelihood ~ Rate_deterministic, log = "xy"))
dev.copy(png, file = "rate_likelihood_vs_deterministic_log.png")
# TODO: this looks odd! - this should be effectively a straight line (more or less)

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
# the physical position in that line (inclusive) and the physical position on the
# next line (exclusive)."
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
#   Chr                            Name bp_position cM_deterministic    recrate bp_distance
# 1   1              ARS-BFGL-NGS-16466      907810       0.00000000 0.00000000 907810 ...
# 2   1             ARS-BFGL-NGS-105096      993060       0.00091345 0.00000913  85250 ...
# 3   1 Hapmap34944-BES1_Contig627_1906     1032564       0.00318736 0.00002274  39504 ...
# The first row above is from 0-907810, meaning we should change to
#   Chr                            Name bp_position cM_deterministic    recrate bp_distance
# 0   1                      Chr1_START           0       0.00000000 0.00000913 907810 ...
# 1   1              ARS-BFGL-NGS-16466      907810       0.00091345 0.00002274  85250 ...
# 2   1             ARS-BFGL-NGS-105096      993060       0.00318736        ...  39504 ...
# 3   1 Hapmap34944-BES1_Contig627_1906     1032564       ...

cols <- c(
  "Chr",
  "Name",
  "bp_position",
  "inter_marker_distance",
  "cM_deterministic",
  "cM_likelihood",
  "recrate_adjacent_deterministic",
  "recrate_adjacent_likelihood",
  "Rate_deterministic",
  "Rate_likelihood"
)
map <- map[, cols]
chrs <- unique(map$Chr)
for (chr in chrs) {
  # chr <- 1
  cat("\n")
  print(paste0("Chromosome: ", chr))
  sel <- map$Chr == chr
  tmp <- map[sel, ]
  tmp0 <- tmp[1, , drop = FALSE]
  tmp0$Name <- paste0("Chr", chr, "_START")
  tmp0$bp_position <- 0
  tmp0$inter_marker_distance <- 0
  tmp0$cM_deterministic <- 0
  tmp0$cM_likelihood <- 0
  tmp0$recrate_adjacent_deterministic <- 0
  tmp0$recrate_adjacent_likelihood <- 0
  tmp0$Rate_deterministic <- 0
  tmp0$Rate_likelihood <- 0
  tmp <- rbind(tmp0, tmp)
  n <- nrow(tmp)
  cols <- c(
    "inter_marker_distance",
    "recrate_adjacent_deterministic",
    "recrate_adjacent_likelihood",
    "Rate_deterministic",
    "Rate_likelihood"
  )
  tmp[1:(n - 1), cols] <- tmp[2:n, cols]
  print(head(tmp))
  # tail(tmp)
  cols <- c(
    "recrate_adjacent_deterministic",
    "recrate_adjacent_likelihood",
    "Rate_deterministic",
    "Rate_likelihood"
  )
  tmp[n, cols] <- 0
  print(tail(tmp))
  cat("\n")
  if (chr == 1) {
    map2 <- tmp
  } else {
    map2 <- rbind(map2, tmp)
  }
}

for (type in c("deterministic", "likelihood")) {
  # type <- "deterministic"
  dir.create(path = type)
  RateMap <- map2[, c("Chr", "bp_position", paste0(c("Rate_", "cM_"), type))]
  colnames(RateMap) <- c("Chromosome", "Position(bp)", "Rate(cM/Mb)", "Map(cM)")
  chrs <- unique(RateMap$Chr)
  for (chr in chrs) {
    # chr <- 1
    sel <- RateMap$Chr == chr
    tmp <- RateMap[sel, ]
    tmp <- tmp[order(tmp$Position), ]
    write.table(
      file = paste0(type, "/chr_", chr, "_ratemap.txt"),
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
    dev.copy(png, file = paste0(type, "/chr_", chr, "_rate.png"))
    dev.off()
    plot(
      tmp$`Map(cM)` ~ tmp$`Position(bp)`,
      type = "l",
      xlab = "Position (bp)",
      ylab = "Map (cM)",
      main = paste0("Chromosome ", chr)
    )
    dev.copy(png, file = paste0(type, "/chr_", chr, "_map.png"))
    dev.off()
  }
}
