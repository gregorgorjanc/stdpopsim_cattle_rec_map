#!/usr/bin/env Rscript

# setwd(dir = "~/Storages/GitBox/popsim/stdpopsim/stdpopsim_cattle_rec_map/Qanbari_and_Wittenburg_2020/")

# This script imports two files and merges them to create a RateMap object and
# plots the rate and map along each chromosome. The inital script. See the other
# script for the final version.

map <- read.csv(
  file = "12711_2020_593_MOESM2_ESM.csv",
  header = TRUE,
  sep = ","
)
head(map)
nrow(map) # 44696

map2 <- read.table(file = "physical_map.txt", header = FALSE)
head(map2)
nrow(map2) # 44696
colnames(map2) <- c("Chr", "Name", "Mbp_position", "Position")
head(map2)

test_order <- map$Name == map2$Name
sum(test_order) # 44696
sum(!test_order) # 0
# --> files are in the same order. Yuhhuu ...;)

map <- cbind(map, Position = map2$Position)
head(map)

# Doerte says: "Column recrate gives you the recombination rate between
# the marker named and its predecessor. (That’s why the first entry is zero.)
# “Mbp_inter_marker_distance” is the physical distance between the marker named
# and its predecessor. To achieve your target format, it should be fine if you
# take recrate_adjacent.... * 100 / Mbp_inter.... to get Rate (cM/Mb). Right?"

# Recalculating the intermarker distance to avoid loosing precision in using
# the rounded Mbp_inter_marker_distance
map$inter_marker_distance <- diff(c(0, map$Position))
head(map)

# Should we use cM_deterministic or cM_likelihood?
# Doert says that likelihood-based approach gives better results in case of huge
# data (as it was for Holstein) but in general that’s not true for smaller data.
map$recrate_adjacent_likelihood <- diff(c(0, map$cM_likelihood))
map$Rate_deterministic <- (map$recrate_adjacent_deterministic * 100) /
  (map$inter_marker_distance / 10^6)
map$Rate_likelihood <- (map$recrate_adjacent_likelihood * 100) /
  (map$inter_marker_distance / 10^6)
head(map)

plot(map$recrate_adjacent_likelihood ~ map$recrate_adjacent_deterministic)
plot(
  map$recrate_adjacent_likelihood ~ map$recrate_adjacent_deterministic,
  xlim = c(0.000001, max(map$recrate_adjacent_deterministic))
)
# TODO: I am getting negatives for recrate_adjacent_likelihood? Why!?

summary(map$recrate_adjacent_deterministic)
#      Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
# 0.0000000 0.0003732 0.0005016 0.0005465 0.0006697 0.0071024
summary(map$recrate_adjacent_likelihood)
#       Min.    1st Qu.     Median       Mean    3rd Qu.       Max.       NA's
# -115.18931    0.00000    0.00000    0.02258    0.07844    1.68733         80
# TODO: OK, something more fundamental is wrong here. What should I do?

RateMap <- map[, c("Chr", "Position", "Rate", "cM_deterministic")]
colnames(RateMap) <- c("Chromosome", "Position(bp)", "Rate(cM/Mb)", "Map(cM)")
chrs <- unique(RateMap$Chr)
for (chr in chrs) {
  # chr <- 1
  sel <- RateMap$Chr == chr
  tmp <- RateMap[sel, ]
  tmp <- tmp[order(tmp$Position), ]
  write.table(
    file = paste0("chr_", chr, "_ratemap.txt"),
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
