#!/usr/bin/env Rscript

# setwd(dir = "~/Storages/GitBox/popsim/stdpopsim/stdpopsim_cattle_rec_map/Qanbari_and_Wittenburg_2020/")

# This script imports two files and merges them to create a RateMap object and
# plots the rate and map along each chromosome.

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
map$Rate <- (map$recrate_adjacent_deterministic * 100) /
  (map$inter_marker_distance / 10^6)
head(map)

RateMap <- map[, c("Chr", "Position", "Rate", "cM_deterministic")]
colnames(RateMap) <- c("Chromosome", "Position(bp)", "Rate(cM/Mb)", "Map(cM)")
chrs <- unique(RateMap$Chr)
for (chr in chrs) {
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
