#!/usr/bin/env Rscript

# setwd(dir = "~/Storages/GitBox/popsim/stdpopsim/cattle_rec_map/Ma_et_al_2015/")

# This script imports the recombination map from Ma et al. (2015), which
# was based on a modified UMD3.1 assembly, and remaps it to the new ARS-UCD1.2.
# As explored in the cattle_rmap_to_RateMap_initial.R, remapping markers leads
# to map shrinkage. Here, we take the recombination rate from the lost markers
# add add it to the previous markers to avoid map shrinkage.

# ---- Import and inspect the map ----

map <- read.table(file = "cattle_rmap.txt", header = TRUE)
nrow(map) # 59309
head(map)
map <- map[order(map$Chr, map$Location), ]
head(map)
tail(map)
summary(map$map_f)
#       Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
# 0.000e+00 2.000e-09 1.280e-04 3.909e-04 4.870e-04 1.316e-02
summary(map$map_m)
#       Min.   1st Qu.    Median      Mean   3rd Qu.     Max.
# 0.000e+00 3.000e-09 1.371e-04 4.304e-04 5.226e-04 2.122e-02
map$map <- (map$map_f + map$map_m) / 2
summary(map$map)
hist(map$map)
# looks like the map column is actually recombination rate between the loci
colnames(map) <- c("Name", "Chr", "Position", "Rate_f", "Rate_m", "Rate")
sum(map$Rate_f)
# 23.18364
sum(map$Rate_m)
# 25.52622
sum(map$Rate)
# 24.35493

summary(map$Position)
# Min.   1st Qu.    Median      Mean   3rd Qu.      Max.
#    0  21393283  43799317  49343799  71485655 158319652

# Let's check this 0 position
sel <- map$Position == 0
sum(sel) # 1
map[sel, ]
head(map[map$Chr == 6, ])
#                        Name Chr Position        Rate_f Rate_m          Rate
# 14821 Hapmap38362-BTA-94562   6        0 9.382746e-105      0 4.691373e-105
# 14822    ARS-BFGL-NGS-87800   6    66096 9.382746e-105      0 4.691373e-105
# Odd - likely the position is not 0, but some missing value?

# Find which are previous marker positions and marker names
map <- map[order(map$Chr, map$Position), ]
chrs <- sort(unique(map$Chr))
map$NamePrevious <- NA
map$PositionPrevious <- NA
map$Rate_fPrevious <- NA
map$Rate_mPrevious <- NA
map$RatePrevious <- NA
map$NameNext <- NA
map$PositionNext <- NA
map$Rate_fNext <- NA
map$Rate_mNext <- NA
map$RateNext <- NA
for (chr in chrs) {
  # chr <- 1
  sel <- map$Chr == chr
  n <- length(map$Position[sel])
  map$NamePrevious[sel] <- c(NA, map$Name[sel][-n])
  map$PositionPrevious[sel] <- c(NA, map$Position[sel][-n])
  map$Rate_fPrevious[sel] <- c(NA, map$Rate_f[sel][-n])
  map$Rate_mPrevious[sel] <- c(NA, map$Rate_m[sel][-n])
  map$RatePrevious[sel] <- c(NA, map$Rate[sel][-n])

  map$NameNext[sel] <- c(map$Name[sel][-1], NA)
  map$PositionNext[sel] <- c(map$Position[sel][-1], NA)
  map$Rate_fNext[sel] <- c(map$Rate_f[sel][-1], NA)
  map$Rate_mNext[sel] <- c(map$Rate_m[sel][-1], NA)
  map$RateNext[sel] <- c(map$Rate[sel][-1], NA)
}

# ---- Import the new assembly and remap the map to it ----
# (ARS-UCD1.2 assembly, the above map is based on the modified UMD3.1 assembly)

assembly <- read.csv(
  file = "../Schnabel_2018/UMC_marker_names_180910/9913_ARS1.2_raw_data_180910.csv",
  header = TRUE
)

# First let's see if any markers miss altogether
test_presence <- map$Name %in% assembly$marker_name
sum(test_presence) # 59305
sum(!test_presence) # 4
map[!test_presence, ] # hmm, probably names have changed ...

move_rec_rate_to_previous_marker <- function(x, markers_to_remove) {
  # We assume we have t-1 (previous) t (current) and t+1 (next) markers and for
  # marker t we will move recombinantion rate to marker t-1. If t-1 is missing,
  # we will ???
  # x is the map data frame
  # markers_to_remove is a vector of marker names to remove
  removed <- NULL
  initialRateSum <- sum(x$Rate, na.rm = TRUE)
  x <- x[order(x$Chr, x$Position), ]
  sel_remove_markers <- x$Name %in% markers_to_remove
  no_previous_marker <- is.na(x$NamePrevious[sel_remove_markers])
  if (any(no_previous_marker)) {
    print(paste(
      sum(no_previous_marker),
      "marker(s) don't have a previous marker on the chromosome!"
    ))
    print(map[sel_remove_markers, ][no_previous_marker, ])
    print("We will drop their recombination rate!")
    remove_markers <- map[sel_remove_markers, ][no_previous_marker, "Name"]
    sel_remove_markers2 <- x$Name %in% remove_markers
    removed <- x[sel_remove_markers2, , drop = FALSE]
    x <- x[!sel_remove_markers2, , drop = FALSE]
    sel_remove_markers <- x$Name %in% markers_to_remove
  }
  # Working one marker at a time and from the end because there is s recursion -
  # a marker might be removed, but the same can happen to a marker before it!
  # This will be slow!!!!
  remove_markers <- x$Name[sel_remove_markers]
  sel_previous_markers <- x$Name %in% x$NamePrevious[sel_remove_markers]
  # x[sel_remove_markers, c("Name", "NamePrevious", "Rate", "RatePrevious", "RateNext")]
  # x[sel_previous_markers, c("Name", "Rate", "RatePrevious", "RateNext")]
  na2zero <- function(x) {
    x[is.na(x)] <- 0
    return(x)
  }
  count <- 0
  for (marker in rev(remove_markers)) {
    # rev() because we will process marker removal and moving rec rate from
    # the end of chromosome towards the beginning
    # marker <- rev(remove_markers)[25]
    # marker <- rev(remove_markers)[26]
    count <- count + 1
    sel <- which(x$NameNext == marker) # find the marker before the one we will remove
    if (length(sel) == 0) {
      # We have hit the first marker on chromosome!
      stop("Have to decide what to do with the first marker on chromosome!")
    }
    sel2 <- which(x$Name == marker)
    # x[sel, ]
    # x[sel2, ]
    x$Rate_f[sel] <- na2zero(x$Rate_f[sel]) + na2zero(x$Rate_f[sel2])
    x$Rate_m[sel] <- na2zero(x$Rate_m[sel]) + na2zero(x$Rate_m[sel2])
    x$Rate[sel] <- na2zero(x$Rate[sel]) + na2zero(x$Rate[sel2])
    if (is.null(removed)) {
      removed <- x[sel2, , drop = FALSE]
    } else {
      removed <- rbind(x[sel2, , drop = FALSE], removed)
    }
    x[sel2, "Rate_fPrevious"] <- x$Rate_f[sel]
    x[sel2, "Rate_mPrevious"] <- x$Rate_m[sel]
    x[sel2, "RatePrevious"] <- x$Rate[sel]
    x[sel2, c("Rate_m", "Rate_f", "Rate")] <- NA
    # x[sel, ]
    # x[sel2, ]
    if (abs(sum(x$Rate, na.rm = TRUE) - initialRateSum) > 1e-10) {
      print(marker)
      print(paste(count, "of", length(remove_markers)))
      print("initialRateSum")
      print(initialRateSum)
      print("sum(x$Rate)")
      print(sum(x$Rate, na.rm = TRUE))
      print("initialRateSum - sum(x$Rate)")
      print(initialRateSum - sum(x$Rate, na.rm = TRUE))
      print("head(removed)")
      print(head(removed))
      print("x[sel, ]")
      print(x[sel, ])
      print("map[map$NameNext == marker, ]")
      sel <- map$NameNext == marker # find the marker before the one we will remove
      sel[is.na(sel)] <- FALSE
      print(map[sel, ])
      stop("The sum of recombination rates is less than in the original map!")
      # which(map$Name == "UA-IFASA-6857")
      # val <- 58969
      # map[c(val - 1, val, val + 1),]
      #                     Name Chr Position       Rate_f       Rate_m         Rate       NamePrevious PositionPrevious Rate_fPrevious Rate_mPrevious
      # 58972 BovineHD2900012528  29 41419488 1.990877e-04 1.629160e-06 1.003584e-04       BTB-01541976         41407366   1.990877e-04   2.257184e-07
      # 58973      UA-IFASA-6857  29 41470852 1.278465e-05 3.657065e-08 6.410609e-06 BovineHD2900012528         41419488   1.990877e-04   1.629160e-06
      # 58974       BTB-01029026  29 41512425 4.990333e-08 1.107515e-03 5.537827e-04      UA-IFASA-6857         41470852   1.278465e-05   3.657065e-08
      #       RatePrevious              NameNext PositionNext   Rate_fNext   Rate_mNext     RateNext
      # 58972 9.965670e-05         UA-IFASA-6857     41470852 1.278465e-05 3.657065e-08 6.410609e-06
      # 58973 1.003584e-04          BTB-01029026     41512425 4.990333e-08 1.107515e-03 5.537827e-04
      # 58974 6.410609e-06 Hapmap43902-BTA-65773     41552472 2.388395e-03 1.781538e-05 1.203105e-03
      # which(x$Name == "UA-IFASA-6857")
      # val <- 58968
      # x[c(val - 1, val, val + 1),]
    }
  }
  # x[sel_remove_markers, c("Name", "NamePrevious", "Rate", "RatePrevious", "RateNext")]
  # x[sel_previous_markers, c("Name", "Rate", "RatePrevious", "RateNext")]
  return(list(kept = x, removed = removed))
}

map2 <- move_rec_rate_to_previous_marker(
  x = map,
  markers_to_remove = map$Name[!test_presence]
)
sum(map$Rate, na.rm = TRUE) # 24.35493
nrow(map) # 59309
sum(is.na(map$Rate)) # 0
sum(map2$kept$Rate, na.rm = TRUE) # 24.35493
nrow(map2$kept) # 59309
sum(is.na(map2$kept$Rate)) # 4
sum(map$Rate, na.rm = TRUE) - sum(map2$kept$Rate, na.rm = TRUE) # 3.552714e-15, so ~0
sum(map2$removed$Rate, na.rm = TRUE) # 0.0004469681
nrow(map2$removed) # 4

sel <- map$Name %in% map2$removed$Name
cols <- c(
  "Name",
  "Position",
  "Rate",
  "NamePrevious",
  "PositionPrevious",
  "RatePrevious"
)
map[sel, cols]
sel2 <- map2$kept$Name %in% map2$removed$NamePrevious
map2$kept[sel2, c("Name", "Rate")]
map <- map2$kept

# Explore more how markers map onto the new assemblt
sel_markers_in_assembly <- assembly$marker_name %in% map$Name
sum(sel_markers_in_assembly) # 507373
# Huh, we have markers multiple times in this assembly mapping file ...
# Let's see what's going on ...
# Here is one marker ...
sel <- assembly$marker_name == "BovineHD0100000035"
assembly[sel, c("marker_name", "assay", "ars120_pos")]
# Right, it's mapping for multiple SNP arrays.
# All values for ars120_pos are the same for this marker.
# Does this hold for all markers? Below is code to check this. The answer is no.
assembly <- assembly[sel_markers_in_assembly, ]
nrow(assembly) # 507373
if (FALSE) {
  for (marker in unique(assembly$marker_name)) {
    # marker = "BovineHD2300015209"
    sel <- assembly$marker_name == marker
    # assembly[sel, c("marker_name", "assay", "ars120_pos")]
    tmp <- assembly$ars120_pos[sel]
    tmp <- tmp[!is.na(tmp)]
    tmp <- tmp[tmp != 0]
    if (length(unique(tmp)) != 1) {
      print(assembly[sel, c("marker_name", "assay", "ars120_pos")])
      print(var(tmp))
    }
  }
}
# Let's take the HD snp array as a "reference" since Ma et al. worked with it in
# addition to all other SNP arrays, but the HD one is the largest so contains
# most/all SNP markers on the other arrays as well.
sel <- assembly$assay == "HD" &
  (assembly$ars120_pos != 0 | !is.na(assembly$ars120_pos))
sum(sel) # 55278
assemblyHD <- assembly[sel, ]
test_assemblyHD_in_map <- assemblyHD$marker_name %in% map$Name
sum(test_assemblyHD_in_map) # 55278
test_map_in_assemblyHD <- map$Name %in% assemblyHD$marker_name
sum(test_map_in_assemblyHD) # 55278
sum(!test_map_in_assemblyHD) # 4027
# OK, we went from 59309/59305 to 55278 markers, so now we have to move rec rate
# for 4027 markers to their previous markers.
markers_not_in_assemblyHD <- map$Name[!test_map_in_assemblyHD]
length(markers_not_in_assemblyHD) # 4027
markers_in_assemblyHD <- map$Name[test_map_in_assemblyHD]
length(markers_in_assemblyHD) # 55278
sum(map$Rate[test_map_in_assemblyHD]) # 22.67355
sum(map$Rate[!test_map_in_assemblyHD]) # 1.68138
# 22.67355 + 1.68138 = 24.35493
map2 <- move_rec_rate_to_previous_marker(
  x = map,
  markers_to_remove = markers_not_in_assemblyHD
)
sum(map$Rate, na.rm = TRUE) # 24.35493
nrow(map) # 59309
sum(is.na(map$Rate)) # 4
sum(map2$kept$Rate, na.rm = TRUE) # 24.35493
nrow(map2$kept) # 59308
sum(is.na(map2$kept$Rate)) # 4030
sum(map$Rate, na.rm = TRUE) - sum(map2$kept$Rate, na.rm = TRUE) # -7.81597e-14
sum(map2$removed$Rate, na.rm = TRUE) # 1.803793 (more than 1.68138 due to recursion!?)
nrow(map2$removed) # 4031
sum(map2$removed$Rate, na.rm = TRUE) -
  sum(map$Rate[!test_map_in_assemblyHD], na.rm = TRUE) # 0.1224109
# 1.803793 - 1.68138 = 0.122413

summary(map2$removed$Position)
summary(map2$kept$Position)
map2$kept[map2$kept$Position == 0, ] # let's see if position moved in new assembly
map <- map2$kept

assemblyHD <- assemblyHD[, c("marker_name", "ars120_pos", "ars120_chrn")]
head(assemblyHD)
colnames(assemblyHD) <- c("Name", "PositionNew", "ChrNew")
map <- merge(x = map, y = assemblyHD, by = "Name", all.x = TRUE)
nrow(map) # 59308
head(map)
sum(map$Rate, na.rm = TRUE) # 24.35493

# Let's now explore how many markers have moved chromosomes
map$ChrDiff <- abs(map$Chr - map$ChrNew)
summary(map$ChrDiff)
head(map[order(map$ChrDiff, decreasing = TRUE), ])
table(map$Chr)
table(map$ChrNew) # 2 markers on chr 30(X) and 6 markers on chr 99(some contig?)
map2 <- move_rec_rate_to_previous_marker(
  x = map,
  markers_to_remove = map$Name[map$ChrNew > 29]
)
sum(map2$kept$Rate, na.rm = TRUE) # 24.35493
nrow(map2$kept) # 59308
sum(is.na(map2$kept$Rate)) # 4038
map <- map2$kept

map$ChrDiff <- abs(map$Chr - map$ChrNew)
summary(map$ChrDiff)
head(map[order(map$ChrDiff, decreasing = TRUE), ])
map[order(map$ChrDiff, decreasing = TRUE), ][1:100, ]
sel <- map$ChrDiff > 0 |
  is.na(map$ChrNew) |
  map$PositionNew == 0 |
  is.na(map$PositionNew)
sum(sel) # 4141
map2 <- move_rec_rate_to_previous_marker(
  x = map,
  markers_to_remove = map$Name[sel]
)
sum(map2$kept$Rate, na.rm = TRUE) # 24.35493
nrow(map2$kept) # 59308
sum(is.na(map2$kept$Rate)) # 4141
map <- map2$kept

# Let's now explore how many markers have moved within chromosomes
map$PositionDiff <- abs(map$Position - map$PositionNew)
summary(map$PositionDiff)
hist(map$PositionDiff)
# some moved a little, some moved much more

# Let's check order changes, but first we assign the order and cumulative map
cumsum2 <- function(x, na.rm = FALSE) {
  y <- x
  y[!is.na(y)] <- cumsum(y[!is.na(y)])
  return(y)
}
chrs <- sort(unique(map$Chr))
map <- map[order(map$Chr, map$Position), ]
map$Order <- NA
map$Map_f <- NA
map$Map_m <- NA
map$Map <- NA
for (chr in chrs) {
  # chr <- 1
  sel <- map$Chr == chr
  map$Order[sel] <- 1:sum(sel)
  map$Map_f[sel] <- cumsum2(map$Rate_f[sel], na.rm = TRUE)
  map$Map_m[sel] <- cumsum2(map$Rate_m[sel], na.rm = TRUE)
  map$Map[sel] <- cumsum2(map$Rate[sel], na.rm = TRUE)
}

map <- map[order(map$ChrNew, map$PositionNew), ]
map$OrderNew <- NA
map$Map_fNew <- NA
map$Map_mNew <- NA
map$MapNew <- NA
for (chr in chrs) {
  sel <- map$Chr == chr
  map$OrderNew[sel] <- 1:sum(sel)
  map$Map_fNew[sel] <- cumsum2(map$Rate_f[sel], na.rm = TRUE)
  map$Map_mNew[sel] <- cumsum2(map$Rate_m[sel], na.rm = TRUE)
  map$MapNew[sel] <- cumsum2(map$Rate[sel], na.rm = TRUE)
}

head(map)
head(map[sel, ])
tail(map)
summary(map$Map_f)
summary(map$Map_fNew)
summary(map$Map_m)
summary(map$Map_mNew)
summary(map$Map)
summary(map$MapNew)

map$OrderDiff <- abs(map$Order - map$OrderNew)
summary(map$OrderDiff)
hist(map$OrderDiff)
hist(log(map$OrderDiff))
sum(map$OrderDiff > 0) # 58856
# Of course the above is a large number because even inserting one marker will
# shift order of all markers after it! The below new-vs-old position plots
# look reasonable - some reordering, but not much. We just need a way to pick out
# most outlying markers.

sel <- !is.na(map$PositionNew) & !is.na(map$Position)
plot(y = map$PositionNew[sel], x = map$Position[sel])
fit <- lm(PositionNew ~ Position, data = map[sel, ])
abline(fit, col = "red")
predictions <- predict(fit, interval = "confidence")
lines(map$Position[sel], predictions[, "lwr"], col = "green", lty = 2)
lines(map$Position[sel], predictions[, "upr"], col = "green", lty = 2)
predictions <- predict(fit, interval = "prediction")
lines(map$Position[sel], predictions[, "lwr"], col = "blue", lty = 2)
lines(map$Position[sel], predictions[, "upr"], col = "blue", lty = 2)

nPoints <- 0
pointedMarkers <- NULL
for (chr in chrs) {
  sel <- !is.na(map$PositionNew) & !is.na(map$Position) & map$Chr == chr
  sum(sel)
  plot(y = map$PositionNew[sel], x = map$Position[sel], main = chr, cex = 0.2)
  fit <- lm(PositionNew ~ Position, data = map[sel, ])
  abline(fit, col = "red")
  predictions <- predict(fit, interval = "confidence")
  lines(map$Position[sel], predictions[, "lwr"], col = "green", lty = 2)
  lines(map$Position[sel], predictions[, "upr"], col = "green", lty = 2)
  predictions <- predict(fit, interval = "prediction")
  lines(map$Position[sel], predictions[, "lwr"], col = "blue", lty = 2)
  lines(map$Position[sel], predictions[, "upr"], col = "blue", lty = 2)

  # Identify points outside the prediction intervals
  outside_bounds <- (map$PositionNew[sel] < predictions[, "lwr"]) |
    (map$PositionNew[sel] > predictions[, "upr"])
  points(
    x = map$Position[sel][outside_bounds],
    y = map$PositionNew[sel][outside_bounds],
    col = "red",
    pch = 19
  )
  print(paste("Chr:", chr))
  nPoints <- nPoints + sum(outside_bounds)
  print(paste("larger changes:", sum(outside_bounds), nPoints))
  pointedMarkers <- c(pointedMarkers, map$Name[sel][outside_bounds])
  dev.copy(png, file = paste0("chr", chr, ".png"))
  dev.off()
}

map2 <- move_rec_rate_to_previous_marker(
  x = map,
  markers_to_remove = pointedMarkers
)
# TODO: this function fails because on chromosome 28 we have a large block at
#       the start of the chromosome and we keep removing that block from right
#       the left (towards the start) and then we hit the first marker - we could
#       remove or keep that maker - have to decide to find another map!
sum(map2$kept$Rate, na.rm = TRUE) # 24.35493
nrow(map2$kept) # 59308
sum(is.na(map2$kept$Rate)) # 4038
map <- map2$kept


for (chr in chrs) {
  sel <- map$Chr == chr
  plot(
    y = map$MapNew[sel],
    x = map$PositionNew[sel],
    main = chr,
    type = "p",
    cex = 0.1
  )
}
for (chr in chrs) {
  sel <- map$Chr == chr
  plot(
    y = map$Rate[sel],
    x = map$Position[sel],
    main = chr,
    type = "p",
    cex = 0.1
  )
  points(
    y = map$Rate[sel],
    x = map$PositionNew[sel],
    main = chr,
    cex = 0.1,
    col = "red"
  )
}

# TODO: How can I pick out markers that have moved a lot?

sum(map$Rate_f, na.rm = TRUE)
# before 23.18364
sum(map$Rate_m, na.rm = TRUE)
# before 25.52622
sum(map$Rate, na.rm = TRUE)
# before 24.35493

# TODO: Export to RateMap format
