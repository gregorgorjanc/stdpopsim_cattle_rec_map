# setwd(dir = "~/Storages/GitBox/popsim/stdpopsim/stdpopsim_cattle_rec_map/Ma_et_al_2015/")

# This script imports the recombination map from Ma et al. (2015), which
# was based on a modified UMD3.1 assembly, and remaps it to the new ARS-UCD1.2.
# During this process we loose ~4K markers, which might not sound too bad, but
# this process shrunks the total genetic length from from 24.3 to 22.6, which
# is not what we want.
# After discussing this independetly with Li Ma (e-mail) and Peter Ralph (GitHub),
# they both suggested to add the recombination rate from the removed markers to
# the marker just before it. This is what we will do in the next script.

# ---- Import and inspect the map ----

map <- read.table(file = "cattle_rmap.txt", header = TRUE)
nrow(map)
head(map)
map <- map[order(map$Chr, map$Location), ]
head(map)
tail(map)
summary(map$map_f)
summary(map$map_m)
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

# ---- Import the new assembly and remap the map to it ----
# (ARS-UCD1.2 assembly, the above map is based on the modified UMD3.1 assembly)

assembly <- read.csv(
  file = "../Schnabel_2018/UMC_marker_names_180910/9913_ARS1.2_raw_data_180910.csv",
  header = TRUE
)

test_presence <- map$Name %in% assembly$marker_name
sum(test_presence) # 59305
sum(!test_presence) # 4
# OK, does not look too bad.
# We will drop the 4 markers that are not present in the new assembly.
# TODO: revise what to do with the "shrunken" map

# Explore more
sel_markers_in_assembly <- assembly$marker_name %in% map$Name
sum(sel_markers_in_assembly) # 507373
# Huh, we have markers multiple times ... let's see what's going on
# Here is one marker ...
sel <- assembly$marker_name == "BovineHD0100000035"
assembly[sel, c("marker_name", "assay", "ars120_pos")]
# and all values for ars120_pos are the same - does this hold for all markers?
assembly <- assembly[sel_markers_in_assembly, ]
nrow(assembly) # 507373
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
# Let's take the HD snp array as a "reference" since Ma et al. worked with it in
# addition to all other SNP arrays, but the HD one is the largest so contains
# most/all SNP markers on the other arrays as well.
sel <- assembly$assay == "HD" & assembly$ars120_pos != 0
assemblyHD <- assembly[sel, ]
sel_markers_in_assemblyHD <- assemblyHD$marker_name %in% map$Name
sum(sel_markers_in_assemblyHD) # 55260
test_assemblyHD_in_markers <- map$Name %in% assemblyHD$marker_name
sum(test_assemblyHD_in_markers) # 55260
# OK, we went from 59309 to 55260 markers, which is not too bad.
# TODO: revise what to do with the "shrunken" map

assemblyHD <- assemblyHD[, c("marker_name", "ars120_pos", "ars120_chrn")]
head(assemblyHD)
colnames(assemblyHD) <- c("Name", "PositionNew", "ChrNew")
map <- merge(x = map, y = assemblyHD, by = "Name")
nrow(map) # 55260
head(map)

# Let's now explore how many markers have moved chromosomes
map$ChrDiff <- abs(map$Chr - map$ChrNew)
summary(map$ChrDiff)
head(map[order(map$ChrDiff, decreasing = TRUE), ])
table(map$Chr)
table(map$ChrNew) # 2 markers on chr 30(X) and 6 markers on chr 99(some contig?)
map <- map[map$ChrNew < 30, ]
map$ChrDiff <- abs(map$Chr - map$ChrNew)
summary(map$ChrDiff)
head(map[order(map$ChrDiff, decreasing = TRUE), ])
map[order(map$ChrDiff, decreasing = TRUE), ][1:100, ]
sum(map$ChrDiff > 0) # 85, not that bad
map <- map[map$ChrDiff < 1, ]
nrow(map) # 55167
# TODO: revise what to do with the "shrunken" map

# Let's now explore how many markers have moved within chromosomes
map$PositionDiff <- abs(map$Position - map$PositionNew)
summary(map$PositionDiff)
hist(map$PositionDiff)
# some moved a little, some moved much more
# TODO: revise what to do with the "shrunken" map

chrs <- sort(unique(map$Chr))
map <- map[order(map$Chr, map$Position), ]
map$Order <- NA
map$Map_f <- NA
map$Map_m <- NA
map$Map <- NA
for (chr in chrs) {
  sel <- map$Chr == chr
  map$Order[sel] <- 1:sum(sel)
  map$Map_f[sel] <- cumsum(map$Rate_f[sel])
  map$Map_m[sel] <- cumsum(map$Rate_m[sel])
  map$Map[sel] <- cumsum(map$Rate[sel])
}

map <- map[order(map$ChrNew, map$PositionNew), ]
map$OrderNew <- NA
map$Map_fNew <- NA
map$Map_mNew <- NA
map$MapNew <- NA
for (chr in chrs) {
  sel <- map$Chr == chr
  map$OrderNew[sel] <- 1:sum(sel)
  map$Map_fNew[sel] <- cumsum(map$Rate_f[sel])
  map$Map_mNew[sel] <- cumsum(map$Rate_m[sel])
  map$MapNew[sel] <- cumsum(map$Rate[sel])
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

plot(y = map$PositionNew, x = map$Position)
for (chr in chrs) {
  sel <- map$Chr == chr
  plot(y = map$PositionNew[sel], x = map$Position[sel], main = chr, cex = 0.2)
}
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

map$OrderDiff <- abs(map$Order - map$OrderNew)
summary(map$OrderDiff)
hist(map$OrderDiff)
hist(log(map$OrderDiff))
sum(map$OrderDiff > 0) # 35172
# Of course the above is a large number because even inserting one marker will
# shift order of all markers after it! The above new-vs-old position plots don't
# look reasonable - some reordering, but not much.

sum(map$Rate_f)
# before 23.18364, now 21.52453
sum(map$Rate_m)
# before 25.52622, now 23.72068
sum(map$Rate)
# before 24.35493, now 22.62261

# TODO: Export to RateMap format
