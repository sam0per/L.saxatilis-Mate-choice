# fix CZA, CZC, CZD
summary(CZ_all$shore)
isl = "CZD"
colnames(CZ_all)
colnames(mfn)
CZ_LCP = read.csv("data/CZD_spatial_LCP_20190110.csv")
colnames(CZ_LCP)
CZ = merge(mfn, CZ_LCP, by = "snail_ID")
plot(CZ$LCmeanDist, CZ$length_mm)
plot(CZ$LCmeanDist, CZ$PC1)
plot(CZ$LCmeanDist, CZ$PC2)
plot(CZ_all[CZ_all$shore==isl,]$LCmeanDist, CZ_all[CZ_all$shore==isl,]$shape)
CZ$shape = CZ$PC1
CZ$shore = isl
CZ$test_sex = CZ$sex
intersect(colnames(CZ), colnames(CZ_all))
setdiff(colnames(CZ_all), colnames(CZ))
subCZ = CZ[, intersect(colnames(CZ), colnames(CZ_all))]
# subCZ$X = CZ$X.y
head(subCZ)
summary(subCZ)
write.table(subCZ, "old/CZD/final_data/CZD_use_cleanup3.csv", row.names = FALSE, col.names = TRUE, sep = ",")

noCZ = CZ_all[!CZ_all$shore==isl, ]
CZ_fin = noCZ[, intersect(colnames(noCZ), colnames(subCZ))]
CZ_fin = arrange(rbind(CZ_fin, subCZ), shore)
summary(CZ_fin)
CZ_fin = CZ_fin[-(which(CZ_fin$ref_size<2)), ]

# write.table(CZ_fin, "data/CZ_all_mating_clean.csv", row.names = FALSE, col.names = TRUE, sep = ",")

# CZ_data = read.csv("data/CZ_all_mating_clean.csv", sep = ",")
write.table(CZ_data, "data/CZ_all_mating_clean.csv", row.names = FALSE, col.names = TRUE, sep = ";")
CZ_mate = read.csv("data/CZ_all_mating_clean.csv", sep = ";")
colnames(CZ_mate)
islands = c("CZA", "CZB", "CZC", "CZD")
# CZ_spatial = read.csv(list.files(path = "data", pattern = paste0("CZD", "_spatial"), full.names = TRUE))
# head(CZ_spatial)


add_spatial_LCP = function(isl, merged_df) {
  cz_spat = read.csv(list.files(path = "data", pattern = paste0(isl, "_spatial"), full.names = TRUE))
  merged_df = merge(merged_df, cz_spat, by = "snail_ID")
  return(merged_df)}

# summary(CZ_mate)
# merged_tmp = add_spatial_LCP(isl = "CZD", merged_df = CZ_mate)
# summary(merge(CZ_mate, CZ_spatial, by = "snail_ID"))

CZ_df_list = lapply(islands, function(x) {
  add_spatial_LCP(isl = x, merged_df = CZ_mate)
})
glimpse(CZ_df_list)

multi_join <- function(list_of_loaded_data, join_func, ...){
  
  require("dplyr")
  
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  
  return(output)
}

CZ_merged = multi_join(list_of_loaded_data = CZ_df_list, join_func = full_join)
(discrep = mapply(setdiff, CZ_mate, CZ_merged[, 1:ncol(CZ_mate)]))
head(CZ_merged)
summary(CZ_merged)

# CZ_data = CZ_merged[is.na(CZ_merged$LCmeanDist)==FALSE & is.na(CZ_merged$size_ratio)==FALSE, ]
# CZ_data = CZ_merged[is.na(CZ_merged$size_ratio)==FALSE, ]
# CZ_tmp = CZ_merged[is.na(CZ_merged$size_ratio)==FALSE, ]
# CZ_tmp[is.na(CZ_tmp$LCmeanDist)==TRUE, ]

CZ_merged$X.1 = NULL
# summary(CZ_data)
# summary(CZ_mate)
# length(unique(CZ_spatial$snail_ID))
# length(unique(CZ_mate[CZ_mate$shore=="CZD",]$snail_ID))
# colnames(CZ_data)
# colnames(CZ_mate)

write.table(CZ_merged, file = "data/CZ_all_mating_clean.csv", sep = ",", row.names = FALSE)


CZ_data = read.csv("data/CZ_all_mating_clean.csv")
(islands = unique(as.character(CZ_data$shore)))

lapply(islands, function(x) {
  with(data = CZ_data[CZ_data$shore==x, ], plot(x = LCmeanDist, y = length_mm))
})

CZ_snail_pos = lapply(islands, function(x) {
  ggplot(data = CZ_data[CZ_data$shore==x, ], aes(x = X, y = Y)) +
    geom_point(col='red', size = 3) +
    labs(title = x) +
    theme(rect = element_rect(fill = "transparent"))
})

ggsave(filename = paste0("figures/CZA_snail_pos.png"), plot = CZ_snail_pos[[1]],bg = "transparent")
ggsave(filename = paste0("figures/CZB_snail_pos.png"), plot = CZ_snail_pos[[2]],bg = "transparent")
ggsave(filename = paste0("figures/CZC_snail_pos.png"), plot = CZ_snail_pos[[3]],bg = "transparent")
ggsave(filename = paste0("figures/CZD_snail_pos.png"), plot = CZ_snail_pos[[4]],bg = "transparent")
