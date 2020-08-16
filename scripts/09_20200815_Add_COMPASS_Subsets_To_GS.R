library(here)
library(tidyverse)
library(flowWorkspace)

# The counts for the COMPASS subsets are only stored in COMPASSResult objects if they were discovered by COMPASS for that stim,
# so we have to manually add boolean gates for each subset to the GatingSets and then extract the count data later.

gsPath <- here::here("out/GatingSets/20200805_HAARVI_ICS_GatingSet_AllBatches")
gs <- load_gs(gsPath)

##################################################

merged_cd4_compass_data <- readRDS("processed_data/20200815_Merged_CD4_COMPASS_Data.rds")
merged_cd8_compass_data <- readRDS("processed_data/20200815_Merged_CD8_COMPASS_Data.rds")

##################################################

# Add the CD4+ COMPASS boolean cytokine subset gates to the GatingSet
mapMarkers <- c("IL2", "IL4/5/13", "IFNg", "TNFa", "IL17a", "CD154", "CD107a")
cd4NodeMarkerMap <- mapMarkers
# NodeMarkerMap names are gating tree paths
names(cd4NodeMarkerMap) <- paste0("4+", "/", c("IL2", "IL4513", "IFNG", "TNF", "IL17", "154", "107a"))

cd4_cats_mod <- as.data.frame(merged_cd4_compass_data$catsMerged) %>% 
  mutate_all(~ as.numeric(as.character(.))) %>% 
  mutate_all(~ recode(., "0" = "!", "1" = "")) %>% 
  dplyr::rename_at(vars(cd4NodeMarkerMap), ~ names(cd4NodeMarkerMap))
cd4_booleanSubsets <- cd4_cats_mod %>% 
  rowwise() %>% 
  do(booleanSubset = paste(paste0(., colnames(cd4_cats_mod)), collapse="&")) %>% 
  ungroup() %>% 
  dplyr::pull(booleanSubset) %>% 
  unlist()
names(cd4_booleanSubsets) <- paste0("CD4_", gsub("4\\+\\/", "", gsub("\\&", "_AND_", gsub("\\!", "NOT_", cd4_booleanSubsets))))
for(booleanSubsetName in names(cd4_booleanSubsets)) {
  # booleanSubset The booleanSubset (a combination of existing gates) in string format, e.g. "8+/GMM+&!8+/GAMMADELTA"
  call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(cd4_booleanSubsets[[booleanSubsetName]])))
  g <- eval(call)
  suppressWarnings(flowWorkspace::gs_pop_add(gs, g, parent = "4+", name=booleanSubsetName))
}

dput(names(cd4_booleanSubsets))
c("CD4_NOT_107a_AND_NOT_154_AND_NOT_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_TNF", 
  "CD4_NOT_107a_AND_NOT_154_AND_NOT_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_IL4513_AND_NOT_TNF", 
  "CD4_NOT_107a_AND_NOT_154_AND_NOT_IFNG_AND_NOT_IL17_AND_IL2_AND_NOT_IL4513_AND_NOT_TNF", 
  "CD4_NOT_107a_AND_NOT_154_AND_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_NOT_TNF", 
  "CD4_NOT_107a_AND_154_AND_NOT_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_NOT_TNF", 
  "CD4_107a_AND_NOT_154_AND_NOT_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_NOT_TNF", 
  "CD4_NOT_107a_AND_NOT_154_AND_NOT_IFNG_AND_NOT_IL17_AND_IL2_AND_IL4513_AND_NOT_TNF", 
  "CD4_NOT_107a_AND_NOT_154_AND_NOT_IFNG_AND_IL17_AND_NOT_IL2_AND_IL4513_AND_NOT_TNF", 
  "CD4_NOT_107a_AND_NOT_154_AND_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_TNF", 
  "CD4_NOT_107a_AND_154_AND_NOT_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_TNF", 
  "CD4_NOT_107a_AND_154_AND_NOT_IFNG_AND_NOT_IL17_AND_IL2_AND_NOT_IL4513_AND_NOT_TNF", 
  "CD4_NOT_107a_AND_154_AND_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_NOT_TNF", 
  "CD4_NOT_107a_AND_NOT_154_AND_IFNG_AND_NOT_IL17_AND_IL2_AND_NOT_IL4513_AND_TNF", 
  "CD4_NOT_107a_AND_NOT_154_AND_IFNG_AND_NOT_IL17_AND_IL2_AND_IL4513_AND_NOT_TNF", 
  "CD4_NOT_107a_AND_154_AND_NOT_IFNG_AND_NOT_IL17_AND_IL2_AND_NOT_IL4513_AND_TNF", 
  "CD4_NOT_107a_AND_154_AND_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_TNF", 
  "CD4_NOT_107a_AND_154_AND_IFNG_AND_NOT_IL17_AND_IL2_AND_NOT_IL4513_AND_NOT_TNF", 
  "CD4_NOT_107a_AND_154_AND_NOT_IFNG_AND_NOT_IL17_AND_IL2_AND_IL4513_AND_TNF", 
  "CD4_NOT_107a_AND_154_AND_IFNG_AND_NOT_IL17_AND_IL2_AND_NOT_IL4513_AND_TNF", 
  "CD4_NOT_107a_AND_154_AND_IFNG_AND_NOT_IL17_AND_IL2_AND_IL4513_AND_TNF", 
  "CD4_107a_AND_154_AND_IFNG_AND_NOT_IL17_AND_IL2_AND_NOT_IL4513_AND_TNF"
)

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                list(v = as.symbol(paste(names(cd4_booleanSubsets), collapse="|"))))),
           parent = "/Time/LD-3+/1419-3+/S/Lymph/4+", name = "CD4_COMPASS_Subsets")

#############################

# Add the CD8 COMPASS boolean cytokine subset gates to the GatingSet
mapMarkers <- c("IL2", "IL4/5/13", "IFNg", "TNFa", "IL17a", "CD154", "CD107a")
cd8NodeMarkerMap <- mapMarkers
# NodeMarkerMap names are gating tree paths
names(cd8NodeMarkerMap) <- paste0("8+", "/", c("IL2", "IL4513", "IFNG", "TNF", "IL17", "154", "107a"))

cd8_cats_mod <- as.data.frame(merged_cd8_compass_data$catsMerged) %>%
  mutate_all(~ as.numeric(as.character(.))) %>% 
  mutate_all(~ recode(., "0" = "!", "1" = "")) %>%
  dplyr::rename_at(vars(cd8NodeMarkerMap), ~ names(cd8NodeMarkerMap))
cd8_booleanSubsets <- cd8_cats_mod %>%
  rowwise() %>%
  do(booleanSubset = paste(paste0(., colnames(cd8_cats_mod)), collapse="&")) %>%
  ungroup() %>%
  dplyr::pull(booleanSubset) %>%
  unlist()
names(cd8_booleanSubsets) <- paste0("CD8_", gsub("8\\+\\/", "", gsub("\\&", "_AND_", gsub("\\!", "NOT_", cd8_booleanSubsets))))
for(booleanSubsetName in names(cd8_booleanSubsets)) {
  # booleanSubset The booleanSubset (a combination of existing gates) in string format, e.g. "8+/GMM+&!8+/GAMMADELTA"
  call <- substitute(flowWorkspace::booleanFilter(v), list(v = as.symbol(cd8_booleanSubsets[[booleanSubsetName]])))
  g <- eval(call)
  suppressWarnings(flowWorkspace::gs_pop_add(gs, g, parent = "8+", name=booleanSubsetName))
}
dput(names(cd8_booleanSubsets))
c("CD8_NOT_107a_AND_NOT_154_AND_NOT_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_IL4513_AND_NOT_TNF", 
  "CD8_NOT_107a_AND_NOT_154_AND_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_NOT_TNF", 
  "CD8_107a_AND_NOT_154_AND_NOT_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_NOT_TNF", 
  "CD8_NOT_107a_AND_NOT_154_AND_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_TNF", 
  "CD8_NOT_107a_AND_NOT_154_AND_IFNG_AND_NOT_IL17_AND_IL2_AND_NOT_IL4513_AND_NOT_TNF", 
  "CD8_107a_AND_NOT_154_AND_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_NOT_TNF", 
  "CD8_107a_AND_NOT_154_AND_IFNG_AND_NOT_IL17_AND_NOT_IL2_AND_NOT_IL4513_AND_TNF"
)

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste(names(cd8_booleanSubsets), collapse="|"))))),
           parent = "/Time/LD-3+/1419-3+/S/Lymph/8+", name = "CD8_COMPASS_Subsets")

#############################

# Define the memory subpopulations under CD4_COMPASS_Subsets and CD8_COMPASS_Subsets
cd4_ccr7_path <- "/Time/LD-3+/1419-3+/S/Lymph/4+/CCR7+"
cd4_cd45ra_path <- "/Time/LD-3+/1419-3+/S/Lymph/4+/CD45RA+"
cd8_ccr7_path <- "/Time/LD-3+/1419-3+/S/Lymph/8+/CCR7+"
cd8_cd45ra_path <- "/Time/LD-3+/1419-3+/S/Lymph/8+/CD45RA+"

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_ccr7_path,
                                                         "&", cd4_cd45ra_path))))),
           parent = "/Time/LD-3+/1419-3+/S/Lymph/4+/CD4_COMPASS_Subsets", name = "Naive")
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd4_ccr7_path,
                                                         "&!", cd4_cd45ra_path))))),
           parent = "/Time/LD-3+/1419-3+/S/Lymph/4+/CD4_COMPASS_Subsets", name = "TCM")
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd4_ccr7_path,
                                                         "&", cd4_cd45ra_path))))),
           parent = "/Time/LD-3+/1419-3+/S/Lymph/4+/CD4_COMPASS_Subsets", name = "TEMRA")
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd4_ccr7_path,
                                                         "&!", cd4_cd45ra_path))))),
           parent = "/Time/LD-3+/1419-3+/S/Lymph/4+/CD4_COMPASS_Subsets", name = "TEM")

gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd8_ccr7_path,
                                                         "&", cd8_cd45ra_path))))),
           parent = "/Time/LD-3+/1419-3+/S/Lymph/8+/CD8_COMPASS_Subsets", name = "Naive")
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("", cd8_ccr7_path,
                                                         "&!", cd8_cd45ra_path))))),
           parent = "/Time/LD-3+/1419-3+/S/Lymph/8+/CD8_COMPASS_Subsets", name = "TCM")
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd8_ccr7_path,
                                                         "&", cd8_cd45ra_path))))),
           parent = "/Time/LD-3+/1419-3+/S/Lymph/8+/CD8_COMPASS_Subsets", name = "TEMRA")
gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                               list(v = as.symbol(paste0("!", cd8_ccr7_path,
                                                         "&!", cd8_cd45ra_path))))),
           parent = "/Time/LD-3+/1419-3+/S/Lymph/8+/CD8_COMPASS_Subsets", name = "TEM")

# Add activation gate (HLA-DR+CD38+) under CD4_COMPASS_Subsets and CD8_COMPASS_Subsets
# First, copy activation gates up to the Live gate.
cd3_hladr_path <- "/Time/LD-3+/1419-3+/S/Lymph/CD38+"
cd3_cd38_path <- "/Time/LD-3+/1419-3+/S/Lymph/HLADR+"

for (parent_path in c("/Time/LD-3+/1419-3+/S/Lymph/4+/CD4_COMPASS_Subsets", "/Time/LD-3+/1419-3+/S/Lymph/8+/CD8_COMPASS_Subsets")) {
  gs_pop_add(gs, eval(substitute(flowWorkspace::booleanFilter(v),
                                 list(v = as.symbol(paste0("", cd3_hladr_path,
                                                           "&", cd3_cd38_path))))),
             parent = parent_path, name = "HLADR+CD38+")
}

#############################

lymph_gate <- "/Time/LD-3+/1419-3+/S/Lymph"
flowWorkspace::recompute(gs, lymph_gate) 

#############################

plot(gs, bool = T, fontsize = 10)

# One last thing before saving: Fix the Batch typo
pData(gs)$Batch <- ifelse(pData(gs)$`EXPERIMENT NAME` == "20200605_COVID_ICS-B3", 3, pData(gs)$Batch)

save_gs(gs, here::here("out/GatingSets/20200815_HAARVI_ICS_GatingSet_AllBatches_with_COMPASS_Subsets"), cdf = "symlink") 
