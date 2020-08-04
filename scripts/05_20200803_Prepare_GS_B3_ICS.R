library(openCyto) # 1.24.0
library(CytoML) # 1.12.0
library(flowCore) # required for description()
library(flowWorkspace) # required for gh_pop_get_data()
library(ggcyto) # devtools::install_github("RGLab/ggcyto", ref="ggplot3.3") for update_theme()
library(here)
library(tidyverse)
library(readxl)

# Read in Batch 3 workspace and prepare GatingSet.

date <- 20200803

xml_path_b3 <- here::here("data/20200605_COVID_ICS-B3_KY3.xml")
fcs_subfolder_b3 <- here::here("data/20200605_COVID_ICS-B3/")

ws_b3 <- open_flowjo_xml(xml_path_b3)
merge(fj_ws_get_sample_groups(ws_b3), fj_ws_get_samples(ws_b3), by = "sampleID")
fj_ws_get_keywords(ws_b3, 46)
names(fj_ws_get_keywords(ws_b3, 46))
keywords2import <- c("EXPERIMENT NAME", "$DATE", "SAMPLE ID", "PATIENT ID", "STIM", "WELL ID", "PLATE NAME")

sampleGroup <- "Samples"
gs_b3 <- flowjo_to_gatingset(ws_b3, name=sampleGroup, keywords=keywords2import,
                             path=fcs_subfolder_b3, extend_val=-10000)

# All samples have the same gating tree
pop_lists <- lapply(gs_b3, gh_get_pop_paths)
unique(pop_lists)

pData(gs_b3)$filename <- sapply(rownames(pData(gs_b3)), function(x) {
  gsub("/.+/", "", description(gh_pop_get_data(gs_b3[[x]]))$FILENAME)
}, USE.NAMES = F)
pData(gs_b3)$rowname <- rownames(pData(gs_b3))
pData(gs_b3)

# Read in the patient manifest
manifest <- read.csv(here::here("data/Seshadri_HAARVI_PBMC_manifest_merged_11June2020.csv"), check.names = F, stringsAsFactors = F)
# Add metadata to pData
pData_tmp <- pData(gs_b3) %>% 
  # GatingSet SAMPLE ID 07B is Record ID 551432 in the manifest. Use 551432 from here on out
  mutate(`SAMPLE ID` = ifelse(`SAMPLE ID` == "07B", "551432", `SAMPLE ID`)) %>% 
  mutate(`SAMPLE ID` = toupper(`SAMPLE ID`)) %>% 
  left_join(manifest %>%
              dplyr::select(`Record ID`, `Sample ID`, "Collection date", "Cell count", 
                            "Cohort", "Age", "Sex", "Race", "Hispanic?", "Days symptom onset to visit 1", 
                            "Pair ID", "Race_v2"),
            by = c("SAMPLE ID" = "Record ID")) %>% 
  mutate(Batch = 2)
rownames(pData_tmp) <- rownames(pData(gs_b3))
pData(gs_b3) <- pData_tmp

pData(gs_b3) %>% dplyr::select("SAMPLE ID", "Cohort", "Age") %>% unique() %>% arrange(Cohort, `SAMPLE ID`)

pData(parameters(gh_pop_get_data(gs_b3[[1]])))[,c(1, 2)]

# Add names to all channels
dput(unname(pData(parameters(gh_pop_get_data(gs_b3[[1]])))[,2]))
markernames_b3<- c("Time", "FSC-A", "FSC-H", "SSC-A", "SSC-H", "CD8b", "TNFa", "CD107a", 
                   "CD154", "CD3 ECD", "IL2", "CD4", "IL17a", 
                   "IL4/5/13", "CD14/CD19", "CCR7", "CD38", 
                   "L/D", "IFNg", "CD45RA", "HLADR")
names(markernames_b3) <- pData(parameters(gh_pop_get_data(gs_b3[[1]])))[,1]
# markernames_b3
markernames(gs_b3) <- markernames_b3
pData(parameters(gh_pop_get_data(gs_b3[[1]])))[,c(1, 2)]
# name      desc
# $P1      Time      Time
# $P2     FSC-A     FSC-A
# $P3     FSC-H     FSC-H
# $P4     SSC-A     SSC-A
# $P5     SSC-H     SSC-H
# $P6  <B710-A>      CD8b
# $P7  <B515-A>      TNFa
# $P8  <G780-A>    CD107a
# $P9  <G660-A>     CD154
# $P10 <G610-A>   CD3 ECD
# $P11 <G575-A>       IL2
# $P12 <R780-A>       CD4
# $P13 <R710-A>     IL17a
# $P14 <R660-A>  IL4/5/13
# $P15 <V780-A> CD14/CD19
# $P16 <V710-A>      CCR7
# $P17 <V610-A>      CD38
# $P18 <V510-A>       L/D
# $P19 <V450-A>      IFNg
# $P20 <U730-A>    CD45RA
# $P21 <U395-A>     HLADR

# Grab the cytokine and memory gates from NOT4+ and add under 8+ gate
gates_to_copy <- c("/Time/LD-3+/1419-3+/S/Lymph/NOT4+/107a", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/154", 
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL2", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL17",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL4513", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/TNF",
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CCR7+", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CD45RA+")
# "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IFNG", already added
for(path in gates_to_copy) {
  gs_pop_add(gs_b3, lapply(gs_b3, gh_pop_get_gate, y=path),
             parent = "/Time/LD-3+/1419-3+/S/Lymph/8+")
}
recompute(gs_b3, "/Time/LD-3+/1419-3+/S/Lymph/8+")

png(here::here(sprintf("out/QC/B3_GatingTree_%s.png", date)), width = 7, height = 5, units = "in", res = 300)
(plot(gs_b3, fontsize=15, bool=T))
dev.off()

save_gs(gs_b3, here::here("out/GatingSets/20200803_HAARVI_ICS_GatingSet_B3"))

# 6 samples for each patient
table(pData(gs_b3)$`SAMPLE ID`)
# 23 unique patients
length(unique(pData(gs_b3)$`SAMPLE ID`))

#####################################################################

# gs_b3 <- load_gs(here::here("out/GatingSets/20200803_HAARVI_ICS_GatingSet_B3"))

dput(gh_get_pop_paths(gs_b3))

# Perform QC on CD3 counts
cd3_path <- "/Time/LD-3+/1419-3+/S/Lymph"
cd4_path <- "/Time/LD-3+/1419-3+/S/Lymph/4+"
cd8_path <- "/Time/LD-3+/1419-3+/S/Lymph/8+"
not4_path <- "/Time/LD-3+/1419-3+/S/Lymph/NOT4+"
cd3_cd4_cd8_counts_b3 <- pData(gs_b3) %>% 
  left_join(gs_pop_get_count_fast(gs_b3, subpopulations = c(cd3_path, cd4_path, cd8_path, not4_path)) %>% 
              pivot_wider(id_cols = name, names_from = "Population", values_from = "Count") %>% 
              rename(CD3 = !!cd3_path,
                     CD4 = !!cd4_path,
                     CD8 = !!cd8_path,
                     NOT4 = !!not4_path),
            by = c("rowname" = "name")) %>% 
  dplyr::select("SAMPLE ID", "Cohort", "Age", "Sex", "Race", "Days symptom onset to visit 1", "Batch", "STIM", "CD3", "CD4", "CD8", "NOT4")

png(here::here(sprintf("out/QC/Counts/B3_CD3_Counts_%s.png", date)),
    width = 10, height = 6, units="in", res=300)
ggplot(cd3_cd4_cd8_counts_b3 %>% 
         mutate(Color = ifelse(CD3 < 10000, "red", "black")),
       aes(x=`SAMPLE ID`, y=CD3)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = Color), position=position_jitter(width = 0.05, height=0)) +
  geom_hline(aes(yintercept = 10000, color="red"), linetype="dashed") +
  scale_color_identity() +
  ggtitle("Batch 3 CD3 Counts") +
  theme(legend.position = "none")
dev.off()

png(here::here(sprintf("out/QC/Counts/B3_CD4_Counts_%s.png", date)),
    width = 10, height = 6, units="in", res=300)
ggplot(cd3_cd4_cd8_counts_b3 %>% 
         mutate(Color = ifelse(CD4 < 3000, "red", "black")),
       aes(`SAMPLE ID`, CD4)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(aes(color = Color), position=position_jitter(width = 0.05, height=0)) +
  geom_hline(aes(yintercept = 3000, color="red"), linetype="dashed") +
  scale_color_identity() +
  ggtitle("Batch 3 CD4 Counts")
dev.off()

png(here::here(sprintf("out/QC/Counts/B3_CD8_Counts_%s.png", date)),
    width = 10, height = 6, units="in", res=300)
ggplot(cd3_cd4_cd8_counts_b3 %>% 
         mutate(Color = ifelse(CD8 < 3000, "red", "black")),
       aes(`SAMPLE ID`, CD8, color = Color)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.05, height=0)) +
  geom_hline(aes(yintercept = 3000, color="red"), linetype="dashed") +
  scale_color_identity() +
  ggtitle("Batch 3 CD8 Counts")
dev.off()

png(here::here(sprintf("out/QC/Counts/B3_NOT4_Counts_%s.png", date)),
    width = 10, height = 6, units="in", res=300)
ggplot(cd3_cd4_cd8_counts_b3 %>% 
         mutate(Color = ifelse(NOT4 < 3000, "red", "black")),
       aes(`SAMPLE ID`, NOT4, color = Color)) +
  geom_boxplot(outlier.shape = NA) +
  geom_point(position=position_jitter(width = 0.05, height=0)) +
  geom_hline(aes(yintercept = 3000, color="red"), linetype="dashed") +
  scale_color_identity() +
  ggtitle("Batch 3 NOT4 Counts")
dev.off()

# Investigate the low count wells
cd3_cd4_cd8_counts_b3 %>% 
  dplyr::filter(CD3 < 10000 | CD4 < 3000 | CD8 < 3000 | NOT4 < 3000)
# SAMPLE ID           Cohort Age  Sex                                   Race Days symptom onset to visit 1 Batch    STIM   CD3   CD4 CD8 NOT4
# 1     BWT22             <NA>  NA <NA>                                   <NA>                          <NA>     2    DMSO 14574  9069 589 5505
# 2     BWT22             <NA>  NA <NA>                                   <NA>                          <NA>     2     SEB 14629  9673 592 4956
# 3     BWT22             <NA>  NA <NA>                                   <NA>                          <NA>     2    VEMP 13972  8739 517 5233
# 4     BWT22             <NA>  NA <NA>                                   <NA>                          <NA>     2 Spike 1 13232  8914 549 4318
# 5     BWT22             <NA>  NA <NA>                                   <NA>                          <NA>     2 Spike 2 13324  9038 491 4286
# 6     BWT22             <NA>  NA <NA>                                   <NA>                          <NA>     2    NCAP 16463 11246 580 5217
# 7      116C Non-hospitalized  68    M American Indian or Alaska Native,White                            55     2    DMSO   454   168 101  286

# 116C DMSO: 27.8% of events in LD-CD3+ gate, unlike the rest of the wells for this patient which are ~40-50%.
# The CD3 and CD4 counts are so low that we have no choice but to drop this patient's wells all-together from COMPASS.

# Also investigate BWT22, which has low-but-passing CD3 and CD4 counts
# BWT22: <2% of events make it through LD-CD3+ gate. Guava Day 1 and Day 2 viabilities were also comparatively low compared to other tubes.

###################################################################

pData(gs_b3)$STIM <- factor(pData(gs_b3)$STIM, levels = c("DMSO", "VEMP", "Spike 1", "Spike 2", "NCAP", "SEB"))

plotter <- function(myGS, currentGates, currentBins=1024, drawMarginalFilter=FALSE) {
  firstGate <- currentGates[[1]]
  currentGateBoundaries <- attributes(gh_pop_get_gate(myGS[[1]], firstGate))$boundaries
  currentXaxis <- colnames(currentGateBoundaries)[[1]]
  currentYaxis <- colnames(currentGateBoundaries)[[2]]
  parentGate <- sub("\\/[^\\/]+$", "", firstGate)
  
  g <- if(drawMarginalFilter) {
    ggcyto(myGS,
           aes(!!currentXaxis, !!currentYaxis),
           subset = if(parentGate == "") "root" else parentGate,
           filter = marginalFilter) 
  } else {
    ggcyto(myGS,
           aes(!!currentXaxis, !!currentYaxis),
           subset = if(parentGate == "") "root" else parentGate) # filter = marginalFilter cuts off some axes
  }
  
  g +
    geom_hex(bins=currentBins) +
    geom_gate(currentGates) +
    axis_x_inverse_trans() + axis_y_inverse_trans() +
    ggcyto_par_set(limits = "instrument") +
    facet_grid(`SAMPLE ID` ~ STIM, switch = "y") +
    theme(#plot.title = element_blank(),
      legend.position = "none",
      strip.text.x = element_text(size = 14, margin = margin(0,0,0,0, "cm")),
      panel.grid.major = ggplot2::element_blank()) +
    theme_bw(base_size=28) +
    geom_stats(size=8, alpha=0.4) +
    labs(title=currentGates[[1]], caption = "Batch 3")
}

gates2draw <- list("/Time", "/Time/LD-3+", "/Time/LD-3+/1419-3+", "/Time/LD-3+/1419-3+/S", 
                   "/Time/LD-3+/1419-3+/S/Lymph", "/Time/LD-3+/1419-3+/S/Lymph/4+",
                   "/Time/LD-3+/1419-3+/S/Lymph/8+",
                   "/Time/LD-3+/1419-3+/S/Lymph/4+/107a", "/Time/LD-3+/1419-3+/S/Lymph/4+/154", 
                   c("/Time/LD-3+/1419-3+/S/Lymph/4+/IFNG", "/Time/LD-3+/1419-3+/S/Lymph/4+/IL2"), 
                   "/Time/LD-3+/1419-3+/S/Lymph/4+/IL17", "/Time/LD-3+/1419-3+/S/Lymph/4+/IL4513", 
                   "/Time/LD-3+/1419-3+/S/Lymph/4+/TNF",
                   c("/Time/LD-3+/1419-3+/S/Lymph/CD38+", "/Time/LD-3+/1419-3+/S/Lymph/HLADR+"), 
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/107a", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/154", 
                   c("/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IFNG", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL2"), 
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL17", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL4513", 
                   "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/TNF",
                   "/Time/LD-3+/1419-3+/S/Lymph/8+/107a", "/Time/LD-3+/1419-3+/S/Lymph/8+/154", 
                   c("/Time/LD-3+/1419-3+/S/Lymph/8+/IFNG", "/Time/LD-3+/1419-3+/S/Lymph/8+/IL2"), 
                   "/Time/LD-3+/1419-3+/S/Lymph/8+/IL17", "/Time/LD-3+/1419-3+/S/Lymph/8+/IL4513", 
                   "/Time/LD-3+/1419-3+/S/Lymph/8+/TNF")
bins2draw <- list(128, 128, 128, 128,
                  128, 128, 128, 1024,
                  1024, 1024, 1024, 1024,
                  1024, 1024, 1024, 1024,
                  1024, 1024, 1024, 1024,
                  1024, 1024, 1024, 1024,
                  1024, 1024)
marginalFilter2draw <- c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, 
                         FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, 
                         FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE)

pmap(list(gates2draw, bins2draw, marginalFilter2draw), function(currentGates, currentBins, drawMarginalFilter) {
  ptids <- unique(pData(gs_b3)$`SAMPLE ID`)
  ptids_per_plot <- 5
  num_plots <- ceiling(length(ptids)/ptids_per_plot)
  parentGate <- sub(".*\\/([^\\/]+)\\/[^\\/]+$", "\\1", currentGates[[1]])
  for(plot_num in seq_along(1:num_plots)) { 
    current_ptids <- ptids[((plot_num-1)*ptids_per_plot + 1):min(((plot_num-1)*ptids_per_plot + ptids_per_plot), length(ptids))]
    png(filename = file.path(here::here("out/QC/FACS_Plots"),
                             sprintf("%s%s_B3_pt%s_%s.png",
                                     if(parentGate %in% c("NOT4+", "4+", "8+")) {paste0(parentGate, "_")} else {""},
                                     sub(".*\\/([^\\/]+$)", "\\1", currentGates[[1]]), plot_num,
                                     format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
        width = 1920, height = 210+174*length(current_ptids), units = "px")
    print(plotter(subset(gs_b3, `SAMPLE ID` %in% current_ptids),
                  currentGates, currentBins, drawMarginalFilter) +
            labs(caption=sprintf("Batch 3, pt %s", plot_num)))
    dev.off()
  }
})

# Memory gates were drawn on histograms and won't work with the above plotter function
plotter_mem <- function(myGS, myGates) {
  firstGate <- myGates[[1]]
  firstGateAttributes <- attributes(gh_pop_get_gate(myGS[[1]], firstGate))
  currentXaxis <- names(firstGateAttributes$min)
  
  secondGate <- myGates[[2]]
  secondGateAttributes <- attributes(gh_pop_get_gate(myGS[[1]], secondGate))
  currentYaxis <- names(secondGateAttributes$min)
  
  parentGate <- sub("\\/[^\\/]+$", "", firstGate)
  
  ggcyto(myGS,
         aes(!!currentXaxis, !!currentYaxis),
         subset = if(parentGate == "") "root" else parentGate,
         filter = marginalFilter) +
    geom_hex(bins=128) +
    geom_gate(myGates) +
    axis_x_inverse_trans() + axis_y_inverse_trans() +
    ggcyto_par_set(limits = "instrument") +
    facet_grid(`SAMPLE ID` ~ STIM, switch = "y") +
    theme(#plot.title = element_blank(),
      legend.position = "none",
      strip.text.x = element_text(size = 14, margin = margin(0,0,0,0, "cm")),
      panel.grid.major = ggplot2::element_blank()) +
    theme_bw(base_size=28) +
    geom_stats(size=8, alpha=0.4) +
    labs(title=myGates[[1]], caption = "Batch 3")
}

cd4_mem_gates <- c("/Time/LD-3+/1419-3+/S/Lymph/4+/CCR7+", "/Time/LD-3+/1419-3+/S/Lymph/4+/CD45RA+")
notcd4_mem_gates <-  c("/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CCR7+", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/CD45RA+")
cd8_mem_gates <- c("/Time/LD-3+/1419-3+/S/Lymph/8+/CCR7+", "/Time/LD-3+/1419-3+/S/Lymph/8+/CD45RA+")

for(currentGates in list(cd4_mem_gates, notcd4_mem_gates, cd8_mem_gates)) {
  ptids <- unique(pData(gs_b3)$`SAMPLE ID`)
  ptids_per_plot <- 5
  num_plots <- ceiling(length(ptids)/ptids_per_plot)
  parentGate <- sub(".*\\/([^\\/]+)\\/[^\\/]+$", "\\1", currentGates[[1]])
  for(plot_num in seq_along(1:num_plots)) { 
    current_ptids <- ptids[((plot_num-1)*ptids_per_plot + 1):min(((plot_num-1)*ptids_per_plot + ptids_per_plot), length(ptids))]
    png(filename = file.path(here::here("out/QC/FACS_Plots"),
                             sprintf("%s%s_B3_pt%s_%s.png",
                                     if(parentGate %in% c("NOT4+", "4+", "8+")) {paste0(parentGate, "_")} else {""},
                                     sub(".*\\/([^\\/]+$)", "\\1", currentGates[[1]]), plot_num,
                                     format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
        width = 1920, height = 210+174*length(current_ptids), units = "px")
    print(plotter_mem(subset(gs_b3, `SAMPLE ID` %in% current_ptids),
                      currentGates) +
            labs(caption=sprintf("Batch 3, pt %s", plot_num)))
    dev.off()
  }
}

###########################################################

# Visualize activation gates on CD3-CD19+

# library(openCyto) # 1.24.0
# library(CytoML) # 1.12.0
# library(flowCore) # required for description()
# library(flowWorkspace) # required for gh_pop_get_data()
# library(ggcyto) # devtools::install_github("RGLab/ggcyto", ref="ggplot3.3") for update_theme()
# library(here)
# library(tidyverse)
# 
# gs_b3 <- load_gs(here::here("out/GatingSets/20200803_HAARVI_ICS_GatingSet_B3"))
# pData(gs_b3)$STIM <- factor(pData(gs_b3)$STIM, levels = c("DMSO", "VEMP", "Spike 1", "Spike 2", "NCAP", "SEB"))

# Activation gates on CD3-CD19+

# Since a CD3-CD19+ gate does not exist, looks like I will have to make one
# First copy singlet gate under Time
S_gate <- gs_pop_get_gate(gs_b3, "/Time/LD-3+/1419-3+/S")
LD_Neg_CD3_Neg_gate <- readRDS(here::here("processed_data/B1_LD_Neg_CD3_Neg_gate.rds")) # from batch 1

# Now try adding the singlet gate again, and then the LD-3- gate
gs_pop_add(gs_b3, S_gate, parent = "/Time")
gs_pop_add(gs_b3, LD_Neg_CD3_Neg_gate, parent = "/Time/S")

plot(gs_b3, bool=T, fontsize=15)

recompute(gs_b3, y="/Time/S")

# # Good enough gate (doesn't need to be accurate for this visualization)
# ggcyto(gs_b3[1], aes("<G610-A>", "<V510-A>"),
#        subset = "/Time/S",
#        filter = marginalFilter) +
#   geom_hex(bins=128) +
#   geom_gate("/Time/S/LD-3-")
# 
# # <V610-A>      CD38 BV605   <U395-A>    HLADR BUV395
# ggcyto(gs_b3[1], aes("<V610-A>", "<U395-A>"),
#        subset = "/Time/S/LD-3-",
#        filter = marginalFilter) +
#   geom_hex(bins=128) +
#   geom_gate(c("/Time/LD-3+/1419-3+/S/Lymph/CD38+", "/Time/LD-3+/1419-3+/S/Lymph/HLADR+"))

currentGates <- c("/Time/LD-3+/1419-3+/S/Lymph/CD38+", "/Time/LD-3+/1419-3+/S/Lymph/HLADR+")
ptids <- unique(pData(gs_b3)$`SAMPLE ID`)
ptids_per_plot <- 5
num_plots <- ceiling(length(ptids)/ptids_per_plot)
parentGate <- "/Time/S/LD-3-"
for(plot_num in seq_along(1:num_plots)) { 
  current_ptids <- ptids[((plot_num-1)*ptids_per_plot + 1):min(((plot_num-1)*ptids_per_plot + ptids_per_plot), length(ptids))]
  png(filename = file.path(here::here("out/QC/FACS_Plots"),
                           sprintf("%s%s_B3_pt%s_%s.png",
                                   paste0("LD-3-", "_"),
                                   sub(".*\\/([^\\/]+$)", "\\1", currentGates[[1]]), plot_num,
                                   format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))),
      width = 1920, height = 210+174*length(current_ptids), units = "px")
  
  myGS <- subset(gs_b3, `SAMPLE ID` %in% current_ptids)
  myGates <- currentGates
  
  firstGate <- myGates[[1]]
  currentGateBoundaries <- attributes(gh_pop_get_gate(myGS[[1]], firstGate))$boundaries
  currentXaxis <- colnames(currentGateBoundaries)[[1]]
  currentYaxis <- colnames(currentGateBoundaries)[[2]]
  
  print(ggcyto(myGS,
         aes(!!currentXaxis, !!currentYaxis),
         subset = "/Time/S/LD-3-",
         filter = marginalFilter) +
    geom_hex(bins=128) +
    geom_gate(myGates) +
    axis_x_inverse_trans() + axis_y_inverse_trans() +
    ggcyto_par_set(limits = "instrument") +
    facet_grid(`SAMPLE ID` ~ STIM, switch = "y") +
    theme(#plot.title = element_blank(),
      legend.position = "none",
      strip.text.x = element_text(size = 14, margin = margin(0,0,0,0, "cm")),
      panel.grid.major = ggplot2::element_blank(),
      plot.title = element_text(size=18)) +
    theme_bw(base_size=28) +
    # geom_stats(size=8, alpha=0.4) +
    labs(title=paste0(myGates[[1]], " on /Time/S/LD-3-"), caption = sprintf("Batch 3, pt %s", plot_num)))

  dev.off()
}
