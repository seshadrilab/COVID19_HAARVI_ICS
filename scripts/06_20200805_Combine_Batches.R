library(openCyto) # 1.24.0
library(CytoML) # 1.12.0
library(flowCore) # required for description()
library(flowWorkspace) # required for gh_pop_get_data()
library(here)
library(tidyverse)

# Note: This script requires a lot of space in the R temporary directory.
# If you need to, you can modify the directory location, e.g.:
# Open /etc/R/Renviron and add a line TMP = 'path/to/large/space'
tempdir() # And make sure it's where you want it to be

gs1 <- load_gs(here::here("out/GatingSets/20200803_HAARVI_ICS_GatingSet_B1/"))
gs2 <- load_gs(here::here("out/GatingSets/20200803_HAARVI_ICS_GatingSet_B2/"))
gs3 <- load_gs(here::here("out/GatingSets/20200803_HAARVI_ICS_GatingSet_B3"))

# Make sure nodes and markers are consistent between the three batches
setdiff(sort(gh_get_pop_paths(gs1)), sort(gh_get_pop_paths(gs2)))
setdiff(sort(gh_get_pop_paths(gs1)), sort(gh_get_pop_paths(gs3)))
all(sort(gh_get_pop_paths(gs1)) == sort(gh_get_pop_paths(gs2)))
all(sort(gh_get_pop_paths(gs1)) == sort(gh_get_pop_paths(gs3)))
all(colnames(pData(gs1)) == colnames(pData(gs2)))
all(colnames(pData(gs1)) == colnames(pData(gs3)))
all(markernames(gs1) == markernames(gs2))
all(markernames(gs1) == markernames(gs3))

pData(parameters(gh_pop_get_data(gs1[[1]])))[,c(1, 2)]

# Combining into a single GatingSet necessary if running t-SNE (we will be running both COMPASS and t-SNE on this dataset)
gs <- gslist_to_gs(GatingSetList(list(gs1, gs2, gs3))) # This takes a while

save_gs(gs, here::here("out/GatingSets/20200805_HAARVI_ICS_GatingSet_AllBatches"))
# gs <- load_gs(here::here("out/GatingSets/20200805_HAARVI_ICS_GatingSet_AllBatches"))

#######################################################

gs_sub <- subset(gs, STIM == "DMSO" &
                   !(`SAMPLE ID` %in% c("37C", "BWT23", "116C", "BWT22")))

# Prior to running COMPASS, plot the DMSO marginal frequencies for the 7 cytokine gates.
# We are looking for differences between hospitalized and not hospitalized that might confound COMPASS.
dput(grep("\\/4\\+\\/|\\/NOT4\\+\\/", gh_get_pop_paths(gs[[1]]), value=T))
pop_paths_of_interest <- c("/Time/LD-3+/1419-3+/S/Lymph/4+", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+",
                           "/Time/LD-3+/1419-3+/S/Lymph/8+",
                           "/Time/LD-3+/1419-3+/S/Lymph/4+/107a", "/Time/LD-3+/1419-3+/S/Lymph/4+/154", 
                           "/Time/LD-3+/1419-3+/S/Lymph/4+/IFNG", "/Time/LD-3+/1419-3+/S/Lymph/4+/IL2", 
                           "/Time/LD-3+/1419-3+/S/Lymph/4+/IL17", "/Time/LD-3+/1419-3+/S/Lymph/4+/IL4513", 
                           "/Time/LD-3+/1419-3+/S/Lymph/4+/TNF", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/107a", 
                           "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/154", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IFNG", 
                           "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL2", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL17", 
                           "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/IL4513", "/Time/LD-3+/1419-3+/S/Lymph/NOT4+/TNF",
                           "/Time/LD-3+/1419-3+/S/Lymph/8+/107a", "/Time/LD-3+/1419-3+/S/Lymph/8+/154", 
                           "/Time/LD-3+/1419-3+/S/Lymph/8+/IFNG", "/Time/LD-3+/1419-3+/S/Lymph/8+/IL2", 
                           "/Time/LD-3+/1419-3+/S/Lymph/8+/IL17", "/Time/LD-3+/1419-3+/S/Lymph/8+/IL4513", 
                           "/Time/LD-3+/1419-3+/S/Lymph/8+/TNF")
pops_of_interest_short <- str_replace(pop_paths_of_interest, "\\/Time\\/LD\\-3\\+\\/1419\\-3\\+\\/S\\/Lymph\\/", "")
pop_dat <- gs_pop_get_count_with_meta(gs_sub, 
                                      subpopulations = pop_paths_of_interest) %>% 
  dplyr::select(Population, Count, "SAMPLE ID", "Cohort", "Age", "Sex", "Race", "Hispanic?", 
                "Days symptom onset to visit 1", "Batch") %>% 
  pivot_wider(names_from = Population, values_from = Count) %>% 
  rename_at(vars(all_of(pop_paths_of_interest)),
            ~ pops_of_interest_short) %>% 
  mutate(Cohort = ifelse(is.na(Cohort) | Cohort == "Healthy control 2017-2018", "Healthy control", Cohort),
         Cohort = factor(Cohort, levels = c("Healthy control", "Non-hospitalized", "Hospitalized")))

pop_dat %>% dplyr::filter(Cohort == "Healthy control") %>% dplyr::pull(`SAMPLE ID`) # two TRIMAs remaining

plot_pop <- function(pop) {
  parent <- sub("(.*)\\/.*", "\\1", pop)
  tmp_dat <- pop_dat %>% 
    mutate(prop = !!as.name(pop) / !!as.name(parent))
  mw_p <- wilcox.test(prop ~ Cohort, data = tmp_dat %>% 
                        dplyr::filter(Cohort %in% c("Non-hospitalized", "Hospitalized")))$p.value
  p.unadj.text <- sprintf("Mann-Whitney U test (Non-hosp vs Hosp):   p-unadj%s",
                          if_else(mw_p < 0.001, "<.001", paste0("=", sub("0.", ".", round(mw_p, 3)))))
  ggplot(tmp_dat,
         aes(Cohort, prop)) +
    geom_boxplot(outlier.shape = NA) +
    theme_bw(base_size = 22) +
    geom_jitter(width=0.15, height=0, pch=21, fill="grey", alpha=0.8) +
    labs(y=sprintf("%% %s of %s", sub(".*\\/(.*)", "\\1", pop), parent),
         caption = paste0("DMSO cytokine frequencies\n", p.unadj.text)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_text(size=22),
          axis.text.y = element_text(color="black", size=22),
          axis.text.x = element_text(color="black", size=22),
          plot.title = element_blank(),
          plot.caption = element_text(size=12),
          panel.grid.minor = element_blank(),
          #legend.position = "none",
          plot.margin = margin(1.3, 0.2, 0, 0.2, "cm")) +
    scale_x_discrete(labels=c("Healthy control" = "Healthy", "Non-hospitalized" = "Conv\nNon-Hosp", "Hospitalized" = "Conv\nHosp")) +
    scale_y_continuous(labels = function(x) paste0(x*100))
}

plot_pop(pops_of_interest_short[[4]])

for(pop in pops_of_interest_short[4:length(pops_of_interest_short)]) {
  png(file=here::here(sprintf("out/QC/DMSO_Cytokine_Signal/20200805_%s_vs_Cohort.png", sub("\\/", "_", pop))),
      width=408, height=265, units = "px")
  print(plot_pop(pop))
  dev.off()
}

# NOT4+IL2+ and NOT4+IL4/5/13 are of note. Enriched in Conv-Hosp group. However, the significant tests would not pass multiple comparisons corrections.