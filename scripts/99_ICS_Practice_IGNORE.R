##### ICS Practice IGNORE

library(openCyto) # 1.24.0
library(CytoML) # 1.12.0
library(flowCore) # required for description()
library(flowWorkspace) # required for gh_pop_get_data()
library(ggcyto) # devtools::install_github("RGLab/ggcyto", ref="ggplot3.3") for update_theme()
library(here)
library(tidyverse)
library(readxl)


##### FCS files for Batch 1  already should already be downloaded into /data/20200528_COVID_ICS-B1
##### These can be obtained from the GDrive folder 20200528_COVID_ICS-B1/20200528_COVD_ICS-B1

date <- 20210826

xml_path_b1 <- here::here("data/20200528_COVID_ICS-B1_CS2_KY.xml")     # Location of XML file
fcs_subfolder_b1 <- here::here("data/20200528_COVID_ICS-B1/")          # Location of .fcs files
ws_b1 <- open_flowjo_xml(xml_path_b1)                                  # Create the workspace

b1_sample_group <- fj_ws_get_sample_groups(ws_b1)                      # Extract data frame of sample groups
head(b1_sample_group)
table(b1_sample_group$groupName)
b1_samples <- fj_ws_get_samples(ws_b1)                                 # Extract data frame of samples
head(b1_samples)
names(fj_ws_get_keywords(ws_b1, 50))                                   # Over 200 keywords tracked!
keywords2import <- c("EXPERIMENT NAME",                                # Let's plan to keep only these keywords
                     "$DATE", 
                     "SAMPLE ID",                                      # PATIENT ID dropped from Malisa's code because one PTID was missing it. 
                     "STIM", 
                     "WELL ID", 
                     "PLATE NAME")

gs_b1 <- flowjo_to_gatingset(ws_b1,                                    # Load Gating Set
                             name=sampleGroup, 
                             keywords=keywords2import,
                             path=fcs_subfolder_b1, 
                             extend_val=-10000)

# Warning message:
# In flowjo_to_gatingset(ws_b1, name = sampleGroup, keywords = keywords2import,  :
# GatingSet contains different gating tree structures and must be cleaned before using it!

# This means that some samples have different gating trees.  This can be fixed in R but best to check this stuff in FlowJo first

gh_get_pop_paths(gs_b1[[1]])                          # Get the gating tree for the first sample in the gating set
pop_lists <- lapply(gs_b1, gh_get_pop_paths)          # Make a list of all the gating trees in the gating set
unique(pop_lists)                                     # 3 unique gating trees!

# It's quickly obvious there is a problem.  The second tree has an errant node at [31] and the third tree at [10]
# Which samples contain these errant nodes?

x <- as.vector(NULL)
for(i in 1:length(pop_lists)){
  x[i] <- "/Time/LD-3-" %in% pop_lists[[i]]
}
x                         # The 6th entry is TRUE
pData(gs_b1[[6]])         # 112580

x <- as.vector(NULL)
for(i in 1:length(pop_lists)){
  x[i] <- "/Time/LD-3+/1419-3+/S/Lymph/4+/<V710-A>, <U730-A> subset" %in% pop_lists[[i]]
}
which(x == TRUE)          # The 30th and 90th entries are TRUE
pData(gs_b1[[c(30,90)]])  # 112540, 112538

# Remove the errant nodes from the gating set
gs_pop_remove(gs_b1["112580.fcs_414334"], "/Time/LD-3-")
gs_pop_remove(gs_b1[c("112538.fcs_166800", "112540.fcs_242350")], "/Time/LD-3+/1419-3+/S/Lymph/4+/<V710-A>, <U730-A> subset")

# Recheck the gating trees
pop_lists <- lapply(gs_b1, gh_get_pop_paths)
unique(pop_lists)         # Now only 1!
plot(gs_b1)               # plot the gating tree - some phantom nodes remaining?

# Save Gating Set (Create directory below only if you haven't done it already)
# dir.create(here::here("out/GatingSets"), recursive = T)
save_gs(gs_b1, here::here("out/GatingSets/20210826_HAARVI_ICS_GatingSet_B1"))

# End Code
#
#
