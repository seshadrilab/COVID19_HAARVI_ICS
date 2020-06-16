#!/usr/bin/Rscript

# 1) Ideally, run this script from the command line, i.e. NOT in RStudio. The rest of the pipeline can and should be run in RStudio.
# One way to do this is to navigate to the top level of the project folder, start up R in interactive mode, and paste in this code.
# Or run: R CMD scripts/20200612_RunCompass.R
# If you are on Unix/macOS, multicore processing (efficient forked/shared memory) only works when run from the command line.
# Windows does not support forking.
# See: future::supportsMulticore()
# 
# 2) It's also a good idea to log the output: type "script COMPASS_log.txt" in your terminal prior to running R.
# The COMPASS stdout and stderr should then get printed into COMPASS_log at the end of the furrr::future_pmap loop.
# Then after you quit R, type "exit" to close the log.

library(furrr) # For parallel computation 
library(COMPASS)
library(grid)
library(flowWorkspace)
source(here::here("scripts/20200604_Helper_Functions.R"))

date <- 20200612

stims_for_compass_runs <- c("VEMP", "Spike 1", "Spike 2", "NCAP") #, "SEB") # The SEB runs will take the longest, skip for now
parent_nodes_for_compass_runs <- c("4+", "NOT4+") # simpler to use shortened path
seeds_for_compass_runs <- as.list(date:(date + length(stims_for_compass_runs)*length(parent_nodes_for_compass_runs) - 1))

stims_for_compass_runs_rep <- rep(stims_for_compass_runs, each = length(parent_nodes_for_compass_runs))
parent_nodes_for_compass_runs_rep <- rep(parent_nodes_for_compass_runs, times = length(stims_for_compass_runs))

gsPath <- here::here("out/GatingSets/20200611_HAARVI_ICS_GatingSet_AllBatches/")
gs <- load_gs(gsPath)
gs2 <- subset(gs, !(`SAMPLE ID` %in% c("37C", "BWT23", "116C", "BWT22")) &
                !(`SAMPLE ID` == "07B" & STIM == "Spike 2"))

# mapMarkers contains output of markernames(gs)
mapMarkers <- list("IL2 PE", "IL4/5/13 APC", "IFNg V450", "TNFa FITC", "IL17a Ax700", "CD154 PE-Cy5", "CD107a PE-Cy7")
# TODO: change markernames to shortened cytokine version

future::supportsMulticore() # Run in terminal to get TRUE
# If you run this script in RStudio, the next line throws the following warning:
# "Warning message:
# [ONE-TIME WARNING] Forked processing ('multicore') is disabled in future (>= 1.13.0) when running R from RStudio,
# because it is considered unstable. Because of this, plan("multicore") will fall back to plan("sequential"),
# and plan("multiprocess") will fall back to plan("multisession") - not plan("multicore") as in the past.
# For more details, how to control forked processing or not, and how to silence this warning in future R sessions,
# see ?future::supportsMulticore "
future::plan(multiprocess(workers = max(1, availableCores() - 2)))

out <- furrr::future_pmap(.l = list(stims_for_compass_runs_rep,
                              parent_nodes_for_compass_runs_rep,
                              seeds_for_compass_runs),
                    .f = function(currentStim, parent, currentSeed) {
                      
                        o <- tryCatch( {
                          gsSub <- subset(gs2, STIM %in% c("DMSO", currentStim))
                          currentNodeMarkerMap <- mapMarkers
                          # currentNodeMarkerMap names are gating tree paths
                          names(currentNodeMarkerMap) <- paste0(parent, "/", c("IL2", "IL4513", "IFNG", "TNF", "IL17", "154", "107a"))
                          outDir <- here::here(sprintf("out/CompassOutput/%s/%s", parent, gsub(" ", "_", currentStim)))
                          if(!dir.exists(outDir)) {
                            dir.create(outDir, recursive = T)
                          }
                          
                          runCompassOnce(gs=gsSub,
                                         seed=currentSeed,
                                         outDir=outDir,
                                         parentNode=parent,
                                         nodeMarkerMap=currentNodeMarkerMap,
                                         uniqueIdentifier="SAMPLE ID",
                                         treatmentCol="STIM",
                                         currentTreatment=currentStim,
                                         currentControl="DMSO",
                                         stratifyBy=NULL, # "Cohort"
                                         iter=40000,
                                         eventCountFilterThreshold=0,
                                         textForRunOutputId=paste0(parent, "_", gsub(" ", "_", currentStim)))
                          gc()
                        }, error = function(e) { print(e) })
                        o
                    },
                    # Progress bar reflects how many COMPASS runs have completed
                    .progress = T)