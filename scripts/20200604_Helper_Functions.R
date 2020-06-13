# t-SNE Helper functions

# Distribute totalEvents as evenly as possible among a vector of different sub-group sizes
# Recursive function
#
# Example:
# cd1c_panel_counts %>%
#   dplyr::select(NHP, TimepointTissue, Batch, rCD1c_GMM_clean) %>% 
#   group_by(TimepointTissue, Batch) %>% 
#   nest() %>% 
#   ungroup() %>% 
#   mutate(nsamp = map2(376, data, function(totalEvents, df) {
#     distributeEvents(totalEvents, df$rCD1c_GMM_clean)
#   })) %>% 
#   unnest(cols = c(data, nsamp))
distributeEvents <- function(totalEvents, subGroupSizes) {
  #print(totalEvents)
  #print(subGroupSizes)
  if(totalEvents == 0) {
    rep(0, length(subGroupSizes))
  } else if(sum(subGroupSizes) < totalEvents) {
    stop("Not enough events to achieve requested sample size")
  } else if(any(subGroupSizes == 0)) {
    # This branch exists in case the first call to this function includes some subGroups of size 0. 
    output <- rep(0, length(subGroupSizes))
    nonZeroIndices <- which(subGroupSizes > 0)
    output[nonZeroIndices] <- distributeEvents(totalEvents, subGroupSizes[nonZeroIndices])
  } else {
    numSubGroups <- length(subGroupSizes)
    total_div_numSubGroups <- totalEvents %/% numSubGroups # ideally, how many events would be sampled from each sub-group?
    remainder <- totalEvents %% numSubGroups
    minGroupSize <- min(subGroupSizes)
    
    if(total_div_numSubGroups == 0) {
      # totalEvents is non-zero but total_div_numSubGroups is zero, meaning there is a non-zero remainder
      stopifnot(remainder != 0)
      # If the function got to this point, we know that all subGroups are non-zero (see previous branch)
      stopifnot(all(subGroupSizes > 0))
      
      # Ok, now time to distribute the remainder among subGroupSizes events
      # We know that remainder is smaller than the number of subGroups
      stopifnot(remainder < numSubGroups)
      # So now we have to evenly and randomly distribute the remainder among the subGroups
      # This is the only random portion of this function
      output <- rep(0, numSubGroups)
      output[sample.int(numSubGroups, remainder)] <- 1
      # This is a base case, so we return output. No more recursive calls.
      output
    } else {
      howMuchToAssign <- min(total_div_numSubGroups, minGroupSize)
      output <- rep(howMuchToAssign, numSubGroups)
      
      #print(output)
      
      remaining_subGroupSizes <- subGroupSizes - howMuchToAssign
      nonZeroIndices <- which(remaining_subGroupSizes > 0)
      if(length(nonZeroIndices)) {
        toAdd <- rep(0, numSubGroups)
        toAdd[nonZeroIndices] <- distributeEvents(totalEvents - sum(output), remaining_subGroupSizes[nonZeroIndices])
        output + toAdd
      } else {
        output
      }
    }
  }
}

# A function to sample n events from a GatingHierarchy containing a single sample
sampleGatingHierarchy <- function(gh, parentGate, n, otherGates = NULL) {
  library(openCyto)
  library(CytoML) # 1.12.0
  library(flowCore) # required for description()
  library(flowWorkspace) # required for gh_pop_get_data()
  
  stopifnot(length(gh) == 1)
  allMarkerNames <- pData(parameters(gh_pop_get_data(gh)))[,c(1,2)] # First column is flow channel, second is marker name
  if (any(is.na(allMarkerNames[,2])) | length(unique(allMarkerNames[,2])) < length(allMarkerNames[,2]))
    stop ("all marker names (even FSC-A and Time) must be assigned and be unique")
  # May want to loosen above requirement, using channel names where marker names are unavailable.
  
  # First take length n sample from all events in parentGate
  availableEvents <- gh_pop_get_count(gh, parentGate)
  nSampled <- sample.int(availableEvents, size = n)
  
  # Then extract boolean gate data and mfi data, and cbind them.
  # For now, in the interest of saving memory, the data is subset to the sampled events as soon as possible.
  parentGateIndices <- gh_pop_get_indices(gh, parentGate) # relative to all data in gh, i.e. root node. a TRUE/FALSE vector. Used for subsetting boolean data
  gates2Extract <- unique(c(parentGate, if(is.null(otherGates)) {gh_get_pop_paths(gh)} else {otherGates}))
  perCellGateMembership <- data.frame(lapply(gates2Extract, function(currentGate) {
    as.integer(gh_pop_get_indices_mat(gh, currentGate)[parentGateIndices,][nSampled]) }))
  colnames(perCellGateMembership) <- gates2Extract
  
  perCellMFIData <- exprs(gh_pop_get_data(gh, parentGate))[nSampled,]
  colnames(perCellMFIData) <- allMarkerNames[match(colnames(perCellMFIData), allMarkerNames[,1]), 2]
  
  # Combine Gate membership data with MFI expression data
  stopifnot(nrow(perCellGateMembership) == nrow(perCellMFIData))
  gateAndMFIData <- cbind(perCellGateMembership, perCellMFIData)
  
  # Add metadata and return 
  cbind(pData(gh), gateAndMFIData, row.names = NULL)
}

# COMPASS Helper functions

# TODO: load libraries. COMPASS, grid, ?
runCompassOnce <- function(gs,
                           seed=NULL,
                           outDir,
                           parentNode,
                           nodeMarkerMap,
                           uniqueIdentifier,
                           treatmentCol="trt",
                           currentTreatment="Treatment",
                           currentControl="Control",
                           stratifyBy=NULL,
                           iter=40000,
                           eventCountFilterThreshold=0,
                           textForRunOutputId=NULL) {
  require(COMPASS)
  require(grid)
  
  currentRunTextForConsole <- paste0(parentNode, " ", currentTreatment)
  currentRunTextForConsole <- if(is.null(textForRunOutputId)) {
    currentRunTextForConsole
  } else {
    sprintf("%s (%s)", currentRunTextForConsole, textForRunOutputId)
  }
  message(sprintf("Running COMPASS for %s", currentRunTextForConsole))
  
  # Set the seed
  if (!is.null(seed)) {
    rngKind <- "L'Ecuyer-CMRG" # This seems to be the recommended kind of RNG for reproducible parallel processing
    RNGkind(kind = rngKind)
    message(sprintf("Setting %s seed to %s", rngKind, seed))
    set.seed(seed)
  }
  
  # Create a COMPASSContainer from the GatingSet or GatingSetList.
  CC <- COMPASS::COMPASSContainerFromGatingSet(gs, node=parentNode, individual_id=uniqueIdentifier,
                                               mp=nodeMarkerMap, countFilterThreshold=eventCountFilterThreshold)
  
  # Run COMPASS
  # The treatment and control arguments for COMPASS::COMPASS() accept expressions, which makes it difficult to use COMPASS::COMPASS() programmatically.
  # This work-around uses bquote() to build up an expression consisting of a call to COMPASS::COMPASS(). eval() then evaluates this expression.
  # bquote() returns an expression consisting of the code you provided, but it first evaluates everything in .() 
  # It is not an ideal solution, but it works for now.
  fit <- eval(bquote(COMPASS::COMPASS( CC,
                                       treatment= .(as.name(treatmentCol)) == .(eval(currentTreatment)), 
                                       control= .(as.name(treatmentCol)) == .(eval(currentControl)),
                                       iterations=iter
  )))
  
  message("COMPASS complete, now saving output")
  
  if(is.null(textForRunOutputId)) {
    textForRunOutputId <- paste0(parentNode, "_", currentTreatment)
  }
  
  # Save the COMPASS run output as an RDS file for safekeeping
  saveRDS(fit, file.path(outDir, sprintf("COMPASSResult_%s.rds", textForRunOutputId)))
  
  # Save the Functionality and Polyfunctionality Scores to a tsv
  # TODO: rewrite this using stack(COMPASS::FunctionalityScore(fit)) instead of data.frame()
  FS <- COMPASS::FunctionalityScore(fit)
  PFS <- COMPASS::PolyfunctionalityScore(fit)
  FS_df <- data.frame(tmp = names(FS), FS = FS)
  colnames(FS_df) <- c(uniqueIdentifier, "FS")
  PFS_df <- data.frame(tmp = names(PFS), PFS = PFS)
  colnames(PFS_df) <- c(uniqueIdentifier, "PFS")
  FS_PFS_df <- merge(FS_df, PFS_df, by = uniqueIdentifier)
  write.table(FS_PFS_df,
              file = file.path(outDir, sprintf("FS_PFS_%s.tsv", textForRunOutputId)),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # Plot a heatmap of the mean probability of response
  plotTitleSuffix <- paste(c(",\n", run, ", ", parentNode, " Cells"), collapse="")
  cytokine_annotation_colors <- c("black", "black", "black", "black", "black", "black", "black")
  
  png(filename=file.path(outDir, sprintf("HeatmapMeanProbResponse_%s.png", textForRunOutputId)),
      width=800, height=650)
  # COMPASS::plot.COMPASSResult implements the generic S3 method graphics::plot(), which uses the object's class (in this case a COMPASSResult) to dispatch to the COMPASS plotting function
  # The COMPASS package doesn't allow me to call COMPASS::plot.COMPASSResult directly
  try(grid::grid.draw(print(graphics::plot(fit, stratifyBy, show_rownames = TRUE,
                                           main = sprintf("Heatmap of Mean Probability of Response %s", textForRunOutputId),
                                           fontsize=14, fontsize_row=13, fontsize_col=11,
                                           cytokine_annotation_colors=cytokine_annotation_colors))))
  dev.off()
  
  message(sprintf("Done with run %s", currentRunTextForConsole))
}

############################################################################################################################

# Dotplot of background-corrected subset percents out of the parent population
# cr is the COMPASSResult object
make_dotplot_for_COMPASS_run <- function(cr, run_name, output_folder=NA, current_ylim=NULL, add_legend=FALSE, legend_position=c(0.16, 0),
                                         save_test_results=TRUE, p_text_size=5, include_0_line=FALSE, zeroed_BgCorr = FALSE, plot_width=7, plot_height=6,
                                         dichotomize_by_cytokine=NA, group_by_colname="Status", group_by_order=c("Persistent_Neg", "Qfn_Pos"),
                                         group_by_colors=c("Persistent_Neg" = "#373db8", "Qfn_Pos" = "#28c914"), return_output=FALSE, parentSubset="CD4+",
                                         point_size=0.3) {
  
  mean_gamma <- cr$fit$mean_gamma
  cats <- as.data.frame(cr$fit$categories[, -ncol(cr$fit$categories),drop=FALSE]) # drop the "Counts" column
  numMarkers <- ncol(cats)
  rownames(cats) <- colnames(mean_gamma)
  
  # Filter the subsets to those where the average mean_gamma is greater than the threshold (default in heatmap and this function is 0.01)
  compassSubsetsFiltered <- names(which(apply(mean_gamma, 2, function(x) { mean(x, na.rm = TRUE) }) > 0.01))
  # Note that we don't need mean_gamma anymore for this task
  # And remove the subset with 0 positive markers
  compassSubsetsFiltered <- compassSubsetsFiltered[lengths(regmatches(compassSubsetsFiltered, gregexpr("!", compassSubsetsFiltered))) != numMarkers]
  
  # Subset the cats rows to compassSubsetsFiltered, and put the columns in the order of cytokine_order_for_annotation
  cats <- cats[compassSubsetsFiltered, cytokine_order_for_annotation]
  
  stim_counts <- as.data.frame(cr$data$n_s) %>%
    mutate(Individual = rownames(cr$data$n_s)) %>%
    dplyr::select(c("Individual", compassSubsetsFiltered)) %>% 
    dplyr::left_join(cr$data$counts_s %>%
                       stack() %>%
                       rename("ParentCount" = "values", "Individual" = "ind"),
                     by = "Individual") %>%
    dplyr::left_join(cr$data$meta %>% 
                       dplyr::select(!!as.symbol(cr$data$individual_id), !!as.symbol(group_by_colname)),
                     by = c("Individual"=cr$data$individual_id)) %>% 
    mutate(Stim = "Dummy_Stim_Name")
  bg_counts <- as.data.frame(cr$data$n_u) %>%
    mutate(Individual = rownames(cr$data$n_u)) %>%
    dplyr::select(c("Individual", compassSubsetsFiltered)) %>% 
    dplyr::left_join(cr$data$counts_u %>%
                       stack() %>%
                       rename("ParentCount" = "values", "Individual" = "ind"),
                     by = "Individual") %>%
    dplyr::left_join(cr$data$meta %>% 
                       dplyr::select(!!as.symbol(cr$data$individual_id), !!as.symbol(group_by_colname)),
                     by = c("Individual"=cr$data$individual_id)) %>% 
    mutate(Stim = "DMSO")
  dat_bgCorr_long <- bind_rows(bg_counts, stim_counts) %>% 
    mutate_at(.vars = compassSubsetsFiltered, `/`, quote(ParentCount)) %>%  # convert counts to proportions
    dplyr::select(-ParentCount) %>% 
    tidyr::pivot_longer(cols = -c("Individual", "Stim", !!as.symbol(group_by_colname)),
                        names_to = "BooleanSubset",
                        values_to = "Proportion") %>% 
    tidyr::pivot_wider(id_cols = c("Individual", !!as.symbol(group_by_colname), "BooleanSubset"),
                       names_from = Stim,
                       values_from = Proportion) %>% 
    drop_na() %>% # Filter out rows with NA
    mutate(BgCorr = if(zeroed_BgCorr) {pmax(0, Dummy_Stim_Name - DMSO)} else {Dummy_Stim_Name - DMSO}) %>% 
    dplyr::select(-c(DMSO, Dummy_Stim_Name))
  dat_bgCorr_wide <- dat_bgCorr_long %>% 
    tidyr::pivot_wider(id_cols = c("Individual", !!as.symbol(group_by_colname)),
                       names_from = BooleanSubset,
                       values_from = BgCorr)
  
  tests <- lapply(compassSubsetsFiltered, function(boolSubset) {
    wilcox.test(as.formula(sprintf("`%s` ~ %s", boolSubset, group_by_colname)), data=dat_bgCorr_wide)
  })
  pvals_df <- data.frame(BooleanSubset = compassSubsetsFiltered,
                         p = unlist(lapply(tests, function(x) {x$p.value}))) %>% 
    mutate(p.adj = p.adjust(p, method = "bonferroni")) %>% # Strict
    mutate(p.adj.text = if_else(p.adj < 0.001, "p<.001", paste0("p=", sub("0.", ".", round(p.adj, 3)))))
  
  # Before plotting, put the categories data frame rows in the desired order (columns were already re-ordered above)
  # This is essentially the same code I added to the pheatmap function
  cats <- cats[rev(do.call(order, cats)),,drop=FALSE]
  # And order the cats df rows by degrees (number of cytokines in subset)
  ckr<-apply(cats,1,function(x)sum(as.numeric(as.character(x))))
  cats = cats[order(ckr),]
  if(!is.na(dichotomize_by_cytokine)) {
    # And then dichotomize the cats df rows so that all subsets containing the cytokine in dichotomize_by_cytokine (e.g. "IFNg") appear last
    cats <- cats[order(cats[,dichotomize_by_cytokine]),]
  }
  
  # Use the cats df row order to order the boolean subsets in dat_bgCorr_long and pvals_df
  dat_bgCorr_long$BooleanSubset <- factor(dat_bgCorr_long$BooleanSubset, levels = rownames(cats))
  pvals_df$BooleanSubset <- factor(pvals_df$BooleanSubset, levels = rownames(cats))
  
  # Calculate medians of each group for each subset
  dat_bgCorr_medians <- dat_bgCorr_long %>%
    dplyr::group_by(!!as.symbol(group_by_colname), BooleanSubset) %>%
    dplyr::summarise(BgCorr = median(BgCorr))
  
  if(save_test_results) {
    # Save some form of the test results to disk so it doesn't only exist in the plot as adjusted p-values
    pvals_df_for_file <- cats %>% rownames_to_column("BooleanSubset") %>% 
      dplyr::left_join(pvals_df) %>% 
      dplyr::left_join(dat_bgCorr_medians %>%
                         pivot_wider(id_cols = "BooleanSubset",
                                     names_from = !!as.symbol(group_by_colname),
                                     values_from = BgCorr,
                                     names_prefix = "med_")) %>% 
      arrange(p.adj)
    
    if(!is.na(output_folder)) {
      test_results_file_path <- file.path(output_folder,
                                          sprintf("%s_%s_BooleanSubsets_BgCorrProps_MannWhitney.tsv",
                                                  run_name, if(zeroed_BgCorr) {"Zeroed"} else {"NotZeroed"}))
      write.table(pvals_df_for_file,
                  file = test_results_file_path,
                  sep = "\t", row.names = FALSE, quote = FALSE)
    }
  }
  
  # Draw the dotplot
  p_dotplot <- ggplot(dat_bgCorr_long, aes(x = !!as.symbol(group_by_colname), y = BgCorr, fill = !!as.symbol(group_by_colname), group = !!as.symbol(group_by_colname)))
  if(include_0_line) {
    p_dotplot <- p_dotplot + geom_hline(yintercept = 0, linetype="dashed", alpha = 0.5)
  }
  p_dotplot <- p_dotplot +
    geom_jitter(aes(color = !!as.symbol(group_by_colname)), width = 0.2, size=point_size) +
    geom_errorbarh(data = dat_bgCorr_medians,
                   aes(y = BgCorr,
                       xmax = 1.5 + 0.6,
                       xmin = 1.5 - 0.6, height = 0),
                   position=position_dodge(width=0.25), color = "black") +
    facet_grid(. ~ BooleanSubset) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text = element_text(color="black", size=12),
          axis.title = element_text(size=18),
          text = element_text(family="Arial"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) + 
    labs(y=sprintf("%% Responding %s T cells", sub("+", "", parentSubset, fixed=T)))
  if(!is.na(group_by_colors)) {
    p_dotplot <- p_dotplot + scale_color_manual(values=group_by_colors)
  }
  
  # Adjust ylim here manually if specified
  if(is.null(current_ylim)) {
    p_dotplot <- p_dotplot +
      scale_y_continuous(labels = function(x) paste0(x*100))
  } else {
    p_dotplot <- p_dotplot +
      scale_y_continuous(labels = function(x) paste0(x*100), limits=current_ylim)
  }
  
  if(add_legend) {
    p_dotplot <- p_dotplot +
      theme(legend.justification = c(1,1),
            legend.position = legend_position,
            legend.text=element_text(size=16)) +
      guides(color=guide_legend(title=NULL), fill = FALSE)
  } else {
    p_dotplot <- p_dotplot +
      theme(legend.position = "none")
  }
  
  showSignificanceBracket <- TRUE
  p_alpha <- 0.05
  onlyShowPBelowAlpha <- TRUE
  if(showSignificanceBracket) {
    
    # TODO rewrite this
    get_y_pos <- function(boolSubsets) {
      sapply(boolSubsets, function(boolSubset) {
        boolSubset <- as.character(boolSubset)
        find_y_max_in_visible_range <- function(x) {ifelse(is.null(current_ylim), max(x), max(subset(x, x < current_ylim[[2]])))}
        y_visible_max <- dat_bgCorr_wide %>%
          group_by(!!as.symbol(group_by_colname)) %>%
          summarise(y_visible_max = find_y_max_in_visible_range(!!as.symbol(boolSubset))) %>% 
          dplyr::pull(y_visible_max) %>% 
          max()
        y_visible_max + if(is.null(current_ylim)) {max(dat_bgCorr_long$BgCorr)/20} else {current_ylim[[2]]/20}
      })
      
    }
    
    annotation_df <- pvals_df %>% 
      mutate(start = group_by_order[[1]],
             end = group_by_order[[2]],
             y_pos = get_y_pos(BooleanSubset))
    if(onlyShowPBelowAlpha) {
      annotation_df <- annotation_df %>% dplyr::filter(p.adj < p_alpha)
    }
    
    # If I don't use the full path for ggsignif::geom_signif, it may try to use a global environment variable GeomSignif and ignore manual = T. Odd.
    p_dotplot <- p_dotplot +
      ggsignif::geom_signif(inherit.aes=F,data=annotation_df,
                            aes_string(xmin="start", xmax="end", annotations="p.adj.text", y_position="y_pos"), # , family="Arial"
                            tip_length = c(0.005, 0.005),
                            textsize=p_text_size,
                            manual = TRUE)
  }
  
  # Now make the categories legend
  
  # Set the order of the cytokines and subsets once the categories df is in long format. Then plot.
  cats_long <- as.data.frame(cats) %>%
    rownames_to_column("BooleanSubset") %>%
    gather(Cytokine, Membership, -BooleanSubset) %>%
    mutate(Membership = dplyr::recode(Membership, "0" = "-", "1" = "+")) %>% 
    mutate(BooleanSubset = factor(BooleanSubset, levels = rownames(cats)),
           Cytokine = factor(Cytokine, levels = colnames(cats)))
  
  cats_plot <- ggplot(cats_long,
                      aes(x = BooleanSubset, y = Cytokine)) +
    geom_tile(fill="white") +
    geom_text(aes(label=Membership), color="black", size=7) +
    theme(axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title = element_blank(),
          axis.text.y = element_text(color="black", size=14),
          panel.border=element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), "cm")) + 
    scale_y_discrete(expand=c(0,0)) + 
    scale_x_discrete(expand=c(0,0))
  
  dotplot_with_cats <- plot_grid(p_dotplot, cats_plot, ncol = 1, axis = "lr", align = "v", rel_heights = c(1, 0.4))
  
  if(!is.na(output_folder)) {
    p_base_path <- file.path(output_folder, sprintf("%s_%s_Dotplot%s%s", run_name,
                                                    if(zeroed_BgCorr) {"Zeroed"} else {"NotZeroed"},
                                                    if(!is.null(current_ylim)) { paste0("_ylim_", paste0(paste0(c("min", "max"), current_ylim), collapse="")) } else {""},
                                                    if(include_0_line) { "_with0line" } else { "" }))
    ggsave(filename=paste0(p_base_path, ".png"), plot=dotplot_with_cats, width=plot_width, height=plot_height, dpi=300)
    ggsave(filename=paste0(p_base_path, ".pdf"), plot=dotplot_with_cats, width=plot_width, height=plot_height, units = "in",
           onefile = TRUE, bg = "transparent", family = "Arial", fonts = "Arial")
  }
  
  if(return_output) {
    to_return <- list("Dotplot" = dotplot_with_cats)
    if(save_test_results) {
      to_return$Test_Results <- pvals_df_for_file
    }
    to_return
  }
}
