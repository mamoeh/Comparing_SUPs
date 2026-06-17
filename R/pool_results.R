#source("R/remove_bogus.R")
#source("R/check_primary.R")

##### Pool results over many instances ####################################### #
# filenames: names of RDS files (including '.rds') to pool
# path:      path to folder where RDS files are stored
############################################################################## #

pool_results <- function(filenames, path = "merged_data") {
 
  n_inst <- length(filenames)
  
  # load RDS's
  inst_list <- vector("list", n_inst)
  for(i in 1:n_inst) { 
    inst_list[[i]] <- readRDS(paste0(path, "/", filenames[i])) 
  }
  
  inst_names <- substr(filenames, 1, nchar(filenames) - 4)
  names(inst_list) <- inst_names
  
  # do the pooling
  inst_pooled <- pool_RDS(inst_list)
  
  
  ## collect new info
  
  # info from info_output()
  info_out <- matrix(nrow = n_inst, ncol = 12)
  for(i in 1:n_inst) {
    
    info <- info_output(inst_list[[i]]$df_merged, 
                        lapply(inst_list[[i]]$hierarchies, remove_bogus),
                        as_char = FALSE)
    
    info_out[i, 1:7]  <- info
    info_out[i, 8:12] <- round(info[2:6] / info[1], 3) * 100
    
    if(i == 1) { 
      nam <- names(info)
      colnames(info_out) <- c(nam, paste0("perc_", substr(nam[2:6], 3, nchar(nam[2:6]))))
    }
  }
  info_out <- as.data.frame(info_out)
  info_out$instance <- inst_names
  
  # share of inner cells
  df_inner <- inst_pooled$df_merged |> 
    dplyr::group_by(instance) |>
    dplyr::summarise(n_inner = sum(inner))
  info_out$perc_inner <- round(df_inner$n_inner / info_out$n_output, 3) * 100
  
  # table dimensionality
  info_out$dim <- sapply(lapply(inst_list, `[[`, "hierarchies"), length)
  
  # max. HiTaS
  df_hitas <- inst_pooled$Suppressions_by_Class |>
    group_by(instance) |>
    summarise(max_HiTaS = max(HiTaS_Class))
  info_out$max_HiTaS <- df_hitas$max_HiTaS
  
  # completeness
  df_complete <- inst_pooled$analysis |>
    dplyr::mutate(method_run = no_secondaries != 0) |>
    dplyr::group_by(instance) |>
    dplyr::summarise(methods_run = sum(method_run)) |>
    dplyr::mutate(complete = methods_run == 4)
  info_out[, c("methods_run", "complete")] <- df_complete[, c("methods_run", "complete")]
  
  inst_pooled$info_instances <- info_out

  inst_pooled
}


pool_RDS <- function(inst_list) {
  
  n_dfm <- length(inst_list)
  
  instance_pooled <- vector("list", 0)
  
  # extract $df_merged and $analysis
  l_dfm   <- lapply(inst_list, `[[`, "df_merged")
  l_other <- lapply(inst_list, `[[`, "analysis")
  l_other <- lapply(l_other, analysis_to_df)
  l_ana <- lapply(l_other, `[[`, 1)
  l_hit <- lapply(l_other, `[[`, 2)
  
  # instance names (if included)
  if(is.null(names(inst_list))) {
    nam <- 1:n_dfm
  } else {
    nam <- names(inst_list)
  }
  # label analysis results with respective instance name
  for(i in 1:n_dfm){ 
    l_dfm[[i]]$instance <- nam[i] 
    l_ana[[i]]$instance <- nam[i]
    l_hit[[i]]$instance <- nam[i]
  }
  
  instance_pooled$df_merged <- dplyr::bind_rows(l_dfm)
  instance_pooled$analysis <- dplyr::bind_rows(l_ana)
  instance_pooled$Suppressions_by_Class <- dplyr::bind_rows(l_hit)

  instance_pooled
}


analysis_to_df <- function(ana) {
  
  df_ana <- as.data.frame(ana[1:5])
  df_ana$method <- c("gauss", "modular", "simpleheuristic", "simpleheuristic_old")
  
  list(df_ana, ana$Suppressions_by_Class)
}

