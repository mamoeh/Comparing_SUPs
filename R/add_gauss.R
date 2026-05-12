


if (FALSE) {  # Example 
  hiers <- list(var1 = read_hier("minimal"), var2 = read_hier("table_header"))
  df <- create_microdata(hiers, n_ids =  100, n_unique = 20)
  initialize_gauss("test", df, hiers, overwrite = TRUE)
  all <- readRDS("merged/test.rds")
  head(all$df_merged)
  add_gauss("test")
  all <- readRDS("merged/test.rds")
  head(all$df_merged)
}

# To add suppressed_gauss to the file using precomputed primary_gauss.
#
# Use output = "df_merged" to return a data frame instead  
# 
add_gauss <- function(filename, 
                      path = "merged", 
                      output = NULL) {
  
  all <- readRDS(file.path(path, paste0(filename, ".rds")))
  
  hrc_GAUSS <- all[["hierarchies"]]
  df_microdata <- all[["df_microdata"]]
  df_merged <- all[["df_merged"]]
  
  
  if("gauss" %in% df_merged$method) {
    rlang::inform(c(
      v = "gauss output already included",
      i = "Skipping computation"
    ))
    return(invisible(NULL))
  }
  
  
  primary <- df_merged$primary_gauss
  
  
  timing <- system.time({
    out  <- GaussSuppression::SuppressDominantCells(data=df_microdata,
                                  numVar = "response",
                                  hierarchies = hrc_GAUSS,
                                  preAggregate = TRUE,
                                  extraAggregate = FALSE,
                                  primary = primary,
                                  singletonMethod = "none",
                                  removeEmpty = TRUE)
  })
  
  hier_names <- names(hrc_GAUSS)
  ok <- all.equal(out[hier_names], df_merged[hier_names])
  
  if (!isTRUE(ok)) {
    print(ok)
    stop("generated table was not identical")
  }
  
  df_merged <- add_info(df_merged, "gauss", timing)
  
  df_merged$suppressed_gauss <- out$suppressed 
  
  if(identical(output,  "df_merged")){
    return(df_merged)
  }
  
  all[["df_merged"]] <- df_merged
  
  saveRDS(all, file.path(path, paste0(filename, ".rds")))
  
}



