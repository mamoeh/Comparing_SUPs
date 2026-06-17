library(dplyr)


##### Matrix of how suppressions by different methods overlap ################ #
# df_m:    df_merged from the RDS result files (single instance df or pooled using pool_RDS)
# dig:     number of decimal digits to show for the share of overlap
# methods: which methods to include
############################################################################## #

make_matrix_commons <- function(df_m, dig = 2, 
                                methods_suppr = c("gauss", "modular", "simpleheuristic", "simpleheuristic_old")) {
  
  # check for available methods in df_merged
  methods_available <- paste0("suppressed_", methods_suppr) %in% names(df_m)
  
  # only analyze secondary suppressions
  for(i in seq(methods_suppr)) {
    if(methods_available[i]) {
      
      df_m[, paste0("secondary_", methods_suppr[i])] <-
        df_m[, paste0("suppressed_", methods_suppr[i])] &
        !df_m[, paste0("primary_", methods_suppr[i])]
    }
  }
  
  pooled <- "instance" %in% names(df_m) # is a pooled df_m?
  
  # span the matrix of method x vs. method y (including x == y)
  matr <- expand.grid(suppr_from = methods_suppr[methods_available],
                      also_in    = methods_suppr[methods_available], 
                      count = NA, count_total = NA, count_share = NA,
                      value = NA, value_total = NA, value_share = NA, ntabs = NA)
  
  # fill the matrix
  for(i in 1:nrow(matr)) {
    
    meth1 <- paste0("secondary_", matr$suppr_from[i])
    meth2 <- paste0("secondary_", matr$also_in[i])
    
    # count only in instances where both methods are available
    df_m_sub <- df_m[!is.na(df_m[, meth1]) & !is.na(df_m[, meth2]), ]
    
    matr$count[i] <- sum(df_m_sub[, meth1] & df_m_sub[, meth2], na.rm = TRUE)
    matr$value[i] <- sum((df_m_sub[, meth1] & df_m_sub[, meth2]) * df_m_sub$response, na.rm = TRUE)
    
    matr$count_total[i] <- sum(df_m_sub[, meth1], na.rm = TRUE)
    matr$value_total[i] <- sum(df_m_sub[, meth1] * df_m_sub$response, na.rm = TRUE)
    
    matr$ntabs[i] <- ifelse(pooled, length(unique(df_m_sub$instance)), 1)
  }
  
  # calculate relative values
  matr$count_share <- round(matr$count / matr$count_total, digits = dig)
  matr$value_share <- round(matr$value / matr$value_total, digits = dig)
  
  matr |>
    dplyr::arrange(suppr_from, also_in)
}

