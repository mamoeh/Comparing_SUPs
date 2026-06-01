# To add analysis to the .rds file.
#
# Use output = "analysis" to return analysis instead  
# 

add_analysis <- function(filename, 
                      path = "merged", 
                      output = NULL) {
  #Load
  all <- readRDS(file.path(path, paste0(filename, ".rds")))
  #initialise analysis list
  analysis <- list()
  
  
  result <- all$df_merged |> fix_empty()
  
  
  
  #add importance score for each cell
  result <- result |> 
    dplyr::mutate(importance_score = response / 10^(5 * mean_n_at))
  
  #identify secondary suppressed cells for each method:
  #GAUSS
  secondary_gauss <- result |>
    dplyr::filter(suppressed_gauss =="TRUE" & primary_gauss=="FALSE")
  #MODULAR
  secondary_modular <- result |>
    dplyr::filter(suppressed_modular=="TRUE" & primary_modular=="FALSE")
  #Simpleheuristic
  secondary_simple <- result |> 
    dplyr::filter(suppressed_simpleheuristic=="TRUE" & primary_simpleheuristic=="FALSE")
  #Simpleheuristic_old
  secondary_simple_old <- result |> 
    dplyr::filter(suppressed_simpleheuristic_old=="TRUE" & primary_simpleheuristic_old=="FALSE")
  #Number of secondaries per method
  secondaries_gauss <- nrow(secondary_gauss)
  secondaries_modular <- nrow(secondary_modular)
  secondaries_simple <- nrow(secondary_simple)
  secondaries_simple_old <- nrow(secondary_simple_old)
  #add to list
  analysis["no_secondaries"] <- list(c(secondaries_gauss,
                                                   secondaries_modular,
                                                   secondaries_simple,
                                                   secondaries_simple_old))
  
  #(absolute) suppressed response value for each method:
  value_surpressed_gauss <- sum(secondary_gauss$response)
  value_surpressed_modular <- sum(secondary_modular$response)
  value_surpressed_simple <- sum(secondary_simple$response)
  value_surpressed_simple_old <- sum(secondary_simple_old$response)
  #add to list
  analysis["absolute_value"] <- list(c(value_surpressed_gauss,
                                        value_surpressed_modular,
                                        value_surpressed_simple,
                                        value_surpressed_simple_old))
  #(relative to response of all cells) suppressed response value for each method:
  rel_value_surpressed_gauss <- sum(secondary_gauss$response)/sum(result$response)
  rel_value_surpressed_modular <- sum(secondary_modular$response)/sum(result$response)
  rel_value_surpressed_simple <- sum(secondary_simple$response)/sum(result$response)
  rel_value_surpressed_simple_old <- sum(secondary_simple_old$response)/sum(result$response)
  
  #add to list
  analysis["relative_value"] <- list(c(rel_value_surpressed_gauss,
                                       rel_value_surpressed_modular,
                                       rel_value_surpressed_simple,
                                       rel_value_surpressed_simple_old))
  
  #Hierarchical damage
  #suppressed importance relative to response of all cells
  importance_value_suppressed_gauss <- sum(secondary_gauss$importance_score)
  importance_value_suppressed_modular <- sum(secondary_modular$importance_score)
  importance_value_suppressed_simple <- sum(secondary_simple$importance_score)
  importance_value_suppressed_simple_old <- sum(secondary_simple_old$importance_score)
  
  analysis["importance_value_surpressed"]<- list(c(importance_value_suppressed_gauss,
                                                   importance_value_suppressed_modular,
                                                   importance_value_suppressed_simple,
                                                   importance_value_suppressed_simple_old))
  
  #Protectedness:
  #For how many primary suppressed cells [lo, up] is OK, i.e  
  #[lomax_response, upmin_respons] is covered by [lo, up]?
    
  #GAUSS
  Protectedness_gauss <- result |>
    dplyr::filter(lomax_response > lo_gauss & upmin_response < up_gauss)
  #in percent
  Protectedness_gauss_percent <- (nrow(Protectedness_gauss)/ sum(!is.na(result$lo_gauss))) *100
  #MODULAR
  Protectedness_modular <- result |>
    dplyr::filter(lomax_response > lo_modular & upmin_response < up_modular)
  #in percent
  Protectedness_modular_percent <- (nrow(Protectedness_modular)/ sum(!is.na(result$lo_modular))) *100
  #simple
  Protectedness_simple <- result |>
    dplyr::filter(lomax_response > lo_simpleheuristic & upmin_response < up_simpleheuristic)
  #in percent
  Protectedness_simple_percent <- (nrow(Protectedness_simple)/ sum(!is.na(result$lo_simpleheuristic ))) *100
  #simple
  Protectedness_simple <<- Protectedness_simple
  result <<- result
  Protectedness_simple_old <- result |>
    dplyr::filter(lomax_response > lo_simpleheuristic_old & upmin_response < up_simpleheuristic_old)
  #in percent
  Protectedness_simple_old_percent <- (nrow(Protectedness_simple_old)/ sum(!is.na(result$lo_simpleheuristic_old))) *100
  
  #add to list
  analysis["Protectedness"] <- list(c(Protectedness_gauss_percent,
                                       Protectedness_modular_percent,
                                       Protectedness_simple_percent,
                                       Protectedness_simple_old_percent))
  
  
  
  # - Percentage disclosures based on exact linear relations.
  # - These disclosures are not reported for gauss (shown as NA), 
  #     since they are exactly what gauss protects against. 
  #     In fact, the calculations are performed using the gauss algorithm.
  # - The lower and upper interval bounds should be identical if the calculations are correct.
  analysis[["exact_linear_reveal_percent"]] <-
    exact_linear_reveal_percent(all$df_merged, 
                                c("gauss", "modular", "simpleheuristic", "simpleheuristic_old"))
  
  
  #suppressions by "Class" to evaluate hierarchical damage
  #For GAUSS
  Class_gauss <- result |> 
    dplyr::filter(primary_gauss==F & suppressed_gauss == T) |> 
    dplyr::count(HiTaS_Class, name = "n_GAUSS")
  #For Modular
  Class_modular <- result|> 
    dplyr::filter(primary_modular ==F & suppressed_modular ==T) |> 
    dplyr::count(HiTaS_Class, name = "n_modular")
  #For Simple
  Class_simple <- result|> 
    dplyr::filter(primary_simpleheuristic ==F & suppressed_simpleheuristic ==T) |> 
    dplyr::count(HiTaS_Class, name = "n_simple")
  #For Simple Old
  Class_simple_old <- result|> 
    dplyr::filter(primary_simpleheuristic_old ==F & suppressed_simpleheuristic_old ==T) |> 
    dplyr::count(HiTaS_Class, name = "n_simple_old")
  #for comparing the two
  Classes_table <- dplyr::full_join(dplyr::full_join(dplyr::full_join(Class_gauss, Class_modular), Class_simple),Class_simple_old)
  
  #add to list
  analysis["Suppressions_by_Class"] <- list(Classes_table)
  
  if (identical(output, "analysis")) {
    return(analysis)
  }
  
  #add to instance and save
  all[["analysis"]] <- analysis
  saveRDS(all, file.path(path, paste0(filename, ".rds")))
  
}


fix_empty <- function(result) {
  check_names <- c("lo_gauss", "lo_modular", "lo_simpleheuristic", "lo_simpleheuristic_old", 
    "primary_modular", "primary_simpleheuristic", "primary_simpleheuristic_old", 
    "suppressed_gauss", "suppressed_modular", "suppressed_simpleheuristic", 
    "suppressed_simpleheuristic_old", "unsafe_modular", "unsafe_simpleheuristic", 
    "unsafe_simpleheuristic_old", "up_gauss", "up_modular", "up_simpleheuristic", 
    "up_simpleheuristic_old")
  for (nam in check_names) {
    if (is.null(result[[nam]])) {
      result[[nam]] <- 0L
    }
  }
  result
}



exact_linear_reveal_percent <- function(df, name) {
  if (length(name) > 1) {
    return(c(exact_linear_reveal_percent(df, name[1]), exact_linear_reveal_percent(df, name[-1])))
  }
  unsafe <- df[[paste0("unsafe_", name)]]
  primary <- df[[paste0("primary_", name)]]
  if (is.null(unsafe)) {
    res <- NA
  } else {
    res <- 100 * sum(unsafe[primary])/sum(primary)
  }
  names(res) <- name
  res
}


