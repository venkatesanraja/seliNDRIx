#' seliNDRIx: A Package for construction of selection index using mixed and
#' random least squares analysis
#'
#' @description
#' This package provides functions for mixed and random least squares analysis
#' to generate the genetic and phenotypic parameters to be used for the
#' construction of selection index
#' Contains two main functions: mixed() and random() that perform
#' the least squares analysis to generate the heritability, genetic and
#' phenotypic correlations of the traits along with their variances and
#' covariances
#'
#' @details
#' The package includes the following main functions:
#' \itemize{
#'   \item \code{mixed()}: Mixed model least squares analysis to estimate the
#'   genetic and phenotypic parameters
#'   \item \code{random()}: Random model least squares analysis to estimate the
#'   genetic and phenotypic parameters
#' }
#'
#' @docType _PACKAGE
#' @name seliNDRIx
#' @aliases seliNDRIx-package
#'#'
#' @examples
#' # Example using mixed function
#' # Read the data
#' data("data", package = "seliNDRIx")
#' # Define your parameters
#' traits <- c("tmy", "py", "fatyield")
#' fixed <- c("farm", "soc", "poc")
#' random <- c("sire")
#' economic_values <- c(1, 0.85, 0.65)

#' # Run the analysis
#' results <- mixed_si(
#' data = data,
#' traits = traits,
#' fixed = fixed,
#' random = random,
#' economic_values = economic_values
#' )
#' results
#' # To calculate the overall selection index for each animal
#' SI <- c(results$SelectionIndex)  # Selection index estimates (weights) for traits
#' traits <- c("tmy", "py", "fatyield")  # Define the trait columns to use
#' overall_index <- function(data, SI, traits) {
#' # Ensure the number of weights matches the number of trait columns
#'  if (length(SI) != length(traits)) {
#'    stop("The number of weights must match the number of trait columns.")
#'  }

#'  # Select only the defined trait columns and calculate the index
#'  data %>%
#'    rowwise() %>%
#'    mutate(Index = sum(c_across(all_of(traits)) * SI)) %>%
#'    ungroup()
#' }

#' # Calculate the selection index
#' result3 <- overall_index(data, SI, traits)

#' # Print the result
#' # print(result3)

#' # Select the top 20% of animals with the highest selection index values
#' top20 <- result3 %>%
#'  arrange(desc(Index)) %>%  # Sort by Index in descending order
#'  slice_head(prop = 0.2)    # Select the top 20%

#' # Example using random function
#' # Read the data
#' data("data", package = "seliNDRIx")
#' # Run the analysis
#' results2 <- random_si(data,
#' traits = c("tmy", "py", "fatyield"),
#' economic_values = c(1, 0.85, 0.65))
#' results2
#' # To calculate the overall selection index for each animal
#' SI <- c(results$SelectionIndex)  # Selection index estimates (weights) for traits
#' traits <- c("tmy", "py", "fatyield")  # Define the trait columns to use
#' overall_index <- function(data, SI, traits) {
#' # Ensure the number of weights matches the number of trait columns
#'  if (length(SI) != length(traits)) {
#'    stop("The number of weights must match the number of trait columns.")
#'  }

#'  # Select only the defined trait columns and calculate the index
#'  data %>%
#'    rowwise() %>%
#'    mutate(Index = sum(c_across(all_of(traits)) * SI)) %>%
#'    ungroup()
#' }

#' # Calculate the selection index
#' result3 <- overall_index(data, SI, traits)

#' # Print the result
#' # print(result3)

#' # Select the top 20% of animals with the highest selection index values
#' top20 <- result3 %>%
#'  arrange(desc(Index)) %>%  # Sort by Index in descending order
#'  slice_head(prop = 0.2)    # Select the top 20%

NULL
