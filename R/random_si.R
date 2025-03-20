#' Title Construction of selection index
#'
#' @param data A data frame containing the fixed, random and traits
#' @param traits The traits for which index values are to be estimated
#' @param random The random effects
#' @param economic_values The relative economic values
#'
#' @returns Results of selection index
#' @export
#' @import dplyr
#' @import psych
#' @examples results <- mixed_si(data = data, traits = traits,
#' random = random, economic_values = economic_values)
#'

random_si <- function(data_input, traits, economic_values) {

  # Input data
  if (is.character(data_input)) {
    # If input is a file path
    if (!file.exists(data_input)) {
      stop("Data file does not exist")
    }
    data <- try(read.csv(data_input), silent = TRUE)
    if (inherits(data, "try-error")) {
      stop("Error reading the data file")
    }
  } else if (is.data.frame(data_input)) {
    # If input is already a data frame
    data <- data_input
  } else {
    stop("data_input must be either a file path or a data frame")
  }

  # Input validation
  if (length(traits) != length(economic_values)) {
    stop("Number of traits must match number of economic values")
  }

  # Validate if all traits exist in the data
  if (!all(traits %in% names(data))) {
    missing_traits <- traits[!traits %in% names(data)]
    stop("Following traits not found in data: ", paste(missing_traits, collapse = ", "))
  }

  # Check if sire column exists
  if (!"sire" %in% names(data)) {
    stop("Data must contain a 'sire' column")
  }

  # Initialize results data frame
  num_traits <- length(traits)
  results <- data.frame(
    Trait = traits,
    Heritability = numeric(num_traits),
    AdditiveGeneticVariance = numeric(num_traits),
    ResidualVariance = numeric(num_traits),
    PhenotypicVariance = numeric(num_traits)
  )

  # To initialize correlation and covariance matrices
  genetic_cor_matrix <- matrix(NA, num_traits, num_traits,
                               dimnames = list(traits, traits))
  phenotypic_cor_matrix <- matrix(NA, num_traits, num_traits,
                                  dimnames = list(traits, traits))
  genetic_cov_matrix <- matrix(NA, num_traits, num_traits,
                               dimnames = list(traits, traits))
  phenotypic_cov_matrix <- matrix(NA, num_traits, num_traits,
                                  dimnames = list(traits, traits))

  # To convert sire to factor
  data$sire <- as.factor(data$sire)
  n <- nrow(data)
  p <- length(unique(data$sire))
  q <- n - p

  # To calculate variance components for a single trait
  calculate_variance_components <- function(y, Z) {
    X <- rep(1, n)
    F_val <- crossprod(y, X) * (solve(crossprod(X)) * crossprod(X, y))
    T_val <- (crossprod(y, Z) %*% solve(crossprod(Z)) %*% crossprod(Z, y)) - F_val
    Residual <- crossprod(y) - (T_val + F_val)

    I <- diag(n)
    QR <- I - X %*% solve(crossprod(X)) %*% t(X)
    QRZ <- t(Z) %*% QR %*% Z
    K <- sum(diag(QRZ)) / (p - 1)

    Sigma_R <- Residual / q
    Sigma_S <- (T_val / (p - 1) - Sigma_R) / K

    return(list(Sigma_S = Sigma_S, Sigma_R = Sigma_R))
  }

  # Calculate Z matrix once
  Z <- model.matrix(~ data$sire - 1)

  # To calculate heritability and variance components for each trait
  for (i in seq_along(traits)) {
    y <- data[[traits[i]]]
    components <- calculate_variance_components(y, Z)

    results$AdditiveGeneticVariance[i] <- components$Sigma_S
    results$ResidualVariance[i] <- components$Sigma_R
    results$PhenotypicVariance[i] <- components$Sigma_S + components$Sigma_R
    results$Heritability[i] <- 4 * (components$Sigma_S /
                                      (components$Sigma_S + components$Sigma_R))
  }

  # To calculate correlations and covariances
  for (i in 1:num_traits) {
    for (j in i:num_traits) {
      if (i == j) {
        genetic_cor_matrix[i, i] <- 1
        phenotypic_cor_matrix[i, i] <- 1
        genetic_cov_matrix[i, i] <- results$AdditiveGeneticVariance[i]
        phenotypic_cov_matrix[i, i] <- results$PhenotypicVariance[i]
      } else {
        y_i <- data[[traits[i]]]
        y_j <- data[[traits[j]]]

        components <- calculate_variance_components(y_i + y_j, Z)
        Sigma_S3 <- (components$Sigma_S -
                       results$AdditiveGeneticVariance[i] -
                       results$AdditiveGeneticVariance[j]) / 2
        Sigma_R3 <- (components$Sigma_R -
                       results$ResidualVariance[i] -
                       results$ResidualVariance[j]) / 2

        # Genetic correlation
        GC <- Sigma_S3 / sqrt(results$AdditiveGeneticVariance[i] *
                                results$AdditiveGeneticVariance[j])

        # Phenotypic correlation and covariance
        PCoV_ij <- Sigma_S3 + Sigma_R3
        PC <- PCoV_ij / sqrt(results$PhenotypicVariance[i] *
                               results$PhenotypicVariance[j])

        # Fill matrices
        genetic_cor_matrix[i, j] <- genetic_cor_matrix[j, i] <- GC
        phenotypic_cor_matrix[i, j] <- phenotypic_cor_matrix[j, i] <- PC
        genetic_cov_matrix[i, j] <- genetic_cov_matrix[j, i] <- Sigma_S3
        phenotypic_cov_matrix[i, j] <- phenotypic_cov_matrix[j, i] <- PCoV_ij
      }
    }
  }

  # Calculate selection index
  REV <- matrix(economic_values, ncol = 1)
  GREV <- genetic_cov_matrix %*% REV
  PCM <- solve(phenotypic_cov_matrix)
  SI <- PCM %*% GREV

  # Accuracy of Index
  SI <- as.matrix(c(SI))
  REV <- (as.matrix(economic_values))
  rti1 <- t(SI) %*% phenotypic_cov_matrix %*% SI
  rti2 <- t(REV) %*%  genetic_cov_matrix %*% REV
  RTI <- sqrt(rti1 / rti2)

  # Reliability of Index
  SI <- as.matrix(c(SI))
  Rti1 <- t(SI) %*% genetic_cov_matrix %*% SI
  Rti2 <- t(SI) %*% phenotypic_cov_matrix %*% SI
  RTI3 <- Rti1 / Rti2

  #Expected genetic Change in Traits
  Nume <- crossprod(SI, genetic_cov_matrix)
  rti3 <- sqrt(rti1)
  Deno <- rep(rti3, each = num_traits)
  Genetic_Gain <- Nume / Deno

  # Results
  return(list(
    Results = results,
    GeneticCorrelationMatrix = genetic_cor_matrix,
    PhenotypicCorrelationMatrix = phenotypic_cor_matrix,
    GeneticCovarianceMatrix = genetic_cov_matrix,
    PhenotypicCovarianceMatrix = phenotypic_cov_matrix,
    SelectionIndex = SI,
    AccuracyofIndex = RTI,
    ReliabilityofIndex = RTI3,
    GeneticGain = Genetic_Gain
  ))
}
# examples
data("data", package = "seliNDRIx")

results <- random_si(data,
                     traits = c("tmy", "py", "fatyield"),
                     economic_values = c(1, 0.85, 0.65))
results[] <- lapply(results, function(x) if(is.numeric(x) || is.matrix(x)) round(x, 3) else x)

results
# To calculate the overall selection index for each animal
SI <- c(results$SelectionIndex)  # Selection index estimates (weights) for traits
traits <- c("tmy", "py", "fatyield")  # Define the trait columns to use
overall_index <- function(data, SI, traits) {
  # Ensure the number of weights matches the number of trait columns
  if (length(SI) != length(traits)) {
    stop("The number of weights must match the number of trait columns.")
  }

  # Select only the defined trait columns and calculate the index
  data %>%
    rowwise() %>%
    mutate(Index = sum(c_across(all_of(traits)) * SI)) %>%
    ungroup()
}

# Calculate the selection index
result3 <- overall_index(data, SI, traits)

# Print the result
print(result3)

# Select the top 20% of animals with the highest selection index values
top20 <- result3 %>%
  arrange(desc(Index)) %>%  # Sort by Index in descending order
  slice_head(prop = 0.2)    # Select the top 20%

write.csv(top20, "D:/Top20.csv", row.names = FALSE)
message("Please ensure that the heritability, phenotypic and genetic correlation estimates are well within their normal range. If any estimate exceeds the range, handle the results with caution.")

