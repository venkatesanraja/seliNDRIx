#' Title Construction of selection index
#'
#' @param data A data frame containing the fixed, random and traits
#' @param traits The traits for which index values are to be estimated
#' @param fixed The fixed effects
#' @param random The random effects
#' @param economic_values The relative economic values
#'
#' @returns Results of selection index
#' @export
#' @import dplyr
#' @import psych
#' @examples results <- mixed_si(data = data, traits = traits,
#' fixed = fixed, random = random, economic_values = economic_values)
#'

mixed_si <- function(data, traits, fixed, random, economic_values) {
  # Input validation
  if (!is.data.frame(data)) {
    stop("Data must be a data frame")
  }
  if (!all(traits %in% names(data))) {
    stop("All traits must be present in the data frame")
  }
  if (!all(fixed %in% names(data))) {
    stop("All fixed effects must be present in the data frame")
  }
  if (!all(random %in% names(data))) {
    stop("All random effects must be present in the data frame")
  }
  if (length(traits) != length(economic_values)) {
    stop("Number of economic values must match number of traits")
  }

  require(dplyr)
  require(psych)
  # Initialize variables
  num_traits <- length(traits)
  RHS_list <- list()

  # Initialize results dataframe
  results <- data.frame(
    Trait = traits,
    Heritability = NA,
    AdditiveGeneticVariance = NA,
    ResidualVariance = NA,
    PhenotypicVariance = NA,
    row.names = traits
  )

  # Convert effects to factors
  for (effect in c(fixed, random)) {
    data[[effect]] <- as.factor(data[[effect]])
  }

  # Initialize all_results list
  all_results <- list()

  # Reduce matrix function
  reduce_matrix <- function(matrix_data, effects) {
    for (effect in effects) {
      effect_rows <- grep(paste0("^", effect), rownames(matrix_data))
      effect_cols <- grep(paste0("^", effect), colnames(matrix_data))

      if (length(effect_rows) > 1) {
        last_row <- effect_rows[length(effect_rows)]
        for (k in effect_rows[-length(effect_rows)]) {
          matrix_data[k, ] <- matrix_data[k, ] - matrix_data[last_row, ]
        }
        matrix_data <- matrix_data[-last_row, ]
      }

      if (length(effect_cols) > 1) {
        last_col <- effect_cols[length(effect_cols)]
        for (j in effect_cols[-length(effect_cols)]) {
          matrix_data[, j] <- matrix_data[, j] - matrix_data[, last_col]
        }
        matrix_data <- matrix_data[, -last_col]
      }
    }
    return(matrix_data)
  }

  # Variance components and heritability calculation
  for (i in seq_along(traits)) {
    trait_name <- traits[i]

    y <- data[[trait_name]]
    n <- nrow(data)

    # Generate design matrices for fixed and random effects
    X1 <- rep(1, nrow(data))  # Intercept
    fixed_matrices <- lapply(fixed, function(effect) model.matrix(~ data[[effect]] - 1))
    random_matrices <- lapply(random, function(effect) model.matrix(~ data[[effect]] - 1))
    X <- do.call(cbind, c(list(X1), fixed_matrices, random_matrices))

    # Create and reduce R matrix
    R <- crossprod(X)
    rownames(R) <- colnames(R) <- c(
      "Intercept",
      unlist(lapply(fixed, function(effect) paste(effect, levels(data[[effect]]), sep = "_"))),
      unlist(lapply(random, function(effect) paste(effect, levels(data[[effect]]), sep = "_")))
    )

    # Reduce matrix
    P <- reduce_matrix(R, c(fixed, random))

    # Generate RHS components for each trait
    RHS <- list(overall = sum(y))
    for (effect in c(fixed, random)) {
      effect_sums <- tapply(y, data[[effect]], sum)
      effect_sums_numeric <- as.numeric(effect_sums)
      if (length(effect_sums_numeric) > 1) {
        effect_sums_numeric <- effect_sums_numeric - effect_sums_numeric[length(effect_sums_numeric)]
        effect_sums_numeric <- effect_sums_numeric[-length(effect_sums_numeric)]
      }
      RHS[[effect]] <- effect_sums_numeric
    }
    RHS_Final <- matrix(c(RHS$overall, unlist(RHS[c(fixed, random)])), ncol = 1)

    # Generate incidence matrix for fixed effects
    X_fixed <- cbind(1, do.call(cbind, lapply(fixed, function(effect) {
      X_effect <- model.matrix(~ data[[effect]] - 1)
      X_effect[, -ncol(X_effect), drop = FALSE]  # Remove last column
    })))

    # Generate incidence matrix for random effect
    random_effect <- random[length(random)]
    n_levels <- length(unique(data[[random_effect]]))
    Z <- model.matrix(~ data[[random_effect]] - 1)

    # Store all components in all_results list
    all_results[[trait_name]] <- list(
      RHS_Final = RHS_Final,
      R = P,
      y = as.numeric(y),
      X_fixed = X_fixed,
      Z = Z,
      L = solve(P)
    )

    # Solve for J
    L <- solve(P)
    J <- RHS_Final

    # Calculate matrix for all fixed effects
    I <- diag(nrow(data))
    X_fixed_proj <- X_fixed %*% solve(crossprod(X_fixed)) %*% t(X_fixed)
    QR <- I - X_fixed_proj

    # Calculate K value
    QRZ <- t(Z) %*% QR %*% Z
    K <- sum(diag(QRZ)) / (n_levels - 1)

    Solution <- L %*% J
    df_fixed <- sum(sapply(fixed, function(e) length(levels(data[[e]])) - 1))
    df_random <- n_levels - 1
    df_residual <- n - df_random - df_fixed - 1

    SSAnova <- t(J) %*% Solution
    SSTotal <- crossprod(y)
    SSResidual <- SSTotal - SSAnova
    SS_Random <- sum((crossprod(y, X_fixed) %*% solve(crossprod(X_fixed)) %*% crossprod(X_fixed, y)))
    SS_Fixed <- sum((crossprod(y, Z) %*% solve(crossprod(Z)) %*% crossprod(Z, y)))
    SSRandom <- SSAnova - SS_Random

    # Variance calculation
    Sigma_R <- max(0, SSResidual / df_residual)
    temp_sigma_s <- ((SSRandom/df_random) - Sigma_R) / K
    Sigma_S <- max(0, temp_sigma_s)
    Sigma_P <- Sigma_R + Sigma_S

    # Heritability estimation
    h2Trait <- ifelse(Sigma_P > 0, 4 * (Sigma_S / Sigma_P), NA)

    # Store results
    results[trait_name, "Heritability"] <- h2Trait
    results[trait_name, "AdditiveGeneticVariance"] <- Sigma_S
    results[trait_name, "ResidualVariance"] <- Sigma_R
    results[trait_name, "PhenotypicVariance"] <- Sigma_P
  }

  # Initialize correlation matrices
  genetic_cor_matrix <- matrix(NA, num_traits, num_traits)
  phenotypic_cor_matrix <- matrix(NA, num_traits, num_traits)
  genetic_cov_matrix <- matrix(NA, num_traits, num_traits)
  phenotypic_cov_matrix <- matrix(NA, num_traits, num_traits)
  rownames(genetic_cor_matrix) <- traits
  colnames(genetic_cor_matrix) <- traits
  rownames(phenotypic_cor_matrix) <- traits
  colnames(phenotypic_cor_matrix) <- traits
  rownames(phenotypic_cov_matrix) <- traits
  colnames(phenotypic_cov_matrix) <- traits
  rownames(genetic_cov_matrix) <- traits
  colnames(genetic_cov_matrix) <- traits

  # Correlation calculations
  for (i in 1:num_traits) {
    for (j in i:num_traits) {
      if (i != j) {
        # Get stored values for traits i and j
        trait_i_data <- all_results[[traits[i]]]
        trait_j_data <- all_results[[traits[j]]]

        y_i <- trait_i_data$y
        y_j <- trait_j_data$y
        RHS_i <- trait_i_data$RHS_Final
        RHS_j <- trait_j_data$RHS_Final
        L <- solve(trait_i_data$R)

        n <- nrow(data)

        # Use stored matrices
        Z_i <- trait_i_data$Z
        Z_j <- trait_j_data$Z
        X_i <- trait_i_data$X_fixed
        X_j <- trait_j_data$X_fixed

        SSAnova3 <- t(RHS_i) %*% L %*% RHS_j
        SSR3 <- crossprod(y_i, y_j)
        SSResidual3 <- SSR3 - SSAnova3
        F3 <- crossprod(y_i, X_i) %*% solve(crossprod(X_i, X_j)) %*% crossprod(X_j, y_j)
        SSRandom3 <- SSAnova3 - F3

        n_levels <- length(unique(data[[random[1]]])) # Using first random effect

        # Calculate K value
        I <- diag(nrow(data))
        X_fixed_proj3 <- X_i %*% solve(crossprod(X_i, X_j)) %*% t(X_j)
        QR3 <- I - X_fixed_proj3

        # Calculate K value
        QRZ3 <- t(Z_i) %*% QR3 %*% Z_j
        K3 <- sum(diag(QRZ3)) / (n_levels - 1)

        Sigma_R3 <- SSResidual3 / df_residual
        Sigma_S3 <- ((SSRandom3 / df_random) - Sigma_R3) / K3

        Sigma_S_i <- results$AdditiveGeneticVariance[i]
        Sigma_S_j <- results$AdditiveGeneticVariance[j]
        Sigma_R_i <- results$ResidualVariance[i]
        Sigma_R_j <- results$ResidualVariance[j]

        GC <- Sigma_S3 / sqrt(Sigma_S_i * Sigma_S_j)
        PV_i <- Sigma_S_i + Sigma_R_i
        PV_j <- Sigma_S_j + Sigma_R_j
        PCoV_ij <- Sigma_S3 + Sigma_R3
        PC <- PCoV_ij / sqrt(PV_i * PV_j)

        genetic_cor_matrix[i, j] <- GC
        phenotypic_cor_matrix[i, j] <- PC
        genetic_cor_matrix[j, i] <- GC
        phenotypic_cor_matrix[j, i] <- PC
        phenotypic_cov_matrix[i, j] <- PCoV_ij
        phenotypic_cov_matrix[j, i] <- PCoV_ij
        genetic_cov_matrix[i, j] <- Sigma_S3
        genetic_cov_matrix[j, i] <- Sigma_S3
      }
    }
  }

  # Update covariance matrices diagonals
  diag(genetic_cor_matrix) <- 1
  diag(phenotypic_cor_matrix) <- 1
  diag(genetic_cov_matrix) <- results$AdditiveGeneticVariance
  diag(phenotypic_cov_matrix) <- results$PhenotypicVariance

  # Selection index calculation
  REV <- t(as.matrix(economic_values))
  GREV <- crossprod(genetic_cov_matrix, t(REV))
  PCM <- solve(phenotypic_cov_matrix)
  SI <- crossprod(PCM, GREV)

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
# Read the data
data("data", package = "seliNDRIx")

# Define your parameters
traits <- c("tmy", "py", "fatyield")
fixed <- c("farm", "soc", "poc")
random <- c("sire")
economic_values <- c(1, 0.85, 0.65)

# Run the analysis
results <- mixed_si(
  data = data,
  traits = traits,
  fixed = fixed,
  random = random,
  economic_values = economic_values
)
#results[] <- lapply(results, function(x) if(is.numeric(x) || is.matrix(x)) round(x, 3) else x)
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
