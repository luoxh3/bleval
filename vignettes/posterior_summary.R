# posterior_summary.R
# This script contains utility functions for calculating summary statistics 
# from Bayesian posterior draws. It is used in the example scripts provided 
# in the `example` folder of this package.

# -----------------------------------------------------------------------------
# Function: calculate_hdi
# Description: 
#   Calculates the Highest Density Interval (HDI) for a vector of MCMC samples.
#   The HDI is the smallest interval that contains a specified probability mass.
#
# Arguments:
#   - sampleVec: A numeric vector of MCMC samples.
#   - credible_mass: The probability mass to include in the HDI (default: 0.95).
#
# Returns:
#   A numeric vector of length 2 containing the lower and upper bounds of the HDI.
# -----------------------------------------------------------------------------

calculate_hdi <- function(sampleVec, credible_mass = 0.95) {
  sortedPts <- sort(sampleVec)
  ciIdxInc <- ceiling(credible_mass * length(sortedPts))
  nCIs <- length(sortedPts) - ciIdxInc
  ciWidth <- rep(0, nCIs)
  
  # Find the interval with the smallest width that contains the credible mass
  for (i in 1:nCIs) {
    ciWidth[i] <- sortedPts[i + ciIdxInc] - sortedPts[i]
  }
  HDImin <- sortedPts[which.min(ciWidth)]
  HDImax <- sortedPts[which.min(ciWidth) + ciIdxInc]
  HDIlim <- c(HDImin, HDImax)
  
  return(HDIlim)
}

# -----------------------------------------------------------------------------
# Function: summarize_posterior
# Description:
#   Calculates summary statistics for Bayesian posterior draws, including mean, 
#   median, mode, standard deviation, HDI, credible intervals, ESS, and R-hat.
#
# Arguments:
#   - samples: An MCMC object or a list of MCMC objects.
#   - credible_mass: The probability mass for the HDI (default: 0.95).
#   - credible_interval: The quantiles for the credible interval (default: c(0.025, 0.975)).
#   - delta: A numeric value for delta (default: 0).
#   - filters: A vector of filters to apply to variable names (default: NULL).
#
# Returns:
#   A data frame containing summary statistics for the posterior draws.
# -----------------------------------------------------------------------------

summarize_posterior <- function(samples, credible_mass = 0.95, 
                                credible_interval = c(0.025, 0.975), 
                                delta = 0, filters = NULL) {
  
  sampleCount <- 1
  
  # If only one model is passed, convert it to a list for consistency
  if (!is.mcmc.list(samples)) {
    sampleCount <- length(samples)
  } else {
    samples <- list(samples)
  }
  
  # Ensure filters are provided as a vector
  if (!is.null(filters) && !is.vector(filters)) {
    filters <- c(filters)
  }
  
  # Pre-calculate the number of expected rows in the result data frame
  rowCount <- 0
  drops <- list()
  for (k in 1:sampleCount) {
    drops[[k]] <- c(0)
    if (is.null(filters)) {
      rowCount <- rowCount + nvar(samples[[k]])
    } else {
      nvar <- nvar(samples[[k]])
      varnames <- varnames(samples[[k]])
      for (l in 1:nvar) {
        varname <- varnames[[l]]
        isOK <- FALSE
        for (f in 1:length(filters)) {
          isOK <- isOK || regexpr(filters[[f]], varname)[1] > 0
        }
        if (isOK) {
          rowCount <- rowCount + 1
        } else {
          drops[[k]] <- c(drops[[k]], l)
        }
      }
    }
  }
  
  columnNames <- c()
  
  # Initialize an empty data frame to store the results
  result <- data.frame(
    mean = rep(NaN, rowCount),
    median = rep(NaN, rowCount),
    mode = rep(NaN, rowCount),
    sd = rep(NaN, rowCount),
    hdiLow = rep(NaN, rowCount),
    hdiHigh = rep(NaN, rowCount),
    quantileLow = rep(NaN, rowCount),
    quantileHigh = rep(NaN, rowCount),
    SS = rep(NaN, rowCount),
    ESS = rep(NaN, rowCount),
    RHAT = rep(NaN, rowCount),
    stringsAsFactors = FALSE
  )
  
  # Track the current row in the result data frame
  currentRow <- 0
  
  # Loop over each model (if multiple models are provided)
  for (k in 1:sampleCount) {
    
    # Add a prefix to variable names if multiple models are present
    prefix <- ""
    if (sampleCount > 1) {
      prefix <- paste(k, ".", sep = "")
    }
    
    # Extract the current MCMC sample
    sample <- samples[[k]]
    
    # Get metadata about the sample
    variables <- nvar(sample)
    varnames <- varnames(sample)
    iterations <- niter(sample)
    chains <- nchain(sample)
    
    # Loop through each variable in the sample
    for (j in 1:variables) {
      
      # Skip variables that are filtered out
      if (!(j %in% drops[[k]])) {
        currentRow <- currentRow + 1
        
        # Extract the samples for the current variable
        uvalue <- unlist(sample[, j])
        value <- sample[, j]
        
        # Store the variable name
        columnNames <- c(columnNames, paste(prefix, varnames[[j]], sep = ""))
        
        # Calculate sample size (SS)
        result[currentRow, "SS"] <- iterations * chains
        
        # Calculate Effective Sample Size (ESS)
        tryCatch({
          result[currentRow, "ESS"] <- as.integer(round(effectiveSize(uvalue), 1))
        }, error = function(e) {
          result[currentRow, "ESS"] <- 333333
        })
        
        # Calculate mean, median, and mode
        result[currentRow, "mean"] <- mean(uvalue, na.rm = TRUE)
        result[currentRow, "median"] <- median(uvalue, na.rm = TRUE)
        mcmcDensity <- density(uvalue)
        result[currentRow, "mode"] <- mcmcDensity$x[which.max(mcmcDensity$y)]
        
        # Calculate HDI
        HDI <- calculate_hdi(uvalue, credible_mass)
        result[currentRow, "hdiLow"] <- HDI[1]
        result[currentRow, "hdiHigh"] <- HDI[2]
        
        # Calculate credible interval (quantile)
        resultCI <- quantile(uvalue, credible_interval)
        result[currentRow, "quantileLow"] <- resultCI[1]
        result[currentRow, "quantileHigh"] <- resultCI[2]
        
        # Calculate standard deviation
        result[currentRow, "sd"] <- sd(uvalue, na.rm = TRUE)
        
        # Calculate R-hat (Gelman-Rubin statistic)
        chainmeans <- c()
        chainvars <- c()
        for (i in 1:chains) {
          sum <- sum(value[[i]])
          var <- var(value[[i]])
          mean <- sum / iterations
          chainmeans <- c(chainmeans, mean)
          chainvars <- c(chainvars, var)
        }
        globalmean <- sum(chainmeans) / chains
        globalvar <- sum(chainvars) / chains
        
        # Compute between- and within-variances
        b <- sum((chainmeans - globalmean)^2) * iterations / (chains - 1)
        varplus <- (iterations - 1) * globalvar / iterations + b / iterations
        
        # Handle edge cases for R-hat calculation
        if (globalvar == 0 | is.infinite(globalvar)) {
          rhat <- 333
        } else {
          rhat <- sqrt(varplus / globalvar)
        }
        result[currentRow, "RHAT"] <- rhat
      }
    }
  }
  
  # Round the results to 4 decimal places
  result <- data.frame(apply(result, 2, function(x) round(x, 4)))
  
  # Rename columns for better readability
  if (length(result) > 0) {
    names(result)[names(result) == 'hdiLow'] <- paste(sprintf("%.0f", round(credible_mass * 100, digits = 2)), "HDI_L", sep = "% ")
    names(result)[names(result) == 'hdiHigh'] <- paste(sprintf("%.0f", round(credible_mass * 100, digits = 2)), "HDI_H", sep = "% ")
    names(result)[names(result) == 'quantileLow'] <- paste("CrI", sprintf("%.1f%%", round(credible_interval[1] * 100, digits = 3)), sep = " ")
    names(result)[names(result) == 'quantileHigh'] <- paste("CrI", sprintf("%.1f%%", round(credible_interval[2] * 100, digits = 3)), sep = " ")
  }
  
  # Set row names to variable names
  row.names(result) <- columnNames
  
  # Return the result data frame
  return(result)
}

