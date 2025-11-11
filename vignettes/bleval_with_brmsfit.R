####################################################################################
# General wrapper for bleval + brms (Gaussian/Poisson/Bernoulli, two-level, any K) #
####################################################################################

# =============================================================================
# Section 1. Metadata extraction and design utilities
# =============================================================================

#' Extract metadata from a brms fit for a two-level model
#'
#' Purpose:
#' - Identify the family/link, the primary grouping factor, random-effect terms,
#'   and collect prior summaries and parameter column names used later.
#'
#' Inputs:
#' - fit: brmsfit object; must contain exactly one grouping factor used for
#'   random effects (two-level model).
#'
#' Outputs:
#' - A list with fields:
#'   - family, link: strings
#'   - group: name of the grouping factor selected (the one with most sd_ cols)
#'   - ranef_terms: character vector of random-effect term names (length K)
#'   - K: integer, number of random-effect terms
#'   - sd_cols: vector of sd_ column names for that group
#'   - cor_cols_mat: KxK upper-triangular matrix of cor_ column names (NA when absent)
#'   - prior_summary: output of brms::prior_summary(fit)
#'
#' Processing:
#' - Parses sd_ and cor_ column names from as_draws_df(fit); picks the dominant
#'   group; builds a mapping from term names to parameter columns.
prepare_brms_meta <- function(fit) {
  fam <- family(fit)
  fam_name <- fam$family
  link_name <- fam$link
  ps <- brms::prior_summary(fit)
  draws_df <- posterior::as_draws_df(fit)
  all_sd_cols  <- grep("^sd_", names(draws_df), value = TRUE)
  all_cor_cols <- grep("^cor_", names(draws_df), value = TRUE)
  if (!length(all_sd_cols)) stop("no sd_ cols detected")
  parse_sd <- do.call(rbind, lapply(all_sd_cols, function(s) {
    m <- regexec("^sd_(.+?)__(.+)$", s)
    p <- regmatches(s, m)[[1]]
    if (length(p) == 3) data.frame(col = s, group = p[2], term = p[3], stringsAsFactors = FALSE) else NULL
  }))
  group_counts <- sort(table(parse_sd$group), decreasing = TRUE)
  target_group <- names(group_counts)[1]
  sd_df <- subset(parse_sd, group == target_group)
  ranef_terms <- sd_df$term
  sd_cols <- sd_df$col
  K <- length(ranef_terms)
  cor_cols <- matrix(NA_character_, nrow = K, ncol = K)
  rownames(cor_cols) <- colnames(cor_cols) <- ranef_terms
  for (i in seq_len(K - 1)) for (j in (i + 1):K) {
    term_i_escaped <- gsub("([.|()\\[\\]^$+*?{}\\\\])", "\\\\\\1", ranef_terms[i], perl = TRUE)
    term_j_escaped <- gsub("([.|()\\[\\]^$+*?{}\\\\])", "\\\\\\1", ranef_terms[j], perl = TRUE)
    pat <- paste0("^cor_", target_group, "__", term_i_escaped, "__", term_j_escaped, "$")
    col_ij <- grep(pat, all_cor_cols, value = TRUE, perl = TRUE)
    if (length(col_ij) == 1) cor_cols[i, j] <- col_ij else cor_cols[i, j] <- NA_character_
  }
  list(
    family = fam_name,
    link = link_name,
    group = target_group,
    ranef_terms = ranef_terms,
    K = K,
    sd_cols = sd_cols,
    cor_cols_mat = cor_cols,
    prior_summary = ps
  )
}

#' Construct the random-effects covariance matrix Sigma_u from a sample
#'
#' Purpose:
#' - Given one posterior draw (samples_s), build Sigma_u = D R D where D = diag(sd)
#'   and R is the correlation matrix based on cor_ parameters; supports any K.
#'
#' Inputs:
#' - samples_s: named list-like (e.g., row from samples matrix) with sd_*/cor_* values.
#' - meta: list from prepare_brms_meta() providing ranef_terms and group name.
#'
#' Outputs:
#' - KxK numeric covariance matrix Sigma_u.
#'
#' Processing:
#' - Collect sd values for each term; build diagonal D; fill correlation matrix R
#'   from cor_* entries when present; default missing correlations to 0.
build_Sigma_u_from_sample <- function(samples_s, meta) {
  sd_vals <- numeric(meta$K)
  for (k in seq_len(meta$K)) {
    nm <- paste0("sd_", meta$group, "__", meta$ranef_terms[k])
    sd_vals[k] <- as.numeric(samples_s[[nm]])
  }
  D <- diag(sd_vals, meta$K, meta$K)
  R <- diag(1, meta$K)
  for (i in seq_len(meta$K - 1)) for (j in (i + 1):meta$K) {
    nm <- paste0("cor_", meta$group, "__", meta$ranef_terms[i], "__", meta$ranef_terms[j])
    if (!is.null(samples_s[[nm]])) rij <- as.numeric(samples_s[[nm]]) else rij <- 0
    R[i, j] <- R[j, i] <- rij
  }
  D %*% R %*% D
}

#' Build the random-effects design matrix Z_i for unit i
#'
#' Purpose:
#' - Create Z_i with columns matching meta$ranef_terms, using fit$data to map
#'   covariate columns; the Intercept term is filled with 1s.
#'
#' Inputs:
#' - fit: brmsfit used only to access original data columns
#' - meta: list from prepare_brms_meta(); provides ranef_terms (column order)
#' - data: list with fields subject (indices), and covariates (e.g., x1, x2)
#' - i: integer unit index in 1..N
#'
#' Outputs:
#' - Numeric matrix Z_i of size T_i x K
#'
#' Processing:
#' - For each term, if term == "Intercept" set column to 1; else copy the
#'   corresponding column from fit$data restricted to observations of unit i.
get_design_for_unit <- function(fit, meta, data, i) {
  idx_i <- which(data$subject == i)
  T_i <- length(idx_i)
  Zi <- matrix(NA_real_, nrow = T_i, ncol = meta$K)
  colnames(Zi) <- meta$ranef_terms
  dat <- fit$data
  for (k in seq_len(meta$K)) {
    term <- meta$ranef_terms[k]
    if (term == "Intercept") Zi[, k] <- 1 else Zi[, k] <- dat[[term]][idx_i]
  }
  Zi
}

#' Check if a matrix is symmetric positive-definite via Cholesky
#'
#' Inputs: m numeric matrix
#' Output: logical TRUE/FALSE
#' Processing: reject non-square, NA-containing; try chol(), catch errors.
is_pos_def <- function(m) {
  if (!is.matrix(m)) return(FALSE)
  if (any(is.na(m))) return(FALSE)
  if (nrow(m) != ncol(m)) return(FALSE)
  ok <- tryCatch({chol(m); TRUE}, error = function(e) FALSE)
  isTRUE(ok)
}

# =============================================================================
# Section 2. Family-specific conditional likelihood functions
# =============================================================================

#' Row-wise conditional Gaussian log-likelihood accumulator
#'
#' Inputs:
#' - y_vec: numeric vector (length T_i)
#' - mu_mat: Q x T_i matrix, each row is the mean vector for one quadrature node
#' - sigma: residual SD (scalar > 0)
#'
#' Output: numeric vector length Q with log-likelihoods per node.
loglik_cond_gaussian <- function(y_vec, mu_mat, sigma) {
  rowSums(dnorm(rep(y_vec, each = nrow(mu_mat)), mean = as.vector(mu_mat), sd = sigma, log = TRUE)
          |> matrix(nrow = nrow(mu_mat), byrow = FALSE))
}

#' Row-wise conditional Bernoulli log-likelihood accumulator
#'
#' Inputs:
#' - y_vec: numeric vector (length T_i)
#' - mu_mat: Q x T_i matrix, each row is the mean vector for one quadrature node
#'
#' Output: numeric vector length Q with log-likelihoods per node.
loglik_cond_bernoulli <- function(y_vec, mu_mat) {
  y_mat <- matrix(y_vec, nrow = nrow(mu_mat), ncol = ncol(mu_mat), byrow = TRUE)
  # Calculate log(1 + exp(mu)) in a stable way.
  log1pexp_mu <- ifelse(mu_mat > 0,
                        mu_mat + log1p(exp(-mu_mat)),
                        log1p(exp(mu_mat)))
  
  ll_mat <- y_mat * mu_mat - log1pexp_mu
  rowSums(ll_mat)
}

#' Row-wise conditional Poisson log-likelihood accumulator
#'
#' Inputs:
#' - y_vec: numeric vector (length T_i)
#' - mu_mat: Q x T_i matrix, each row is the mean vector for one quadrature node
#'
#' Output: numeric vector length Q with log-likelihoods per node.
loglik_cond_poisson <- function(y_vec, mu_mat) {
  lambda_mat <- exp(mu_mat)  # mu_mat is linear predictor (on log scale), converted to rate
  rowSums(dpois(rep(y_vec, each = nrow(mu_mat)), lambda = as.vector(lambda_mat), log = TRUE)
          |> matrix(nrow = nrow(mu_mat), byrow = FALSE))
}

# =============================================================================
# Section 3. Prior parsing and evaluation utilities
# =============================================================================

#' Parse a brms prior string like "normal(0,5)" or "lkj(2)"
#'
#' Inputs: prior_str (character)
#' Output: list(dist=<name>, pars=numeric vector) or dist=NA when flat/empty.
#' Processing: normalizes LKJ naming, regex-splits name and parameter list.
parse_prior <- function(prior_str) {
  s <- as.character(prior_str)
  s <- trimws(s)
  if (is.na(s) || s == "" || s == "(flat)") return(list(dist = NA_character_, pars = numeric(0)))
  s2 <- s
  s2 <- sub("^lkj_corr_cholesky", "lkj", s2, ignore.case = TRUE)
  s2 <- sub("^lkj_corr", "lkj", s2, ignore.case = TRUE)
  m <- regexec("^([a-zA-Z_]+)\\s*\\((.*)\\)$", s2)
  r <- regmatches(s2, m)[[1]]
  if (length(r) != 3) return(list(dist = NA_character_, pars = numeric(0)))
  dist <- tolower(r[2])
  pars_raw <- gsub("\\s", "", r[3])
  pars <- if (nchar(pars_raw) == 0) numeric(0) else as.numeric(strsplit(pars_raw, ",")[[1]])
  list(dist = dist, pars = pars)
}

#' Convert optional bound to numeric with defaults
#'
#' Inputs: x (possibly NA/empty/character), default numeric
#' Output: numeric bound (default when invalid)
safe_bound <- function(x, default) {
  if (is.null(x)) return(default)
  xch <- as.character(x)
  if (length(xch) == 0) return(default)
  if (is.na(xch) || !nzchar(trimws(xch))) return(default)
  val <- suppressWarnings(as.numeric(xch))
  if (is.na(val)) return(default)
  val
}

#' Find a matching prior row in prior_summary
#'
#' Purpose: select the appropriate non-flat prior row for a class (b/sd/sigma/cor)
#'          optionally matching by coef and/or group, preferring vectorized entries.
#'
#' Inputs:
#' - ps: data.frame from brms::prior_summary
#' - class_name: e.g., "b", "sd", "sigma", "cor" or "L"
#' - coef_name: optional coefficient name
#' - group_name: optional grouping factor name
#'
#' Output: single-row data.frame or NULL if none found.
find_prior_row <- function(ps, class_name, coef_name = NULL, group_name = NULL) {
  candidates <- ps[as.character(ps$class) == class_name, , drop = FALSE]
  if (nrow(candidates) == 0) return(NULL)
  ac <- function(x) if (is.null(x)) NA_character_ else as.character(x)
  rows <- candidates
  coef_txt <- trimws(ac(rows$coef))
  group_txt <- trimws(ac(rows$group))
  tag_txt <- tolower(trimws(ac(rows$tag)))
  src_txt <- trimws(ac(rows$source))
  prior_txt <- trimws(ac(rows$prior))
  is_nonflat <- !(is.na(prior_txt) | prior_txt == "" | tolower(prior_txt) == "(flat)")
  if (!any(is_nonflat)) return(NULL)
  rows_nf <- rows[is_nonflat, , drop = FALSE]
  coef_nf <- coef_txt[is_nonflat]
  group_nf <- group_txt[is_nonflat]
  tag_nf <- tag_txt[is_nonflat]
  src_nf <- src_txt[is_nonflat]
  if (!is.null(coef_name) && !is.null(group_name)) {
    idx <- which(coef_nf == coef_name & group_nf == group_name)
    if (length(idx) >= 1) return(rows_nf[idx[1], , drop = FALSE])
  }
  if (!is.null(coef_name)) {
    idx <- which(coef_nf == coef_name)
    if (length(idx) >= 1) return(rows_nf[idx[1], , drop = FALSE])
  }
  if (!is.null(group_name)) {
    idx <- which(grepl("vectorized", tag_nf, fixed = TRUE) & group_nf == group_name)
    if (length(idx) >= 1) return(rows_nf[idx[1], , drop = FALSE])
  }
  idx <- which(grepl("vectorized", tag_nf, fixed = TRUE))
  if (length(idx) >= 1) return(rows_nf[idx[1], , drop = FALSE])
  idx <- which(src_nf == "user")
  if (length(idx) >= 1) return(rows_nf[idx[1], , drop = FALSE])
  return(rows_nf[1, , drop = FALSE])
}

#' Log-density of truncated distributions used for priors
#'
#' Supported dist: normal, cauchy, student_t. Applies [lb, ub] truncation.
#'
#' Inputs: x scalar, dist string, pars numeric vector, lb/ub bounds
#' Output: scalar log-density or -Inf on invalid inputs.
log_dtrunc <- function(x, dist, pars, lb = -Inf, ub = Inf) {
  if (any(is.na(x)) || any(is.na(pars))) return(-Inf)
  safe_p <- function(expr) tryCatch(expr, error = function(e) NA_real_)
  if (dist == "normal") {
    mu <- pars[1]; sd <- pars[2]
    if (!is.finite(mu) || !is.finite(sd) || sd <= 0) return(-Inf)
    logf <- dnorm(x, mu, sd, log = TRUE)
    Z <- safe_p(pnorm(ub, mu, sd) - pnorm(lb, mu, sd))
    if (!is.finite(Z) || Z <= 0) return(-Inf)
    return(logf - log(Z))
  }
  if (dist == "cauchy") {
    loc <- pars[1]; sc <- pars[2]
    if (!is.finite(loc) || !is.finite(sc) || sc <= 0) return(-Inf)
    logf <- dcauchy(x, loc, sc, log = TRUE)
    Z <- safe_p(pcauchy(ub, loc, sc) - pcauchy(lb, loc, sc))
    if (!is.finite(Z) || Z <= 0) return(-Inf)
    return(logf - log(Z))
  }
  if (dist == "student_t") {
    df <- pars[1]; mu <- pars[2]; sc <- pars[3]
    if (!is.finite(df) || !is.finite(mu) || !is.finite(sc) || sc <= 0) return(-Inf)
    z <- (x - mu) / sc
    logf <- dt(z, df = df, log = TRUE) - log(sc)
    Z <- safe_p(pt((ub - mu) / sc, df) - pt((lb - mu) / sc, df))
    if (!is.finite(Z) || Z <= 0) return(-Inf)
    return(logf - log(Z))
  }
  return(-Inf)
}

#' LKJ correlation prior log-density for correlation matrix R
#'
#' Inputs: R KxK correlation matrix, eta > 0
#' Output: scalar log-density; returns -Inf if R invalid.
#' Processing: computes normalizing constant and (eta-1)*log(det(R)).
lkj_log_density <- function(R, eta) {
  if (!is.matrix(R) || nrow(R) != ncol(R) || any(is.na(R))) return(-Inf)
  K <- nrow(R)
  detR <- tryCatch(det(R), error = function(e) NA_real_)
  if (!is.finite(detR) || detR <= 0) return(-Inf)
  logC <- (K * (K - 1) / 2) * log(2)
  for (i in 1:(K - 1)) logC <- logC + lbeta(eta + (K - 1 - i) / 2, i / 2)
  logC + (eta - 1) * log(detR)
}

#' Build lv_mu and lv_cov from a brms fit for each group level
#'
#' This function reconstructs per-group posterior mean (lv_mu) and covariance
#' (lv_cov) of the random-effect coefficients theta = b + u from posterior draws
#' of coef(fit). It supports an arbitrary number K of random effects.
#'
#' @param fit brmsfit; fitted model with a single grouping factor.
#' @param meta list; metadata returned by prepare_brms_meta().
#' @return list with elements: lv_mu (list of length N), lv_cov (list of length N).
build_lv_from_fit <- function(fit, meta) {
  co_draws <- coef(fit, summary = FALSE)[[meta$group]]
  term_names <- dimnames(co_draws)[[3]]
  reorder_idx <- match(meta$ranef_terms, term_names)
  N <- dim(co_draws)[2]
  lv_mu <- vector("list", N)
  lv_cov <- vector("list", N)
  for (i in seq_len(N)) {
    Ti <- co_draws[, i, reorder_idx, drop = FALSE]
    Ti_mat <- matrix(Ti, ncol = meta$K)
    lv_mu[[i]] <- colMeans(Ti_mat)
    lv_cov[[i]] <- stats::cov(Ti_mat)
  }
  list(lv_mu = lv_mu, lv_cov = lv_cov)
}

#' Build samples matrix for bleval from a brms fit
#'
#' Keeps only relevant columns: fixed effects b_\*, group-specific sd_\*, cor_\*
#' (and sigma for Gaussian family). Works for any K implied by meta.
#'
#' @param fit brmsfit; fitted model.
#' @param meta list; metadata returned by prepare_brms_meta().
#' @return numeric matrix of posterior draws for required parameters.
build_samples_matrix <- function(fit, meta) {
  draws_df <- posterior::as_draws_df(fit)
  b_cols <- grep("^b_", names(draws_df), value = TRUE)
  sigma_cols <- if (meta$family == "gaussian") grep("^sigma$", names(draws_df), value = TRUE) else character(0)
  sd_cols_group <- meta$sd_cols
  cor_cols_group <- as.vector(meta$cor_cols_mat[upper.tri(meta$cor_cols_mat)])
  cor_cols_group <- cor_cols_group[!is.na(cor_cols_group)]
  keep_cols <- c(b_cols, sd_cols_group, cor_cols_group, sigma_cols)
  keep_cols <- unique(keep_cols[keep_cols %in% names(draws_df)])
  as.matrix(draws_df[, keep_cols, drop = FALSE])
}

#' Build default bounds for bridge sampling
#'
#' Sets bounds for sd_\* and sigma to [0, +Inf), for cor_\* to [-1, 1], and
#' leaves others unconstrained.
#'
#' @param samps_mat numeric matrix; posterior draws as returned by build_samples_matrix().
#' @param meta list; metadata returned by prepare_brms_meta().
#' @return list with lb and ub named numeric vectors.
default_bounds <- function(samps_mat, meta) {
  lb <- rep(-Inf, ncol(samps_mat))
  ub <- rep( Inf, ncol(samps_mat))
  names(lb) <- names(ub) <- colnames(samps_mat)
  sd_idx <- grepl(paste0("^sd_", meta$group, "__"), colnames(samps_mat))
  lb[sd_idx] <- 0
  if (meta$family == "gaussian") {
    sig_idx <- colnames(samps_mat) == "sigma"
    lb[sig_idx] <- 0
  }
  cor_idx <- grepl(paste0("^cor_", meta$group, "__"), colnames(samps_mat))
  lb[cor_idx] <- -1
  ub[cor_idx] <-  1
  list(lb = lb, ub = ub)
}

# =============================================================================
# Section 4. Joint density and prior in brms parameterization
# =============================================================================

#' Per-unit joint log-density for Gaussian models with random effects
#'
#' Purpose:
#' - Evaluate log p(y_i | theta_i, sigma) + log p(theta_i | b, Sigma_u)
#'   at a grid of nodes (rows of 'nodes'), where theta_i corresponds to the
#'   random-effect coefficients (b + u) in brms parameterization.
#'
#' Inputs:
#' - samples_s: one posterior draw (named vector/list with b_, sd_, cor_, sigma)
#' - data: list with y, subject, covariates; i: group index
#' - Ngrid, nodes: quadrature settings (nodes is Q x K)
#' - fit, meta: brms fit and metadata
#' - ...: unused
#'
#' Output: numeric vector length Q with joint log-density values per node.
log_joint_i <- function(samples_s, data, i, Ngrid, nodes, fit, meta, ...) {
  idx_i <- which(data$subject == i)
  y_i <- data$y[idx_i]
  Zi <- get_design_for_unit(fit, meta, data, i)

  b_center <- numeric(meta$K)
  for (k in seq_len(meta$K)) {
    term <- meta$ranef_terms[k]
    bname <- if (term == "Intercept") "b_Intercept" else paste0("b_", term)
    b_center[k] <- as.numeric(samples_s[[bname]])
  }
  
  Sigma_u <- build_Sigma_u_from_sample(samples_s, meta)
  if (!is_pos_def(Sigma_u)) return(rep(-Inf, nrow(nodes)))
  if (any(is.na(nodes)) || !is.numeric(nodes)) return(rep(-Inf, nrow(nodes)))
  
  # Add fixed-effects contribution for observations of unit i
  fixed_part <- rep(0, nrow(Zi))
  b_names <- grep("^b_", names(samples_s), value = TRUE)
  if (length(b_names) > 0) {
    for (bname in b_names) {
      coef_name <- sub("^b_", "", bname)
      # skip coefficients that are represented among ranef terms
      if (coef_name %in% meta$ranef_terms) next
      val_b <- as.numeric(samples_s[[bname]])
      if (!is.finite(val_b)) return(rep(-Inf, nrow(nodes)))
      if (coef_name == "Intercept") {
        fixed_part <- fixed_part + val_b
      } else {
        # covariate must exist in original data
        if (!(coef_name %in% names(fit$data))) return(rep(-Inf, nrow(nodes)))
        fixed_part <- fixed_part + val_b * fit$data[[coef_name]][idx_i]
      }
    }
  }

  Q <- nrow(nodes)
  mu_mat <- matrix(NA_real_, nrow = Q, ncol = nrow(Zi))
  for (q in seq_len(Q)) mu_mat[q, ] <- fixed_part + as.vector(Zi %*% nodes[q, ])
  
  # Family-specific conditional likelihood
  if (meta$family == "gaussian") {
    sigma <- if (!is.null(samples_s[["sigma"]])) as.numeric(samples_s[["sigma"]]) else stop("sigma missing")
    if (!is.finite(sigma) || sigma <= 0) return(rep(-Inf, nrow(nodes)))
    ll_cond <- loglik_cond_gaussian(y_i, mu_mat, sigma)
  } else if (meta$family == "bernoulli") {
    ll_cond <- loglik_cond_bernoulli(y_i, mu_mat)
  } else if (meta$family == "poisson") {
    ll_cond <- loglik_cond_poisson(y_i, mu_mat)
  } else {
    stop("Unsupported family: ", meta$family)
  }
  ll_theta <- tryCatch(mvtnorm::dmvnorm(nodes, mean = b_center, sigma = Sigma_u, log = TRUE), 
                       error = function(e) rep(-Inf, nrow(nodes)))
  ll_cond + ll_theta
}

#' Log-prior under brms parameterization (b_, sd_, sigma, cor_)
#'
#' Purpose:
#' - Evaluate the joint log-prior for one posterior draw by parsing
#'   brms::prior_summary and applying appropriate distributions with bounds.
#'
#' Inputs:
#' - samples_s: named list-like draw
#' - fit, meta: provide prior_summary, ranef structure, and family info
#'
#' Output: scalar log-prior; returns -Inf if any component invalid.
log_prior <- function(samples_s, fit, meta, ...) {
  ps <- meta$prior_summary
  lp_total <- 0
  for (bname in grep("^b_", names(samples_s), value = TRUE)) {
    coef_name <- sub("^b_", "", bname)
    row <- find_prior_row(ps, "b", coef_name, NULL)
    if (is.null(row) || !nzchar(as.character(row$prior[1])) ) next
    pr <- parse_prior(row$prior[1])
    if (is.na(pr$dist)) next
    val_b <- as.numeric(samples_s[[bname]])
    if (!is.finite(val_b)) return(-Inf)
    lb <- safe_bound(row$lb[1], -Inf)
    ub <- safe_bound(row$ub[1], Inf)
    contrib <- log_dtrunc(val_b, pr$dist, pr$pars, lb, ub)
    if (!is.finite(contrib)) return(-Inf)
    lp_total <- lp_total + contrib
  }
  for (k in seq_len(meta$K)) {
    nm <- paste0("sd_", meta$group, "__", meta$ranef_terms[k])
    row <- find_prior_row(ps, "sd", meta$ranef_terms[k], meta$group)
    if (is.null(row) || !nzchar(as.character(row$prior[1]))) next
    pr <- parse_prior(row$prior[1])
    if (is.na(pr$dist)) next
    val_sd <- as.numeric(samples_s[[nm]])
    if (!is.finite(val_sd) || val_sd < 0) return(-Inf)
    lb <- safe_bound(row$lb[1], 0)
    ub <- safe_bound(row$ub[1], Inf)
    contrib <- log_dtrunc(val_sd, pr$dist, pr$pars, lb, ub)
    if (!is.finite(contrib)) return(-Inf)
    lp_total <- lp_total + contrib
  }
  if (meta$family == "gaussian") {
    row <- find_prior_row(ps, "sigma", NULL, NULL)
    if (!is.null(row) && nzchar(as.character(row$prior[1]))) {
      pr <- parse_prior(row$prior[1])
      if (!is.na(pr$dist)) {
        val_sig <- as.numeric(samples_s[["sigma"]])
        if (!is.finite(val_sig) || val_sig <= 0) return(-Inf)
        lb <- safe_bound(row$lb[1], 0)
        ub <- safe_bound(row$ub[1], Inf)
        contrib <- log_dtrunc(val_sig, pr$dist, pr$pars, lb, ub)
        if (!is.finite(contrib)) return(-Inf)
        lp_total <- lp_total + contrib
      }
    }
  }
  cor_row <- NULL
  cor_row <- find_prior_row(ps, "cor", NULL, meta$group)
  if (is.null(cor_row)) cor_row <- find_prior_row(ps, "L", NULL, meta$group)
  if (!is.null(cor_row) && nzchar(as.character(cor_row$prior[1])) && grepl("lkj", tolower(as.character(cor_row$prior[1])))) {
    pr <- parse_prior(as.character(cor_row$prior[1]))
    if (!is.na(pr$dist)) {
      eta <- pr$pars[1]
      R <- diag(1, meta$K)
      for (i in seq_len(meta$K - 1)) for (j in (i + 1):meta$K) {
        nm <- paste0("cor_", meta$group, "__", meta$ranef_terms[i], "__", meta$ranef_terms[j])
        if (!is.null(samples_s[[nm]])) {
          rij <- as.numeric(samples_s[[nm]])
          if (!is.finite(rij)) return(-Inf)
          R[i, j] <- R[j, i] <- rij
        }
      }
      lkjv <- lkj_log_density(R, eta)
      if (!is.finite(lkjv)) return(-Inf)
      lp_total <- lp_total + lkjv
    }
  }
  lp_total
}

#' Run bleval ICs and fully marginal likelihood from a brms fit
#'
#' Applicability:
#' - Two-level models with a single grouping factor (one group variable).
#' - Supported family-link combinations:
#'   * Gaussian with identity link (Gaussian linear mixed-effects models).
#'   * Bernoulli with logit link (logistic mixed-effects models).
#'   * Poisson with log link (Poisson mixed-effects models).
#' - Single (univariate) response variable.
#' - Any number K of random effects under that grouping factor (e.g.,
#'   random intercept only, random intercept + one or more random slopes, or
#'   random-slope-only models without a random intercept).
#'
#' The function auto-extracts all required inputs for bleval from a brms fit:
#' samples, lv_mu, lv_cov, log_joint_i, meta. Users provide data (as a list) and
#' optional computational controls. Priors are taken from fit$prior using the
#' same parameterization as brms (b_\*, sd_\*, sigma, cor_\* with LKJ).
#'
#' @param fit brmsfit; fitted brms model.
#' @param data list; named list with elements:
#'   - Nobs: integer, total number of observations.
#'   - subject: integer vector (length Nobs), group index in 1..N.
#'   - N: integer, number of groups.
#'   - response column (e.g., y) and any covariates referenced by random-effect terms.
#'   This is the same structure you used for bleval examples.
#' @param Ngrid integer; number of quadrature nodes per dimension (default 3).
#' @param parallel logical; whether to parallelize per-unit quadrature (default FALSE).
#' @param n_cores integer; cores for parallel execution (default 1).
#' @param lb,ub optional named numeric vectors of lower/upper bounds for bridge sampling;
#'   if NULL, sensible defaults are derived from the parameter names.
#' @return list with components:
#'   - lm_res: result from bleval::log_marglik
#'   - ic_res: result from bleval::calc_IC
#'   - fml_res: result from bleval::log_fmarglik
#'   - meta: metadata list (family, group, ranef terms, etc.)
#'   - samples: posterior samples matrix used by bleval
#'   - lv_mu, lv_cov: lists used for quadrature centers and covariances
#'   - controls: list of settings used (Ngrid, parallel, n_cores)
bleval_with_brmsfit <- function(
  fit,
  data,
  Ngrid = 3,
  parallel = FALSE,
  n_cores = 1,
  lb = NULL,
  ub = NULL
) {
  # Family-link combination check
  fam <- family(fit)
  supported_combinations <- list(gaussian = "identity", bernoulli = "logit", poisson = "log")
  if (!fam$family %in% names(supported_combinations)) {
    stop("Unsupported family: ", fam$family, ". Supported families: ",
         paste(names(supported_combinations), collapse = ", "), ".")
  }
  expected_link <- supported_combinations[[fam$family]]
  if (fam$link != expected_link) {
    stop("The implementation of family '", fam$family, "' assumes link function '",
         expected_link, "', but found '", fam$link, "'.")
  }
  # Extract meta information (group, ranef terms, prior summary, etc.)
  meta <- prepare_brms_meta(fit)
  # Build lv_mu/lv_cov from posterior draws of coef(fit)
  lv_parts <- build_lv_from_fit(fit, meta)
  lv_mu <- lv_parts$lv_mu
  lv_cov <- lv_parts$lv_cov
  # Build posterior samples matrix
  samps_mat <- build_samples_matrix(fit, meta)
  # Default bounds for bridge sampling
  if (is.null(lb) || is.null(ub)) {
    bnd <- default_bounds(samps_mat, meta)
    if (is.null(lb)) lb <- bnd$lb
    if (is.null(ub)) ub <- bnd$ub
  }
  helper_exports_list <- helper_exports()
  # Run bleval log marginal likelihood (integrated out latent variables)
  lm_res <- bleval::log_marglik(
    samples = samps_mat,
    data = data,
    Ngrid = Ngrid,
    lv_mu = lv_mu,
    lv_cov = lv_cov,
    log_joint_i = log_joint_i,
    parallel = parallel,
    n_cores = n_cores,
    exports = helper_exports_list,
    fit = fit,
    meta = meta
  )
  # Information criteria
  ic_res <- bleval::calc_IC(lm_res, p_dic_version = 1)
  # Fully marginal likelihood via bridge sampling
  fml_res <- bleval::log_fmarglik(
    samples = samps_mat,
    data = data,
    Ngrid = Ngrid,
    lv_mu = lv_mu,
    lv_cov = lv_cov,
    log_joint_i = log_joint_i,
    log_prior = function(samples_s) log_prior(samples_s, fit = fit, meta = meta),
    fit = fit,
    meta = meta,
    lb = lb,
    ub = ub
  )
  list(
    lm_res = lm_res,
    ic_res = ic_res,
    fml_res = fml_res,
    meta = meta,
    samples = samps_mat,
    lv_mu = lv_mu,
    lv_cov = lv_cov,
    controls = list(Ngrid = Ngrid, parallel = parallel, n_cores = n_cores)
  )
}

#' Get list of helper functions for parallel export
#'
#' This function returns a character vector of all helper function names
#' that need to be exported to parallel workers, with optional filtering.
helper_exports <- function(){
  c("prepare_brms_meta", "build_Sigma_u_from_sample", "get_design_for_unit",
    "is_pos_def", "build_lv_from_fit", "build_samples_matrix", "default_bounds", 
    "loglik_cond_gaussian", "loglik_cond_bernoulli", "loglik_cond_poisson", 
    "parse_prior", "safe_bound", "find_prior_row", "log_dtrunc", "lkj_log_density", 
    "log_joint_i", "log_prior")
}
