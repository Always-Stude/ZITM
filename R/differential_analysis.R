#' Differential Expression Analysis Based on Zero-Inflated Telegraph Model
#'
#' Performs differential expression analysis for single-cell RNA-seq data
#' using Zero-Inflated Telegraph Model.
#'
#' @param counts A count matrix (genes x cells)
#' @param group Factor vector specifying sample groups
#' @param parallel Logical, whether to use parallel computation
#' @param seed_value Random seed
#' @param repetitive Logical, Set seed to ensure reproducibility of parameter estimation
#' @param BPPARAM BiocParallel parameters
#' @return A data frame with differential analysis results. Columns:
#' \itemize{
#'   \item{PVAL}{: Raw p-value of differential expression}
#'   \item{Ex_group1}{: Mean expression in group 1}
#'   \item{Ex_group2}{: Mean expression in group 2}
#'   \item{Var_group1}{: Expression variance in group 1}
#'   \item{Var_group2}{: Expression variance in group 2}
#'   \item{padj}{: Adjusted p-value (FDR correction)}
#' }
#' @examples
#' # Create a simulated count matrix
#' mat <- matrix(
#'   sample(0:50, size = 100*100, replace = TRUE),
#'   nrow = 100,
#'   ncol = 100,
#'   dimnames = list(
#'     paste0("gene", 1:100),
#'     paste0("cell", 1:100)
#'   )
#' )
#'
#' # Define grouping information
#' group <- factor(c(rep(1,50), rep(2,50)))
#'
#' # Run ZITM differential analysis
#' results <- Differential_analysis(counts=mat,group=group)
#'
#' # View the first 6 rows of results
#' head(results)
#' @export
#' @importFrom pbapply pblapply
#' @importFrom Matrix Matrix
#' @importFrom MASS glm.nb
#' @importFrom dplyr %>%
#' @importFrom statmod gauss.quad
#' @importFrom stats glm p.adjust pchisq optim wilcox.test median ppois pnorm make.link
#' @importFrom parallel detectCores
#' @importFrom BiocParallel bplapply bpparam

Differential_analysis <- function(counts, group, parallel = FALSE, seed_value = 123, repetitive = FALSE, BPPARAM = bpparam()) {
  if (sum(is.na(counts)) > 0)
    stop("NA detected in 'counts'")

  if (sum(counts < 0) > 0)
    stop("Negative value detected in 'counts'")

  if (all(counts == 0))
    stop("All elements of 'counts' are zero")

  if (any(colSums(counts) == 0))
    warning("Library size of zero detected in 'counts'")

  if (!is.factor(group))
    stop("Data type of 'group' is not factor")

  counts <- round(as.matrix(counts))
  storage.mode(counts) <- "integer"
  if (any(rowSums(counts) == 0))
    message("Removing ", sum(rowSums(counts) == 0), " rows of genes with all zero counts")
  counts <- counts[rowSums(counts) != 0, ]
  geneNum <- nrow(counts)
  sampleNum <- ncol(counts)
  gc()

  message("Normalizing the data")
  GEOmean <- rep(NA, geneNum)
  for (i in 1:geneNum) {
    gene_NZ <- counts[i, counts[i, ] > 0]
    GEOmean[i] <- exp(sum(log(gene_NZ), na.rm = TRUE)/length(gene_NZ))
  }
  S <- rep(NA, sampleNum)
  counts_norm <- counts
  for (j in 1:sampleNum) {
    sample_j <- counts[, j]/GEOmean
    S[j] <- median(sample_j[which(sample_j != 0)])
    counts_norm[, j] <- counts[, j]/S[j]
  }
  norm_count3 <- ceiling(counts_norm)
  norm_count3 <- Matrix(norm_count3, sparse = TRUE)
  remove(GEOmean, gene_NZ, S, sample_j, i, j, geneNum, sampleNum, counts, counts_norm)
  gc()

  unique_conditions <- unique(as.character(group))
  cat("Normalization completed!\n")
  N3 <- ncol(norm_count3)

  Cell_DE <- function(i){
    Calculate_gene_frequency <- function(gene_expression) {
      if (!is.matrix(gene_expression) && !is.data.frame(gene_expression) && !inherits(gene_expression, "dgCMatrix")) {
        stop("Input must be a matrix, data frame, or dgCMatrix")
      }

      if (nrow(gene_expression) != 1) {
        stop("Input must have exactly one row")
      }

      expression_values <- as.numeric(gene_expression)
      max_value <- max(expression_values)
      freq_table <- table(factor(expression_values, levels = 0:max_value))
      freq_table <- freq_table / sum(freq_table)

      freq_df <- data.frame(matrix(NA, nrow = 1, ncol = max_value + 1))
      colnames(freq_df) <- paste0("expr", 0:max_value)
      freq_df[1, as.numeric(names(freq_table)) + 1] <- freq_table
      freq_df[is.na(freq_df)] <- 0
      row.names(freq_df) <- row.names(gene_expression)
      return(freq_df)
    }

    BP_prob <- function (n = NULL, kon = NULL, koff = NULL, ksyn = NULL, theta = NULL){
      fn3 <- function(x1, x2, m) {
        if (max(m) < 1e+05)
          res <- ppois(x2, m) - ppois(x1, m)
        else res <- pnorm(x2, m, sqrt(m)) - pnorm(x1, m, sqrt(m))
        return(res)
      }

      w <- gauss.quad(10, "jacobi", alpha = koff - 1, beta = kon - 1)
      gs <- sum(w$weight * fn3(x1 = n-1, x2 = n, m = ksyn * (1 + w$node)/2))

      if(n == 0){
        prob <- theta + 1/beta(kon, koff) * 2^(-kon - koff + 1) * gs*(1-theta)
      } else {
        prob <- 1/beta(kon, koff) * 2^(-kon - koff + 1) * gs*(1 - theta)
      }

      return(prob)
    }

    Zero_BPfam <- function(kon = NULL, koff = NULL, ksyn = NULL, theta = NULL, link = "log"){
      linktemp <- substitute(link)
      if (!is.character(linktemp))
        linktemp <- deparse(linktemp)
      okLinks <- c("log", "identity", "sqrt")
      if (linktemp %in% okLinks)
        stats <- make.link(link)
      else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
      } else {
        if (inherits(link, "link-glm")) {
          stats <- link
          if (!is.null(stats$name))
            linktemp <- stats$name
        } else {
          stop(gettextf("link \"%s\" not available for betaPoisson family; available links are %s",
                        linktemp, paste(sQuote(okLinks), collapse = ", ")),
               domain = NA)
        }
      }

      initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the 'betaPoisson' family")
        n <- rep.int(1, nobs)
        mustart <- y + (y == 0)/1e+09
      })

      aic <- function(y, n, mu, wt, dev) NA
      dev.resids <- function(y, mu, wt) {
        r <- mu * wt
        p <- which(y > 0)
        r[p] <- (wt * (y * log(y/mu) - (y - mu)))[p]
        2 * r
      }

      validmu <- function(mu) all(is.finite(mu)) && all(mu > 0)
      famname <- paste("Zero_Beta Poisson(", format(round(kon, 4)),
                       format(round(koff, 4)), format(round(ksyn, 4)),
                       format(round(theta, 4)), ")", sep = "")
      chiq2 <- koff/(kon*(kon+koff+1))
      variance <- function(mu) {
        mu + mu*mu*chiq2
      }

      stats <- make.link(link)
      structure(list(family = famname, link = link, linkfun = stats$linkfun,
                     aic = aic, dev.resids = dev.resids, linkinv = stats$linkinv,
                     variance = variance, mu.eta = stats$mu.eta, validmu = validmu,
                     initialize = initialize, valideta = stats$valideta),
                class = "family")
    }

    Zero_BP <- function(param = NULL, data = NULL, pmdata = NULL, m = NULL){
      ksyn <- exp(param[3]) + 1
      kon <- exp(param[1])
      koff <- exp(param[2])
      theta <- exp(param[4])

      cellnumber <- m*t(pmdata)

      k <- length(pmdata)
      for (i in 1:length(pmdata)) {
        if (pmdata[k] > 0) {
          break
        } else {
          k <- length(pmdata) - i
        }
      }
      N <- 3 * (k + 1)

      dist <- sapply(0:max(data), function(n) {BP_prob(kon = kon, koff = koff, ksyn = ksyn, n = n, theta = theta)})
      dist <- dist + 1e-10
      NN <- min(length(pmdata), N)
      dist_b <- dist[1:NN]

      S1 <- 0
      if (any(is.infinite(dist_b), na.rm = TRUE) || any(is.na(dist_b)) || any(dist_b < 0, na.rm = TRUE) || (sum(dist_b) < 0.95) || (sum(dist_b) > 1.05)) {
        S1 <- 1e+30
      } else {
        for (i in 1:length(dist_b)) {
          ratio <- pmdata[i] / dist_b[i]
          if (pmdata[i] == 0) {
            next
          }
          if (ratio > 0) {
            S1 <- S1 + cellnumber[i] * log(ratio)
          }
        }
      }
      if(theta > 1) S1 <- 1e+40
      if(kon > 1000) S1 <- 1e+40
      if(koff > 1000) S1 <- 1e+40
      return(S1)
    }

    Zero_BPest <- function(name = NULL, data = NULL, N = NULL, count_freq = NULL, seed_value = NULL, repetitive = NULL){
      pmdata <- t(count_freq)

      kon_e <- numeric(0)
      koff_e <- numeric(0)
      ksyn_e <- numeric(0)
      theta_e <- numeric(0)
      Smin_TM <- numeric(0)

      if (repetitive) set.seed(seed_value)

      n_intervals <- 10
      param_range <- c(0.1, 5)
      ksyn_range <- c(5, 50)

      intervals <- seq(param_range[1], param_range[2], length.out = n_intervals + 1)
      ksyn_intervals <- seq(ksyn_range[1], ksyn_range[2], length.out = n_intervals + 1)

      for (i in 1:n_intervals) {
        kon_interval <- sample(1:n_intervals, 1)
        kon <- runif(1, intervals[kon_interval], intervals[kon_interval+1])

        available_intervals <- setdiff(1:n_intervals, kon_interval)
        koff_interval <- sample(available_intervals, 1)
        koff <- runif(1, intervals[koff_interval], intervals[koff_interval+1])

        available_intervals <- setdiff(1:n_intervals, c(kon_interval, koff_interval))
        ksyn_interval <- sample(available_intervals, 1)
        ksyn <- runif(1, ksyn_intervals[ksyn_interval], ksyn_intervals[ksyn_interval+1])

        theta <- runif(1, min = max(0, pmdata[1,1]-0.0001), max = min(pmdata[1,1]+0.0001, 1))
        param <- log(c(kon, koff, ksyn, theta))
        result <- optim(param, Zero_BP, m = N, pmdata = pmdata, data = data)
        bminTM <- result$par
        SminTM <- result$value

        kon_e <- c(kon_e, exp(bminTM[1]))
        koff_e <- c(koff_e, exp(bminTM[2]))
        ksyn_e <- c(ksyn_e, exp(bminTM[3]) + 1)
        theta_e <- c(theta_e, exp(bminTM[4]))
        Smin_TM <- c(Smin_TM, SminTM)
      }

      Smin_minTM_location <- which(Smin_TM == min(Smin_TM))
      kon <- kon_e[Smin_minTM_location[1]]
      koff <- koff_e[Smin_minTM_location[1]]
      ksyn <- ksyn_e[Smin_minTM_location[1]]
      theta <- theta_e[Smin_minTM_location[1]]

      log_probs <- sapply(data, function(n) {
        prob <- BP_prob(n = n, kon = kon, koff = koff, ksyn = ksyn, theta = theta)
        if (prob > 0) return(log(prob)) else return(NA)
      })
      L <- sum(log_probs, na.rm = TRUE)
      EX <- sum(count_freq * (0:(length(count_freq) - 1)), na.rm = TRUE)
      EX2 <- sum(count_freq * ((0:(length(count_freq) - 1))^2), na.rm = TRUE)
      Var <- EX2 - EX^2
      T <- min(Smin_TM)
      df <- length(count_freq) - 5
      pval <- ifelse(df > 0, pchisq(2*T, df = df, lower.tail = FALSE), 0.9999)

      if (name != 'Total'){
        result_df <- data.frame(kon = kon, koff = koff, ksyn = ksyn, theta = theta, L = L, EX = EX, Var = Var, pval = pval)
        colnames(result_df) <- paste0(colnames(result_df), "_", name)
      } else {
        result_df <- data.frame(kon_Total = kon, koff_Total = koff, ksyn_Total = ksyn, pval_Total = pval,
                                theta_Total = theta, L_Total = L, EX_Total = EX, Var_Total = Var)
      }
      return(result_df)
    }

    Nonzero_BP <- function(param = NULL, data = NULL, pmdata = NULL, m = NULL){
      ksyn <- exp(param[3]) + 1
      kon <- exp(param[1])
      koff <- exp(param[2])
      cellnumber <- m*t(pmdata)

      k <- length(pmdata)
      for (i in 1:length(pmdata)) {
        if (pmdata[k] > 0) {
          break
        } else {
          k <- length(pmdata) - i
        }
      }
      N <- 3 * (k + 1)

      dist <- sapply(0:max(data), function(n) {BP_prob(kon = kon, koff = koff, ksyn = ksyn, n = n, theta = 0)})
      dist <- dist + (1e-10)
      NN <- min(length(pmdata), N)
      dist_b <- dist[1:NN]

      S1 <- 0
      if (any(is.infinite(dist_b), na.rm = TRUE) || any(is.na(dist_b)) || any(dist_b < 0, na.rm = TRUE) || (sum(dist_b) < 0.95) || (sum(dist_b) > 1.05)) {
        S1 <- 1e+30
      } else {
        for (i in 1:length(dist_b)) {
          ratio <- pmdata[i] / dist_b[i]
          if (pmdata[i] == 0) {
            next
          }
          if (ratio > 0) {
            S1 <- S1 + cellnumber[i] * log(ratio)
          }
        }
      }
      if(kon > 1000) S1 <- 1e+40
      if(koff > 1000) S1 <- 1e+40
      return(S1)
    }

    Nonzero_BPest <- function(name = NULL, data = NULL, N = NULL, count_freq = NULL, seed_value = NULL, repetitive = NULL){
      pmdata <- t(count_freq)

      kon_e <- numeric(0)
      koff_e <- numeric(0)
      ksyn_e <- numeric(0)
      Smin_TM <- numeric(0)

      if (repetitive) set.seed(seed_value)

      n_intervals <- 10
      param_range <- c(0.1, 5)
      ksyn_range <- c(5, 50)

      intervals <- seq(param_range[1], param_range[2], length.out = n_intervals + 1)
      ksyn_intervals <- seq(ksyn_range[1], ksyn_range[2], length.out = n_intervals + 1)

      for (i in 1:n_intervals) {
        kon_interval <- sample(1:n_intervals, 1)
        kon <- runif(1, intervals[kon_interval], intervals[kon_interval+1])

        available_intervals <- setdiff(1:n_intervals, kon_interval)
        koff_interval <- sample(available_intervals, 1)
        koff <- runif(1, intervals[koff_interval], intervals[koff_interval+1])

        available_intervals <- setdiff(1:n_intervals, c(kon_interval, koff_interval))
        ksyn_interval <- sample(available_intervals, 1)
        ksyn <- runif(1, ksyn_intervals[ksyn_interval], ksyn_intervals[ksyn_interval+1])

        param <- log(c(kon, koff, ksyn))
        result <- optim(param, Nonzero_BP, m = N, pmdata = pmdata, data = data)
        bminTM <- result$par
        SminTM <- result$value

        kon_e <- c(kon_e, exp(bminTM[1]))
        koff_e <- c(koff_e, exp(bminTM[2]))
        ksyn_e <- c(ksyn_e, exp(bminTM[3]) + 1)
        Smin_TM <- c(Smin_TM, SminTM)
      }

      Smin_minTM_location <- which(Smin_TM == min(Smin_TM))
      kon <- kon_e[Smin_minTM_location[1]]
      koff <- koff_e[Smin_minTM_location[1]]
      ksyn <- ksyn_e[Smin_minTM_location[1]]

      log_probs <- sapply(data, function(n) {
        prob <- BP_prob(n = n, kon = kon, koff = koff, ksyn = ksyn, theta = 0)
        if (prob > 0) return(log(prob)) else return(NA)
      })
      L <- sum(log_probs, na.rm = TRUE)
      EX <- sum(count_freq * (0:(length(count_freq) - 1)), na.rm = TRUE)
      EX2 <- sum(count_freq * ((0:(length(count_freq) - 1))^2), na.rm = TRUE)
      Var <- EX2 - EX^2
      T <- min(Smin_TM)
      df <- length(count_freq) - 4
      pval <- ifelse(df > 0, pchisq(2*T, df = df, lower.tail = FALSE), 0.999)

      if (name != 'Total'){
        result_df <- data.frame(kon = kon, koff = koff, ksyn = ksyn, theta = 0, L = L, EX = EX, Var = Var, pval = pval)
        colnames(result_df) <- paste0(colnames(result_df), "_", name)
      } else {
        result_df <- data.frame(kon_Total = kon, koff_Total = koff, ksyn_Total = ksyn, pval_Total = pval,
                                theta_Total = 0, L_Total = L, EX_Total = EX, Var_Total = Var)
      }
      return(result_df)
    }

    norm_count1 <- norm_count3[, group == levels(group)[1]]
    N1 <- ncol(norm_count1)
    norm_count2 <- norm_count3[, group == levels(group)[2]]
    N2 <- ncol(norm_count2)

    row1 <- as.matrix(norm_count1[i, , drop = FALSE])
    row_freq1 <- Calculate_gene_frequency(row1)
    Ex1 <- sum(row_freq1 * (0:(length(row_freq1) - 1)), na.rm = TRUE)
    VAR1 <- (sum(row_freq1 * ((0:(length(row_freq1) - 1))^2), na.rm = TRUE) - (sum(row_freq1 * (0:(length(row_freq1) - 1)), na.rm = TRUE)^2))

    row <- as.matrix(norm_count2[i, , drop = FALSE])
    row_freq <- Calculate_gene_frequency(row)
    Ex <- sum(row_freq * (0:(length(row_freq) - 1)), na.rm = TRUE)
    VAR <- (sum(row_freq * ((0:(length(row_freq) - 1))^2), na.rm = TRUE) - (sum(row_freq * (0:(length(row_freq) - 1)), na.rm = TRUE)^2))

    df1 <- tryCatch({
      if (row_freq[1,1] > 0 && row_freq[1,1] < 1) {
        Zero_BPest(data = row, name = unique_conditions[2],
                   N = N2, count_freq = row_freq, seed_value = seed_value, repetitive = repetitive)
      } else if (row_freq[1,1] == 0) {
        Nonzero_BPest(data = row, name = unique_conditions[2],
                      N = N2, count_freq = row_freq, seed_value = seed_value, repetitive = repetitive)
      } else {
        data.frame(kon = 0, koff = 1e+40, ksyn = 0, theta = 1, L = 0, EX = 0, Var = 0, pval = 1, check.names = FALSE) %>%
          setNames(paste0(c("kon","koff","ksyn","theta","L","EX","Var","pval"), "_", unique_conditions[2]))
      }
    }, error = function(e) {
      message(sprintf("Error in %s %s: %s", unique_conditions[2], row.names(row), e$message))
      return(data.frame(kon = NA, koff = NA, ksyn = NA, theta = NA, L = 1e+40, pval = 1,
                        EX = Ex, Var = VAR, check.names = FALSE) %>%
               setNames(paste0(c("kon", "koff", "ksyn", "theta", "L","EX","Var","pval"), "_", unique_conditions[2])))
    })
    row.names(df1) <- row.names(row)

    row_dense <- as.matrix(norm_count3[i, , drop = FALSE])
    fdata <- as.data.frame(t(row_dense))
    fdata$group <- ifelse(row.names(fdata) %in% colnames(norm_count1), unique_conditions[1], unique_conditions[2])
    fdata$group <- as.factor(fdata$group)
    colnames(fdata) <- c("expression", "group")
    i.pval <- NA
    i.tval <- NA
    i.converged <- NA

    if (!all(row == 0)){
      try({
        kon <- df1[[paste0("kon_", unique_conditions[2])]]
        koff <- df1[[paste0("koff_", unique_conditions[2])]]
        ksyn <- df1[[paste0("ksyn_", unique_conditions[2])]]
        theta <- df1[[paste0("theta_", unique_conditions[2])]]

        fam0 <- do.call("Zero_BPfam", list(kon = kon, koff = koff, ksyn = ksyn, theta = theta, link = "log"))
        fit <- glm(expression ~ ., data = fdata, family = fam0)

        i.pval <- summary(fit)$coefficients[2,4]
        i.tval <- summary(fit)$coefficients[2,3]
        i.converged <- fit$converged
        if (!i.converged) {
          fit <- glm(expression ~ ., data = fdata, family = quasipoisson)
          i.pval <- summary(fit)$coefficients[2,4]
          i.tval <- summary(fit)$coefficients[2,3]
          i.converged <- fit$converged
          if(!i.converged){
            gw_test <- wilcox.test(row, row1)
            i.pval <- gw_test$p.value
            i.tval <- gw_test$statistic
            i.converged <- NA
          }
        }
      }, silent = FALSE)
    } else {
      if(all(row == 0) && all(row1 == 0)){
        i.pval <- NA
        i.tval <- NA
        i.converged <- NA
      } else {
        gw_test <- wilcox.test(row, row1)
        i.pval <- gw_test$p.value
        i.tval <- gw_test$statistic
        i.converged <- NA
      }
    }
    return(data.frame(PVAL = i.pval, TVAL = i.tval, CONVERGED = i.converged, Ex_cond1 = Ex1, Ex_cond2 = Ex,
                      Var_cond1 = VAR1, Var_cond2 = VAR, p = df1[[paste0("pval_", unique_conditions[2])]],
                      row.names = rownames(norm_count3)[i]))
  }

  if(!parallel){
    results <- do.call(rbind, pblapply(1:nrow(norm_count3), Cell_DE))
  } else {
    results <- do.call(rbind, bplapply(1:nrow(norm_count3), Cell_DE, BPPARAM = BPPARAM))
  }
  colnames(results)[colnames(results) == "Ex_cond1"] <- paste0("Ex_", unique_conditions[1])
  colnames(results)[colnames(results) == "Ex_cond2"] <- paste0("Ex_", unique_conditions[2])
  colnames(results)[colnames(results) == "Var_cond1"] <- paste0("Var_", unique_conditions[1])
  colnames(results)[colnames(results) == "Var_cond2"] <- paste0("Var_", unique_conditions[2])

  results <- results[!is.na(results$PVAL), ]
  results$padj <- p.adjust(results$PVAL, method = 'fdr')
  results <- results[, -c(2, 3, 8)]
  return(results)
}
