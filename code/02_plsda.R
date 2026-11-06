# ============================================================================
# 02_plsda.R  —  PLS-DA per genotype, model validation, and key cytokines.
#
# Fits a two-component PLS-DA (ropls; unit-variance scaling) of flicker group on
# the sex-z-scored 28-analyte panel, separately for each genotype, and defines
#
#   key cytokines = top 3 analytes by |loading| on LV1  +  top 3 on LV2
#
# which the downstream figure stages use for analyte selection. The model is
# fitted only here; figure stages read the cached result, so every panel uses
# the same fit as the reported statistics.
#
# MODEL VALIDATION. Leave-one-out cross-validated (LOOCV) classification
# accuracy is compared with a null distribution from 1,000 permutations of the
# group labels: p = (#{permuted accuracy >= observed} + 1) / (B + 1).
# Loading uncertainty is the SD of each loading across 1,000 refits, each
# leaving out one randomly chosen animal (loadings are sign-aligned to the full
# fit first, since PLS components are defined only up to sign).
#
# NO HYPOTHESIS TESTS ARE RUN ON LATENT-VARIABLE SCORES. The latent variables
# are constructed to maximise separation between the groups, so testing group
# differences on them is circular. Score plots are descriptive; inference on
# the model as a whole comes from the permutation test above.
#
# CACHING. The permutation and LOOCV steps take ~4 min per genotype, so they are
# cached in stats/plsda_<genotype>.rds, keyed on a fingerprint of the input
# matrix. Any change to the data, panel or transform invalidates the cache.
# Delete the .rds files to force recomputation.
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
source(file.path(.here(), "00_setup.R"))
suppressPackageStartupMessages({ library(ropls); library(matrixStats) })

N_PERM  <- 1000      # label permutations
N_LOOCV <- 1000      # leave-one-out resamples for the loading error bars
N_KEY   <- 3         # top-N per latent variable -> the key cytokines
set.seed(42)

fit <- function(X, y) suppressWarnings(opls(X, y, predI = 2, permI = 0,
                       scaleC = "standard", fig.pdfC = "none", info.txtC = "none"))

loocv_acc <- function(X, y) {
  n <- nrow(X)
  mean(vapply(seq_len(n), function(i) {
    ok <- colSds(X[-i, , drop = FALSE]) > 1e-10
    m  <- tryCatch(fit(X[-i, ok, drop = FALSE], y[-i]), error = function(e) NULL)
    if (is.null(m)) return(NA)
    p <- tryCatch(as.character(predict(m, newdata = X[i, ok, drop = FALSE])),
                  error = function(e) NA_character_)
    if (is.na(p)) NA else p == y[i]
  }, logical(1)), na.rm = TRUE)
}

# Loading stability. Each resample's loadings are sign-aligned to the full-data
# loadings before the SD is taken -- PLS-DA components are only defined up to
# sign, so unaligned resamples would cancel and report a spuriously huge SD.
loocv_loadings <- function(X, y, Lref, runs = N_LOOCV) {
  n <- nrow(X)
  P <- array(NA_real_, c(runs, ncol(X), 2), list(NULL, colnames(X), c("LV1","LV2")))
  for (i in seq_len(runs)) {
    out <- sample.int(n, 1); Xi <- X[-out, , drop = FALSE]
    ok <- colSds(Xi) > 1e-10
    m  <- tryCatch(fit(Xi[, ok, drop = FALSE], y[-out]), error = function(e) NULL)
    if (is.null(m)) next
    for (k in 1:2) {
      v <- rep(0, ncol(X)); v[ok] <- m@loadingMN[, k]
      P[i, , k] <- v * sign(sum(v * Lref[, k]))
    }
  }
  cbind(LV1 = colSds(P[,,1], na.rm = TRUE), LV2 = colSds(P[,,2], na.rm = TRUE))
}

# Model matrix for one genotype, from the sex-z-scored table.
plsda_input <- function(g) {
  d <- ANIMALS[ANIMALS$Genotype == g, , drop = FALSE]
  d <- d[order(d$Flicker, d$Sex), , drop = FALSE]
  X <- as.matrix(d[, PANEL]); X <- X[, apply(X, 2, sd) > 0, drop = FALSE]
  list(d = d, X = X, y = as.character(d$Flicker))
}

for (g in GENOTYPES) {
  io <- plsda_input(g); X <- io$X; y <- io$y
  m  <- fit(X, y); L <- m@loadingMN

  # Orient each LV so its largest loading is positive, so the sign of "LV1" is
  # reproducible across runs and the two genotypes' panels read the same way.
  sgn <- vapply(1:2, function(k) if (L[which.max(abs(L[, k])), k] < 0) -1 else 1, numeric(1))

  # --- key cytokines --------------------------------------------------------
  key_by_lv <- lapply(1:2, function(k) names(sort(abs(L[, k]), decreasing = TRUE))[seq_len(N_KEY)])
  key <- unique(unlist(key_by_lv))

  fingerprint <- paste(
    g, ncol(X), nrow(X), sprintf("%.10f", sum(X)), sprintf("%.10f", sum(X^2)),
    paste(colnames(X), collapse = ","), sep = "|")

  cache <- file.path(STATS, sprintf("plsda_%s.rds", g))
  cc <- if (file.exists(cache)) readRDS(cache) else NULL
  if (!is.null(cc) && identical(cc$fingerprint, fingerprint)) {
    cat("[02]", g, "- cache valid, reusing permutation + LOOCV\n")
  } else {
    if (!is.null(cc)) cat("[02]", g, "- input changed, recomputing (~4 min)\n")
    else              cat("[02]", g, "- no cache, computing (~4 min)\n")
    obs   <- loocv_acc(X, y)
    null  <- replicate(N_PERM, loocv_acc(X, sample(y)))
    pperm <- (sum(null >= obs, na.rm = TRUE) + 1) / (sum(!is.na(null)) + 1)
    cc <- list(fingerprint = fingerprint, obs = obs, pperm = pperm,
               sdL = loocv_loadings(X, y, L))
  }

  # Refreshed every run; only the expensive resampling is cached.
  cc$loadings  <- L
  cc$sign      <- sgn
  cc$scores    <- cbind(LV1 = m@scoreMN[, 1] * sgn[1], LV2 = m@scoreMN[, 2] * sgn[2])
  cc$flicker   <- y
  cc$sex       <- as.character(io$d$Sex)
  cc$r2x       <- m@modelDF$R2X
  cc$key       <- key
  cc$key_by_lv <- key_by_lv
  cc$analytes  <- colnames(X)
  saveRDS(cc, cache)

  cat(sprintf("[02] %-6s LOOCV %.3f vs %.3f chance | permutation p = %.3f%s\n",
              g, cc$obs, 1/3, cc$pperm,
              if (cc$pperm >= 0.05) "  (not significant)" else ""))
  cat(sprintf("[02]   key LV1: %s\n", paste(key_by_lv[[1]], collapse = ", ")))
  cat(sprintf("[02]   key LV2: %s\n", paste(key_by_lv[[2]], collapse = ", ")))
}

# --- exported for downstream stages -----------------------------------------
# Downstream stages use KEY[[g]] / read_plsda(g) rather than refitting.
read_plsda <- function(g) readRDS(file.path(STATS, sprintf("plsda_%s.rds", g)))
KEY <- lapply(setNames(GENOTYPES, GENOTYPES), function(g) read_plsda(g)$key)

kt <- do.call(rbind, lapply(GENOTYPES, function(g) {
  cc <- read_plsda(g)
  data.frame(genotype = g,
             lv = rep(c("LV1", "LV2"), each = N_KEY),
             rank = rep(seq_len(N_KEY), 2),
             analyte = unlist(cc$key_by_lv),
             loading = c(cc$loadings[cc$key_by_lv[[1]], 1] * cc$sign[1],
                         cc$loadings[cc$key_by_lv[[2]], 2] * cc$sign[2]))
}))
write.csv(kt, file.path(STATS, "plsda_key_cytokines.csv"), row.names = FALSE)
cat("[02] key cytokines -> stats/plsda_key_cytokines.csv\n")
