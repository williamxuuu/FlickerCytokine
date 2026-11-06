# ============================================================================
# 05_correlations.R  —  cytokine ~ amyloid and cytokine ~ morphology.
#
# Spearman correlation of every panel analyte with each target, within
# genotype, on the sex-z-scored table (both sides of every pair, since several
# targets -- e.g. astrocyte cell number -- are themselves sex-dimorphic).
#
# PERMUTATION P-VALUES. Permuting the target against the analytes holds both
# marginal distributions fixed, so the null distribution reflects the exact tie
# structure of the data. This avoids cor.test()'s asymptotic approximation,
# which is unreliable with ties at n ~ 20. B = 50,000.
#
# FREE VS SEX-STRATIFIED PERMUTATION. Because sex has already been removed by
# the within-sex z-score, permutations are unrestricted. A sex-stratified
# permutation (labels shuffled only within sex) is computed as a sensitivity
# analysis and reported as `p_perm_sexstrat`.
#
# MULTIPLE TESTING. A family is one genotype x one target (28 analytes). The
# same permutations are applied to every analyte in a family, preserving the
# correlation structure among cytokines. Two corrections are reported:
#   p_adj  Benjamini-Hochberg FDR on the permutation p-values (primary).
#   p_wy   Westfall-Young single-step maxT family-wise error rate, computed on
#          studentised rho (each analyte's |rho| divided by its null SD), so
#          that analytes with many ties -- whose null distribution of rho is
#          wider -- do not dominate the maximum.
# Pairs that pass BH but not Westfall-Young are treated as weaker evidence and
# shown in supplementary panels; low-resolution analytes are a common cause.
#
# RESOLUTION. `n_distinct` (number of distinct values) and `prop_tied`
# (fraction of tied pairs) are computed on the raw MFI, since measurement
# resolution is a property of the assay rather than of the transform.
#
# Outputs: stats/corr_cytokine_amyloid.csv, stats/corr_cytokine_morphology.csv.
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
source(file.path(.here(), "00_setup.R"))

B <- 50000
set.seed(42)

free_perm <- function(n, B) replicate(B, sample.int(n))
strat_perm <- function(sex, B) {
  n <- length(sex); out <- matrix(seq_len(n), n, B)
  for (g in split(seq_len(n), sex)) if (length(g) > 1)
    out[g, ] <- apply(out[g, , drop = FALSE], 2, sample)
  out
}
scale_cols <- function(M) {
  M <- sweep(M, 2, colMeans(M), "-")
  s <- sqrt(colSums(M^2)); s[s == 0] <- NA_real_
  sweep(M, 2, s, "/")
}

one_family <- function(d, raw, target, subgroup) {
  ok <- !is.na(d[[target]]) & complete.cases(d[, PANEL, drop = FALSE])
  d <- d[ok, , drop = FALSE]; raw <- raw[ok, , drop = FALSE]; n <- nrow(d)
  if (n < 5) return(NULL)

  Rx  <- scale_cols(sapply(PANEL, function(a) rank(d[[a]])))
  ry  <- rank(d[[target]])
  obs <- as.vector(crossprod(Rx, scale_cols(matrix(ry, ncol = 1))))

  ND  <- crossprod(Rx, scale_cols(matrix(ry[free_perm(n, B)],                     nrow = n)))
  NDs <- crossprod(Rx, scale_cols(matrix(ry[strat_perm(as.character(d$Sex), B)], nrow = n)))

  p_perm <- (1 + rowSums(abs(ND)  >= abs(obs) - 1e-12)) / (B + 1)
  p_strat<- (1 + rowSums(abs(NDs) >= abs(obs) - 1e-12)) / (B + 1)

  sdj <- apply(ND, 1, sd); sdj[sdj == 0] <- NA_real_
  maxZ <- apply(abs(sweep(ND, 1, sdj, "/")), 2, max, na.rm = TRUE)
  p_wy <- (1 + sapply(abs(obs) / sdj, function(z) sum(maxZ >= z - 1e-12))) / (B + 1)

  data.frame(subgroup, target, analyte = PANEL, rho = obs, n = n,
             p_perm, p_wy, p_perm_sexstrat = p_strat,
             # resolution measured on the RAW assay values -- see header
             n_distinct = sapply(PANEL, function(a) length(unique(raw[[a]]))),
             prop_tied  = sapply(PANEL, function(a) {
               k <- table(raw[[a]]); sum(k * (k - 1) / 2) / (n * (n - 1) / 2) }),
             row.names = NULL, stringsAsFactors = FALSE)
}

run <- function(targets, genos, label) {
  res <- list()
  for (g in genos) for (tg in targets) {
    r <- one_family(ANIMALS    [ANIMALS$Genotype     == g, , drop = FALSE],
                    ANIMALS_RAW[ANIMALS_RAW$Genotype == g, , drop = FALSE], tg, g)
    if (!is.null(r)) res[[paste(g, tg)]] <- r
  }
  out <- do.call(rbind, res)
  out$p_adj <- ave(out$p_perm, paste(out$subgroup, out$target),
                   FUN = function(p) p.adjust(p, "BH"))
  out <- out[order(out$subgroup, out$target, out$p_perm), ]
  write.csv(out, file.path(STATS, paste0("corr_", label, ".csv")), row.names = FALSE)
  cat(sprintf("[05] %-22s %4d pairs | BH q<.05: %2d | WY p<.05: %d\n",
              label, nrow(out), sum(out$p_adj < .05), sum(out$p_wy < .05)))
  d <- sum(xor(out$p_perm < .05, out$p_perm_sexstrat < .05))
  cat(sprintf("[05]   free vs sex-stratified null disagree at p<.05 on %d of %d pairs\n",
              d, nrow(out)))
  out
}

invisible(run(AMYLOID_TARGETS, "5XFAD", "cytokine_amyloid"))   # WT has no amyloid
invisible(run(MORPH_TARGETS, GENOTYPES, "cytokine_morphology"))
cat("[05] ->", STATS, "\n")
