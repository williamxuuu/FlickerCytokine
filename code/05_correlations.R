# ============================================================================
# 05_correlations.R  —  cytokine <-> amyloid and cytokine <-> morphology.
#
# INPUT: the sex-z-scored table (00_setup.R), for BOTH sides of every pair. A
# z-scored cytokine against a raw morphology measure would reintroduce the
# confound through the target -- astrocyte cell number is itself sex-dimorphic.
#
# WHY A PERMUTATION NULL, NOT cor.test(). Permuting one variable against the
# other holds BOTH marginals fixed, so the null is conditioned on the exact tie
# structure observed and no tie-correction formula is needed. cor.test() cannot
# do this: it warns "cannot compute exact p-value with ties" and falls back to
# an asymptotic t, which is unreliable at n ~ 20.
#
# WHY THE PERMUTATION IS FREE, NOT SEX-STRATIFIED. Stratifying holds sex fixed
# under the null, which is what raw values require. The sex z-score already
# removes the sex effect from the data itself, so a free permutation is the
# matching null. The stratified p is still computed and written as
# `p_perm_sexstrat` -- a sensitivity check, not the headline. Where the two
# disagree, say so.
#
# Permutations are SHARED across analytes within each (genotype, target) family
# so the joint null keeps the correlation structure among cytokines -- required
# for the Westfall-Young column.
#
# TWO multiplicity columns, same permutations:
#   p_adj = BH on the permutation p-values, within row. Primary.
#   p_wy  = Westfall-Young single-step maxT on STUDENTIZED rho (FWER). Strict.
#           Studentized because the null SD of rho differs with tie burden;
#           max over raw |rho| would favour the coarse analytes.
#   Report both. Coarse analytes pass BH and fail WY -- that gap IS the signal
#   that they are quantisation artefacts.
#
# RESOLUTION COLUMNS ARE MEASURED ON RAW MFI, not on the z-scored values. Ties
# are a property of the assay: an analyte read at four distinct levels is coarse
# no matter what is done to it afterwards, and z-scoring within sex would split
# those ties and flatter the analyte. n_distinct and prop_tied therefore travel
# with the result and describe the instrument, not the transform.
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
  cat(sprintf("[05]   %sfree vs sex-stratified null disagree at p<.05 on %d of %d pairs\n",
              "", d, nrow(out)))
  out
}

invisible(run(AMYLOID_TARGETS, "5XFAD", "cytokine_amyloid"))   # WT has no amyloid
invisible(run(MORPH_TARGETS, GENOTYPES, "cytokine_morphology"))
cat("[05] ->", STATS, "\n")
