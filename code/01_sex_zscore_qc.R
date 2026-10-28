# ============================================================================
# 01_sex_zscore_qc.R  —  does the sex z-score hold up? Run it; read the output.
#
# Z-scoring within sex and pooling buys statistical power, but it is only valid
# if removing each sex's mean does not also remove part of the flicker effect.
# That holds when sex is BALANCED across flicker arms within a genotype. If one
# arm were mostly female, the female mean would partly BE that arm's mean, and
# subtracting it would erase the very signal the study is looking for.
#
# This stage therefore checks three things and writes them to stats/:
#
#   1. BALANCE. Sex composition per flicker arm, and the largest deviation from
#      even. Hard stop if any arm is worse than 80/20 -- past that the z-score
#      is removing flicker, not sex.
#   2. WHAT IT REMOVED. Per analyte, the sex difference before and after. After
#      z-scoring the sex effect must be ~0 by construction; this confirms the
#      arithmetic did what the Methods sentence claims.
#   3. WHAT IT KEPT. Per analyte, the flicker omnibus p before and after. These
#      should track each other. Analytes where they diverge sharply are where
#      the sex correction mattered most -- worth knowing, not necessarily wrong.
#
# This stage writes tables and prints a verdict. It draws nothing.
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
source(file.path(.here(), "00_setup.R"))

# ---------------------------------------------------------------------------
# 1. Balance
# ---------------------------------------------------------------------------
bal <- do.call(rbind, lapply(GENOTYPES, function(g) {
  d <- ANIMALS_RAW[ANIMALS_RAW$Genotype == g, ]
  do.call(rbind, lapply(FLICKER_LEVELS, function(f) {
    i <- d$Flicker == f
    nF <- sum(i & d$Sex == "F"); nM <- sum(i & d$Sex == "M")
    data.frame(genotype = g, flicker = f, n_F = nF, n_M = nM, n = nF + nM,
               frac_F = nF / (nF + nM))
  }))
}))
bal$deviation <- abs(bal$frac_F - 0.5)
write.csv(bal, file.path(STATS, "qc_sex_balance.csv"), row.names = FALSE)

cat("=== 1. SEX BALANCE ACROSS FLICKER ARMS ===\n")
print(bal, row.names = FALSE)
worst <- max(bal$frac_F[bal$n > 0], 1 - bal$frac_F[bal$n > 0])
cat(sprintf("\n  worst arm is %.0f%% one sex (hard limit 80%%)\n", 100 * worst))
if (worst > 0.80)
  stop("SEX BALANCE FAILED: an arm exceeds 80/20. The within-sex z-score would ",
       "remove flicker signal, not just sex. Re-decide the sex handling.")
cat("  -> balance OK; pooling after a within-sex z-score is safe here.\n\n")

# ---------------------------------------------------------------------------
# 2. What the z-score removed (the sex effect) and 3. what it kept (flicker)
# ---------------------------------------------------------------------------
sex_p <- function(d, v) {
  x <- d[[v]]; s <- d$Sex
  if (sum(!is.na(x)) < 4 || length(unique(s[!is.na(x)])) < 2) return(NA_real_)
  tryCatch(stats::wilcox.test(x ~ s, exact = FALSE)$p.value, error = function(e) NA_real_)
}
flick_p <- function(d, v) {
  x <- d[[v]]; f <- droplevels(d$Flicker)
  if (sum(!is.na(x)) < 6 || nlevels(f) < 2) return(NA_real_)
  tryCatch(stats::kruskal.test(x ~ f)$p.value, error = function(e) NA_real_)
}

eff <- do.call(rbind, lapply(GENOTYPES, function(g) {
  dr <- ANIMALS_RAW[ANIMALS_RAW$Genotype == g, ]
  dz <- ANIMALS    [ANIMALS$Genotype     == g, ]
  do.call(rbind, lapply(PANEL, function(v) data.frame(
    genotype = g, analyte = v,
    sex_p_raw = sex_p(dr, v), sex_p_z = sex_p(dz, v),
    flicker_p_raw = flick_p(dr, v), flicker_p_z = flick_p(dz, v))))
}))
eff$flicker_p_shift <- eff$flicker_p_z - eff$flicker_p_raw
write.csv(eff, file.path(STATS, "qc_sex_zscore_effect.csv"), row.names = FALSE)

cat("=== 2. SEX EFFECT, BEFORE vs AFTER THE Z-SCORE ===\n")
cat(sprintf("  analytes with a sex effect at p < 0.05, before : %d of %d\n",
            sum(eff$sex_p_raw < 0.05, na.rm = TRUE), nrow(eff)))
cat(sprintf("                                          after  : %d of %d\n",
            sum(eff$sex_p_z   < 0.05, na.rm = TRUE), nrow(eff)))
top <- eff[order(eff$sex_p_raw), ][1:6, ]
cat("  strongest sex effects that were removed:\n")
for (i in seq_len(nrow(top))) cat(sprintf("    %-6s %-12s p %8.2g -> %.2f\n",
    top$genotype[i], top$analyte[i], top$sex_p_raw[i], top$sex_p_z[i]))

cat("\n=== 3. FLICKER EFFECT, BEFORE vs AFTER (should track) ===\n")
ok <- !is.na(eff$flicker_p_raw) & !is.na(eff$flicker_p_z)
cat(sprintf("  Spearman correlation of the two p-value sets: %.3f\n",
            suppressWarnings(cor(eff$flicker_p_raw[ok], eff$flicker_p_z[ok], method = "spearman"))))
cat(sprintf("  reach p < 0.05 before: %d   after: %d\n",
            sum(eff$flicker_p_raw < 0.05, na.rm = TRUE),
            sum(eff$flicker_p_z   < 0.05, na.rm = TRUE)))
mv <- eff[ok, ][order(-abs(eff$flicker_p_shift[ok])), ][1:5, ]
cat("  largest movers (where correcting for sex mattered most):\n")
for (i in seq_len(nrow(mv))) cat(sprintf("    %-6s %-12s flicker p %.3f -> %.3f\n",
    mv$genotype[i], mv$analyte[i], mv$flicker_p_raw[i], mv$flicker_p_z[i]))

cat("\n[01] wrote qc_sex_balance.csv + qc_sex_zscore_effect.csv ->", STATS, "\n")
