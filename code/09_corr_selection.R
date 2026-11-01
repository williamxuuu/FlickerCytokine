# ============================================================================
# 09_corr_selection.R  —  cytokine ~ morphology, by HOW the cytokines are chosen.
#
# Stage 05 already correlates every one of the 28 analytes with every
# morphology target, so choosing a different analyte set does NOT change any
# rho or p. What it can change is the multiplicity family. This stage puts two
# selection rules side by side, on the four morphology targets the figures use:
#
#   key         the PLS-DA key cytokines (stage 02): top 3 |loading| on LV1
#               and on LV2. What Fig 5C/D and Supp 7 show.
#   responsive  the N_RESP analytes with the smallest flicker omnibus p in that
#               genotype (stage 06). In WT the top five reach p < 0.05; in
#               5xFAD none do, so the set there is "least null", not
#               "significant". `omnibus_p` travels with every row.
#
# For each row it reports:
#   q_panel   BH across all 28 analytes for that (genotype, target) -- the
#             primary q, identical to stage 05's p_adj, whatever the selection.
#   q_set     BH within the selected set (6 analytes) for that target. Smaller
#             family, so smaller q. Only defensible if the set is fixed BEFORE
#             the correlations are looked at -- both rules here select on a
#             separate analysis (PLS-DA / group test), never on the correlation.
#   p_wy      Westfall-Young across the 28 analytes, unchanged from stage 05.
#
# Read q_set as a sensitivity view, not as a replacement for q_panel: picking
# whichever family gives the smaller q after seeing both is not a correction.
# ============================================================================

.here <- function() {
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(normalizePath(f)) else getwd()
}
source(file.path(.here(), "00_setup.R"))

N_RESP <- 6        # same size as the key set (3 + 3), so the families match
FIG_TARGETS <- c("AST mean_Distance to Nearest Neighbor", "MG mean_Distance to Nearest Neighbor",
                 "AST mean_Number of Cells",              "MG mean_Number of Cells")

CM <- read.csv(file.path(STATS, "corr_cytokine_morphology.csv"), check.names = FALSE)
OM <- read.csv(file.path(STATS, "group_omnibus_cytokines.csv"),  check.names = FALSE)

responsive <- function(g) {
  o <- OM[OM$subgroup == g & OM$variable %in% PANEL, ]
  o$variable[order(o$p)][seq_len(N_RESP)]
}
key <- function(g) readRDS(file.path(STATS, sprintf("plsda_%s.rds", g)))$key

out <- do.call(rbind, lapply(GENOTYPES, function(g) {
  sets <- list(key = key(g), responsive = responsive(g))
  do.call(rbind, lapply(names(sets), function(s) do.call(rbind, lapply(FIG_TARGETS, function(tg) {
    r <- CM[CM$subgroup == g & CM$target == tg & CM$analyte %in% sets[[s]], ]
    data.frame(selection = s, subgroup = g, target = tg, analyte = r$analyte,
               omnibus_p = OM$p[match(paste(g, r$analyte), paste(OM$subgroup, OM$variable))],
               rho = r$rho, n = r$n, p_perm = r$p_perm,
               q_panel = r$p_adj, q_set = p.adjust(r$p_perm, "BH"), p_wy = r$p_wy,
               n_distinct = r$n_distinct, row.names = NULL)
  }))))
}))
out <- out[order(out$subgroup, out$selection, out$target, out$p_perm), ]
write.csv(out, file.path(STATS, "corr_morphology_by_selection.csv"), row.names = FALSE)

for (g in GENOTYPES) for (s in c("key", "responsive")) {
  o <- out[out$subgroup == g & out$selection == s, ]
  cat(sprintf("[09] %-5s %-10s %s\n", g, s, paste(unique(o$analyte), collapse = ", ")))
  cat(sprintf("[09]       %2d pairs | p < .05: %d | q_panel < .05: %d | q_set < .05: %d | WY < .05: %d\n",
              nrow(o), sum(o$p_perm < .05), sum(o$q_panel < .05),
              sum(o$q_set < .05), sum(o$p_wy < .05)))
}

# For contrast: what clears q_panel in each genotype, whatever the selection.
best <- CM[CM$target %in% FIG_TARGETS & CM$p_adj < 0.05, ]
cat("[09] pairs with q_panel < .05 on the figure targets (any analyte):\n")
for (i in seq_len(nrow(best)))
  cat(sprintf("       %-5s %-8s ~ %-38s rho = %5.2f  q = %.3f  WY = %.3f\n", best$subgroup[i],
              best$analyte[i], best$target[i], best$rho[i], best$p_adj[i], best$p_wy[i]))
cat("[09] -> stats/corr_morphology_by_selection.csv\n")
