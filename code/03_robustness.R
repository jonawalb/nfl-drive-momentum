## 03_robustness.R
## Pre-submission robustness gauntlet for the JQAS revision.
## Seven analyses: 5 CRITICAL pre-submission checks + 2 MAJOR suggested.
##
## CRITICAL:
##   R1. Flexible field-position and score-differential controls
##   R2. Placebo: post-FG vs post-TD (kickoff regime held constant)
##   R3. Logit / probit comparison + average marginal effects
##   R4. Permutation inference for H3 (within-game shuffle)
##   R5. Post-kickoff vs post-punt subsamples
## MAJOR:
##   R6. Heterogeneity by score margin and quarter for H3
##   R7. Rule-era split (2010 / 2011-2017 / 2018-2023 / 2024)

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(modelsummary)
  library(ggplot2)
  library(here)
})

set.seed(20260508)

DATA_DIR <- here::here("data")
TAB_DIR  <- here::here("tables")
FIG_DIR  <- here::here("figures")

d <- readRDS(file.path(DATA_DIR, "drives.rds"))
setDT(d)
d <- d[!is.na(off_success) & !is.na(yardline_100) & !is.na(score_differential)]
d[, def_success := fifelse(is.na(def_success), 0L, def_success)]

d_main <- d[!is.na(prior_own_def_success) & !is.na(prior_own_off_success) &
            !is.na(prior_opp_off_success)]
cat("Robustness sample:", nrow(d_main), "drives\n\n")

results <- list()

## ============================================================================
## R1. Flexible functional form for field position and score differential
## ============================================================================
cat("===== R1: Flexible field-position and score-differential controls =====\n")

d_main[, fp_decile := factor(cut(yardline_100, breaks = seq(0, 100, 10),
                                  include.lowest = TRUE, labels = FALSE))]
d_main[, sd_bin := cut(score_differential,
                       breaks = c(-Inf, -21, -14, -8, -3, 0, 3, 8, 14, 21, Inf),
                       include.lowest = TRUE)]

m_flexFP <- feols(off_success ~ prior_own_def_success + prior_own_off_success +
                    prior_opp_off_success +
                    sd_bin + qtr + half_seconds_remaining |
                    game_id + posteam + fp_decile,
                  data = d_main, cluster = ~game_id + posteam)

cat("  H1 flexible:", round(coef(m_flexFP)["prior_own_def_success"], 4), "\n")
cat("  H2 flexible:", round(coef(m_flexFP)["prior_own_off_success"], 4), "\n")
cat("  H3 flexible:", round(coef(m_flexFP)["prior_opp_off_success"], 4), "\n")

results$R1 <- list(
  beta_h1 = coef(m_flexFP)["prior_own_def_success"],
  beta_h2 = coef(m_flexFP)["prior_own_off_success"],
  beta_h3 = coef(m_flexFP)["prior_opp_off_success"],
  se_h1   = sqrt(diag(vcov(m_flexFP)))["prior_own_def_success"],
  se_h2   = sqrt(diag(vcov(m_flexFP)))["prior_own_off_success"],
  se_h3   = sqrt(diag(vcov(m_flexFP)))["prior_opp_off_success"],
  n_obs   = nobs(m_flexFP)
)

## ============================================================================
## R2. Placebo: post-FG vs post-TD (both lead to kickoff regime)
## ============================================================================
cat("\n===== R2: Placebo (post-FG vs post-TD; kickoff regime held constant) =====\n")

d_main[, prev_was_td := as.integer(prev_drive_result == "Touchdown" &
                                     !is.na(prev_posteam) &
                                     prev_posteam != posteam)]
d_main[, prev_was_fg := as.integer(prev_drive_result == "Field goal" &
                                     !is.na(prev_posteam) &
                                     prev_posteam != posteam)]

placebo_sample <- d_main[prev_was_td == 1 | prev_was_fg == 1]

m_placebo <- feols(off_success ~ prev_was_td +
                     yardline_100 + score_differential + qtr +
                     half_seconds_remaining |
                     game_id + posteam,
                   data = placebo_sample, cluster = ~game_id + posteam)

cat("  Sample size (post-TD or post-FG only):", nrow(placebo_sample), "\n")
cat("  Coefficient on prev_was_td (relative to prev_was_fg):",
    round(coef(m_placebo)["prev_was_td"], 4),
    "(s.e.", round(sqrt(diag(vcov(m_placebo)))["prev_was_td"], 4), ")\n")
cat("  Interpretation: a non-zero coefficient indicates additional H3 signal\n")
cat("  beyond what kickoff-regime field position can absorb.\n")

results$R2 <- list(
  beta_td_vs_fg = coef(m_placebo)["prev_was_td"],
  se_td_vs_fg   = sqrt(diag(vcov(m_placebo)))["prev_was_td"],
  n_obs         = nobs(m_placebo)
)

## ============================================================================
## R3. Logit comparison + average marginal effects + LPM out-of-bound check
## ============================================================================
cat("\n===== R3: Logit/probit comparison and OOB check =====\n")

m_lpm <- feols(off_success ~ prior_own_def_success + prior_own_off_success +
                 prior_opp_off_success +
                 yardline_100 + score_differential + qtr +
                 half_seconds_remaining |
                 game_id + posteam,
               data = d_main, cluster = ~game_id + posteam)

lpm_pred <- predict(m_lpm)
oob_low  <- mean(lpm_pred < 0)
oob_high <- mean(lpm_pred > 1)
cat("  LPM predicted-probability OOB share: <0:",
    sprintf("%.2f%%", 100 * oob_low),
    " >1:", sprintf("%.2f%%", 100 * oob_high), "\n")

m_logit <- feglm(off_success ~ prior_own_def_success + prior_own_off_success +
                   prior_opp_off_success +
                   yardline_100 + score_differential + qtr +
                   half_seconds_remaining |
                   posteam,
                 data = d_main, family = binomial("logit"),
                 cluster = ~game_id + posteam)

ame_logit <- tryCatch({
  beta <- coef(m_logit)
  X    <- model.matrix(m_logit)
  eta  <- as.vector(X %*% beta)
  dens <- exp(eta) / (1 + exp(eta))^2
  list(
    h1 = mean(dens) * beta["prior_own_def_success"],
    h2 = mean(dens) * beta["prior_own_off_success"],
    h3 = mean(dens) * beta["prior_opp_off_success"]
  )
}, error = function(e) {
  cat("  AME computation failed:", e$message, "\n"); NULL
})

if (!is.null(ame_logit)) {
  cat("  Logit AMEs:  H1 =", round(ame_logit$h1, 4),
      "  H2 =", round(ame_logit$h2, 4),
      "  H3 =", round(ame_logit$h3, 4), "\n")
  cat("  LPM coefs:   H1 =", round(coef(m_lpm)["prior_own_def_success"], 4),
      "  H2 =", round(coef(m_lpm)["prior_own_off_success"], 4),
      "  H3 =", round(coef(m_lpm)["prior_opp_off_success"], 4), "\n")
}

results$R3 <- list(
  oob_low  = oob_low,
  oob_high = oob_high,
  ame      = ame_logit,
  lpm_h1   = coef(m_lpm)["prior_own_def_success"],
  lpm_h2   = coef(m_lpm)["prior_own_off_success"],
  lpm_h3   = coef(m_lpm)["prior_opp_off_success"]
)

## ============================================================================
## R4. Permutation inference for H3 (shuffle within team-quarter cells)
## ============================================================================
cat("\n===== R4: Permutation inference for H3 =====\n")

n_perm <- 1000
beta_h3_perm <- numeric(n_perm)
beta_h3_obs  <- coef(m5 <- feols(
  off_success ~ prior_own_def_success + prior_own_off_success +
    prior_opp_off_success +
    yardline_100 + score_differential + qtr + half_seconds_remaining |
    game_id + posteam,
  data = d_main, cluster = ~game_id + posteam))["prior_opp_off_success"]

cat("  Observed H3:", round(beta_h3_obs, 4), "\n")
cat("  Running", n_perm, "permutations within team-quarter cells...\n")

d_perm_base <- copy(d_main)
for (b in seq_len(n_perm)) {
  d_perm_base[, prior_opp_off_perm := sample(prior_opp_off_success), by = .(posteam, qtr)]
  m_perm <- feols(off_success ~ prior_own_def_success + prior_own_off_success +
                    prior_opp_off_perm +
                    yardline_100 + score_differential + qtr +
                    half_seconds_remaining |
                    game_id + posteam,
                  data = d_perm_base, se = "iid")
  beta_h3_perm[b] <- coef(m_perm)["prior_opp_off_perm"]
  if (b %% 100 == 0) cat("    perm", b, "/", n_perm,
                          "  current beta:", round(beta_h3_perm[b], 4), "\n")
}

p_value_perm <- mean(abs(beta_h3_perm) >= abs(beta_h3_obs))
cat("  Permutation p-value (two-sided): ", p_value_perm, "\n")
cat("  Range of permuted H3 estimates: [",
    round(min(beta_h3_perm), 4), ",", round(max(beta_h3_perm), 4), "]\n")

results$R4 <- list(
  beta_h3_obs   = beta_h3_obs,
  perm_min      = min(beta_h3_perm),
  perm_max      = max(beta_h3_perm),
  perm_mean     = mean(beta_h3_perm),
  perm_sd       = sd(beta_h3_perm),
  p_value_perm  = p_value_perm,
  n_perm        = n_perm
)

## ============================================================================
## R5. Post-kickoff vs post-punt subsamples
## ============================================================================
cat("\n===== R5: Post-kickoff vs post-punt subsamples =====\n")

d_main[, prev_was_score := as.integer(prev_drive_result %in% c("Touchdown","Field goal") &
                                        !is.na(prev_posteam) &
                                        prev_posteam != posteam)]
d_main[, prev_was_punt  := as.integer(prev_drive_result == "Punt" &
                                        !is.na(prev_posteam) &
                                        prev_posteam != posteam)]
d_main[, prev_was_to    := as.integer(prev_drive_result %in%
                                        c("Turnover","Turnover on downs","Safety") &
                                        !is.na(prev_posteam) &
                                        prev_posteam != posteam)]

post_kickoff <- d_main[prev_was_score == 1]
post_punt    <- d_main[prev_was_punt == 1]
post_to      <- d_main[prev_was_to == 1]

m_kickoff <- feols(off_success ~ prior_own_def_success + prior_own_off_success +
                     yardline_100 + score_differential + qtr +
                     half_seconds_remaining |
                     game_id + posteam,
                   data = post_kickoff, cluster = ~game_id + posteam)

m_punt <- feols(off_success ~ prior_own_def_success + prior_own_off_success +
                  yardline_100 + score_differential + qtr +
                  half_seconds_remaining |
                  game_id + posteam,
                data = post_punt, cluster = ~game_id + posteam)

cat("  Post-kickoff (opp scored) sample:", nrow(post_kickoff),
    "  P(score):", round(mean(post_kickoff$off_success), 3), "\n")
cat("  Post-punt sample:                ", nrow(post_punt),
    "  P(score):", round(mean(post_punt$off_success), 3), "\n")
cat("  Post-turnover sample:            ", nrow(post_to),
    "  P(score):", round(mean(post_to$off_success), 3), "\n")

results$R5 <- list(
  n_post_kickoff = nrow(post_kickoff),
  n_post_punt    = nrow(post_punt),
  n_post_to      = nrow(post_to),
  pscore_post_kickoff = mean(post_kickoff$off_success),
  pscore_post_punt    = mean(post_punt$off_success),
  pscore_post_to      = mean(post_to$off_success)
)

## ============================================================================
## R6. Heterogeneity: H3 by score margin and quarter
## ============================================================================
cat("\n===== R6: Heterogeneity by score margin and quarter =====\n")

d_main[, abs_sd := abs(score_differential)]
d_main[, close_game := as.integer(abs_sd <= 7)]
d_main[, late_game  := as.integer(qtr >= 4)]

m_het <- feols(off_success ~ prior_own_def_success + prior_own_off_success +
                 prior_opp_off_success * close_game +
                 prior_opp_off_success * late_game +
                 yardline_100 + score_differential + qtr +
                 half_seconds_remaining |
                 game_id + posteam,
               data = d_main, cluster = ~game_id + posteam)

cat("  H3 main effect:                       ",
    round(coef(m_het)["prior_opp_off_success"], 4), "\n")
cat("  H3 x close_game (within 7 pts):       ",
    round(coef(m_het)["prior_opp_off_success:close_game"], 4), "\n")
cat("  H3 x late_game (Q4+):                 ",
    round(coef(m_het)["prior_opp_off_success:late_game"], 4), "\n")
cat("  Psychological story predicts: H3 stronger in close games (negative interaction)\n")
cat("  Mechanical story predicts:    H3 weaker in close games (positive interaction)\n")

results$R6 <- list(
  beta_h3_main      = coef(m_het)["prior_opp_off_success"],
  beta_h3_x_close   = coef(m_het)["prior_opp_off_success:close_game"],
  beta_h3_x_late    = coef(m_het)["prior_opp_off_success:late_game"],
  se_h3_x_close     = sqrt(diag(vcov(m_het)))["prior_opp_off_success:close_game"],
  se_h3_x_late      = sqrt(diag(vcov(m_het)))["prior_opp_off_success:late_game"]
)

## ============================================================================
## R7. Rule-era split (2010 / 2011-2017 / 2018-2023 / 2024)
## ============================================================================
cat("\n===== R7: Rule-era split =====\n")

d_main[, era := fcase(
  season == 2010,                 "2010 (pre-touchback move)",
  season >= 2011 & season <= 2017, "2011-2017 (touchback at 20)",
  season >= 2018 & season <= 2023, "2018-2023 (touchback at 25)",
  season == 2024,                 "2024 (dynamic kickoff)",
  default = NA_character_
)]

eras <- unique(d_main$era)
era_results <- list()
for (e in eras) {
  sub <- d_main[era == e]
  if (nrow(sub) < 100) next
  m_e <- tryCatch(
    feols(off_success ~ prior_own_def_success + prior_own_off_success +
            prior_opp_off_success +
            yardline_100 + score_differential + qtr +
            half_seconds_remaining |
            game_id + posteam,
          data = sub, cluster = ~game_id + posteam),
    error = function(err) { cat("  era", e, "fit failed:", err$message, "\n"); NULL }
  )
  if (is.null(m_e)) next
  era_results[[e]] <- list(
    beta_h1 = coef(m_e)["prior_own_def_success"],
    beta_h2 = coef(m_e)["prior_own_off_success"],
    beta_h3 = coef(m_e)["prior_opp_off_success"],
    se_h3   = sqrt(diag(vcov(m_e)))["prior_opp_off_success"],
    n_obs   = nobs(m_e)
  )
  cat(sprintf("  %s  (N=%d)  H1=%6.3f  H2=%6.3f  H3=%6.3f\n",
              e, nobs(m_e),
              coef(m_e)["prior_own_def_success"],
              coef(m_e)["prior_own_off_success"],
              coef(m_e)["prior_opp_off_success"]))
}

results$R7 <- era_results

## ============================================================================
## Save robustness results
## ============================================================================
saveRDS(results, file.path(DATA_DIR, "robustness.rds"))

## Build a single robustness summary table
build_robustness_table <- function(results) {
  rows <- c()
  rows <- c(rows, sprintf("R1 Flexible FP+SD controls   N=%d  H1=%.3f  H2=%.3f  H3=%.3f",
                          results$R1$n_obs, results$R1$beta_h1, results$R1$beta_h2, results$R1$beta_h3))
  rows <- c(rows, sprintf("R2 Placebo (post-TD vs FG)   N=%d  TD vs FG coef=%.3f (se=%.3f)",
                          results$R2$n_obs, results$R2$beta_td_vs_fg, results$R2$se_td_vs_fg))
  rows <- c(rows, sprintf("R3 LPM OOB share  <0: %.2f%%  >1: %.2f%%",
                          100*results$R3$oob_low, 100*results$R3$oob_high))
  if (!is.null(results$R3$ame)) {
    rows <- c(rows, sprintf("R3 Logit AMEs     H1=%.3f  H2=%.3f  H3=%.3f",
                            results$R3$ame$h1, results$R3$ame$h2, results$R3$ame$h3))
  }
  rows <- c(rows, sprintf("R4 Permutation H3  obs=%.3f  null range=[%.3f, %.3f]  p=%.4f",
                          results$R4$beta_h3_obs, results$R4$perm_min, results$R4$perm_max,
                          results$R4$p_value_perm))
  rows <- c(rows, sprintf("R5 Post-kickoff (opp scored) N=%d  P(score)=%.3f",
                          results$R5$n_post_kickoff, results$R5$pscore_post_kickoff))
  rows <- c(rows, sprintf("R5 Post-punt              N=%d  P(score)=%.3f",
                          results$R5$n_post_punt, results$R5$pscore_post_punt))
  rows <- c(rows, sprintf("R5 Post-turnover          N=%d  P(score)=%.3f",
                          results$R5$n_post_to, results$R5$pscore_post_to))
  rows <- c(rows, sprintf("R6 H3 x close_game        coef=%.3f (se=%.3f)",
                          results$R6$beta_h3_x_close, results$R6$se_h3_x_close))
  rows <- c(rows, sprintf("R6 H3 x late_game         coef=%.3f (se=%.3f)",
                          results$R6$beta_h3_x_late, results$R6$se_h3_x_late))
  for (e in names(results$R7)) {
    r <- results$R7[[e]]
    rows <- c(rows, sprintf("R7 era %s  N=%d  H3=%.3f (se=%.3f)",
                            e, r$n_obs, r$beta_h3, r$se_h3))
  }
  rows
}

cat("\n\n========== ROBUSTNESS SUMMARY ==========\n")
cat(paste(build_robustness_table(results), collapse = "\n"), "\n")

cat("\nSaved: data/robustness.rds\n")
cat("Done.\n")
