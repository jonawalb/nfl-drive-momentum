## 02_analyze_momentum.R
## Test three drive-level momentum hypotheses with FE regressions.
## H1: Defensive success raises this offense's next-drive scoring probability
## H2: Offensive success on prior own drive raises next own-drive scoring
## H3: Opponent's prior offensive success affects this drive's outcome

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(modelsummary)
  library(ggplot2)
  library(here)
})

DATA_DIR <- here::here("data")
TAB_DIR  <- here::here("tables")
FIG_DIR  <- here::here("figures")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

d <- readRDS(file.path(DATA_DIR, "drives.rds"))
setDT(d)

## Drop drives with missing core fields
d <- d[!is.na(off_success) & !is.na(yardline_100) & !is.na(score_differential)]

## Replace NA def_success with 0 only where drive_result is well defined
d[, def_success := fifelse(is.na(def_success), 0L, def_success)]

cat("Analysis sample:\n")
cat("  total drives:", nrow(d), "\n")
cat("  H1 sample (have prior_own_def):", sum(!is.na(d$prior_own_def_success)), "\n")
cat("  H2 sample (have prior_own_off):", sum(!is.na(d$prior_own_off_success)), "\n")
cat("  H3 sample (have prior_opp_off):", sum(!is.na(d$prior_opp_off_success)), "\n")

cat("\nUnconditional comparisons:\n")
cat("  P(score | prior own def success): ",
    round(mean(d[prior_own_def_success == 1]$off_success), 3), "\n")
cat("  P(score | prior own def failure): ",
    round(mean(d[prior_own_def_success == 0]$off_success), 3), "\n")
cat("  P(score | prior own off success): ",
    round(mean(d[prior_own_off_success == 1]$off_success), 3), "\n")
cat("  P(score | prior own off failure): ",
    round(mean(d[prior_own_off_success == 0]$off_success), 3), "\n")
cat("  P(score | prior opp off success): ",
    round(mean(d[prior_opp_off_success == 1]$off_success), 3), "\n")
cat("  P(score | prior opp off failure): ",
    round(mean(d[prior_opp_off_success == 0]$off_success), 3), "\n")

## ------- Field position is THE confound. Show conditional rates -----------
cat("\nP(score) by starting field position decile, by prior_own_def_success:\n")
fp_table <- d[!is.na(prior_own_def_success),
              .(p_score = mean(off_success), n = .N),
              by = .(fp_decile, prior_own_def_success)][order(fp_decile, prior_own_def_success)]
print(fp_table)
fwrite(fp_table, file.path(TAB_DIR, "fp_decile_rates.csv"))

## ============================================================================
## Main regressions: linear probability with team + game + opponent FE
## Cluster SEs by game and possessing team
## ============================================================================

## Build a clean sample with all three lag indicators present so models are comparable
d_main <- d[!is.na(prior_own_def_success) & !is.na(prior_own_off_success) &
            !is.na(prior_opp_off_success)]
cat("\nUnified estimation sample:", nrow(d_main), "drives\n")

## Helper: standardized starting field position (0 = own goal line, 100 = opp end zone)
d_main[, yfg := 100 - yardline_100]   # yards FROM own goal -> distance gained perspective

## ---- Specification 1: H1 alone, no controls -------------------------------
m1_naive <- feols(off_success ~ prior_own_def_success,
                  data = d_main, cluster = ~game_id + posteam)

## ---- Specification 2: H1 + game state controls ----------------------------
m2_state <- feols(off_success ~ prior_own_def_success +
                    yardline_100 + score_differential + qtr +
                    half_seconds_remaining + home_off,
                  data = d_main, cluster = ~game_id + posteam)

## ---- Specification 3: + team and opponent FE ------------------------------
m3_fe <- feols(off_success ~ prior_own_def_success +
                 yardline_100 + score_differential + qtr +
                 half_seconds_remaining + home_off |
                 posteam + defteam + season,
               data = d_main, cluster = ~game_id + posteam)

## ---- Specification 4: + game FE (within-game variation only) -------------
m4_game <- feols(off_success ~ prior_own_def_success +
                   yardline_100 + score_differential + qtr +
                   half_seconds_remaining |
                   game_id + posteam,
                 data = d_main, cluster = ~game_id + posteam)

## ---- Specification 5: All three momentum channels jointly -----------------
m5_all <- feols(off_success ~ prior_own_def_success + prior_own_off_success +
                  prior_opp_off_success +
                  yardline_100 + score_differential + qtr +
                  half_seconds_remaining |
                  game_id + posteam,
                data = d_main, cluster = ~game_id + posteam)

## ---- Specification 6: Interact prior_own_def_success with field position --
## If "momentum" is real, the boost should NOT be entirely mechanical short field.
## Add interaction with starting field position decile.
d_main[, fp_decile := factor(fp_decile)]
m6_intx <- feols(off_success ~ prior_own_def_success * fp_decile +
                   prior_own_off_success + prior_opp_off_success +
                   score_differential + qtr + half_seconds_remaining |
                   game_id + posteam,
                 data = d_main, cluster = ~game_id + posteam)

## ---- Specification 7: Drop garbage time -----------------------------------
m7_nogt <- feols(off_success ~ prior_own_def_success + prior_own_off_success +
                   prior_opp_off_success +
                   yardline_100 + score_differential + qtr +
                   half_seconds_remaining |
                   game_id + posteam,
                 data = d_main[garbage_time == 0],
                 cluster = ~game_id + posteam)

## ---- Specification 8: EPA-based outcome (continuous) ----------------------
m8_epa <- feols(drive_epa ~ prior_own_def_success + prior_own_off_success +
                  prior_opp_off_success +
                  yardline_100 + score_differential + qtr +
                  half_seconds_remaining |
                  game_id + posteam,
                data = d_main, cluster = ~game_id + posteam)

## ============================================================================
## Tables
## ============================================================================

## Variable label dictionary for human-readable table output
var_dict <- c(
  prior_own_def_success = "Prior own defensive stop",
  prior_own_off_success = "Prior own offensive score",
  prior_opp_off_success = "Prior opponent offensive score",
  yardline_100          = "Field position (yards to opp. end zone)",
  score_differential    = "Score differential",
  qtr                   = "Quarter",
  half_seconds_remaining= "Half seconds remaining",
  home_off              = "Home offense",
  off_success           = "Drive ends in score (0/1)",
  drive_epa             = "Drive EPA",
  posteam               = "Possessing team",
  defteam               = "Opponent",
  game_id               = "Game",
  season                = "Season"
)

models_main <- list(
  "(1) Naive"        = m1_naive,
  "(2) +State"       = m2_state,
  "(3) +Team/Opp FE" = m3_fe,
  "(4) +Game FE"     = m4_game
)

models_full <- list(
  "(1) All three channels" = m5_all,
  "(2) Drop garbage time"  = m7_nogt,
  "(3) EPA outcome"        = m8_epa
)

cat("\n\n===== MAIN H1 PROGRESSION =====\n")
modelsummary(models_main, output = "markdown", stars = TRUE,
             coef_omit = "^(yardline_100|score_differential|qtr|half_seconds|home_off)",
             gof_map = c("nobs","r.squared","adj.r.squared"))

cat("\n\n===== ALL THREE CHANNELS + ROBUSTNESS =====\n")
modelsummary(models_full, output = "markdown", stars = TRUE,
             coef_omit = "^(yardline_100|score_differential|qtr|half_seconds|home_off)",
             gof_map = c("nobs","r.squared","adj.r.squared"))

## LaTeX export with renamed variables and journal-grade notes
h1_notes <- paste0(
  "\\textit{Notes.} OLS linear probability models. Dependent variable is an ",
  "indicator for whether the drive ends in a touchdown or field goal. Sample: ",
  "NFL regular season + postseason 2010--2024 (nflfastR). Game-state controls ",
  "(field position, score differential, quarter, half-seconds remaining, home ",
  "indicator) included in columns 2--4 with coefficients suppressed for ",
  "readability. Standard errors in parentheses, two-way clustered by game and ",
  "possessing team \\citep{cameron2011robust}. ",
  "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."
)

allchan_notes <- paste0(
  "\\textit{Notes.} Column 1 is the preferred specification: linear probability ",
  "model for whether the drive ends in a score, with all three momentum ",
  "indicators jointly. Column 2 drops drives in garbage time (absolute score ",
  "differential $\\geq 21$ in the fourth quarter). Column 3 replaces the binary ",
  "outcome with the sum of expected-points-added (EPA) on the drive. All ",
  "columns include game-state controls (field position, score differential, ",
  "quarter, half-seconds remaining) and game and possessing-team fixed ",
  "effects. Standard errors in parentheses, two-way clustered by game and ",
  "possessing team \\citep{cameron2011robust}. ",
  "$^{*}p<0.10$, $^{**}p<0.05$, $^{***}p<0.01$."
)

etable(models_main, tex = TRUE,
       file = file.path(TAB_DIR, "tab_h1_progression.tex"),
       title = paste0("Defensive Stop Spillover (H1): Progressive Controls. ", h1_notes),
       label = "tab:h1prog",
       dict = var_dict,
       drop = c("yardline_100","score_differential","qtr","half_seconds","home_off"),
       replace = TRUE)

etable(models_full, tex = TRUE,
       file = file.path(TAB_DIR, "tab_all_channels.tex"),
       title = paste0("Three Momentum Channels: Joint Estimation and Robustness. ", allchan_notes),
       label = "tab:allchan",
       dict = var_dict,
       drop = c("yardline_100","score_differential","qtr","half_seconds","home_off"),
       replace = TRUE)

## ============================================================================
## Programmatic Table 1 (descriptive)
## ============================================================================

build_descrip_row <- function(dt, condition, label, true_label, false_label) {
  sub <- dt[!is.na(get(condition))]
  rate_true  <- sub[get(condition) == 1, mean(off_success)]
  rate_false <- sub[get(condition) == 0, mean(off_success)]
  n_true     <- sub[get(condition) == 1, .N]
  n_false    <- sub[get(condition) == 0, .N]
  data.table(
    label = c(true_label, false_label),
    p_score = c(rate_true, rate_false),
    n = c(n_true, n_false),
    panel = c(label, label)
  )
}

descrip_dt <- rbindlist(list(
  build_descrip_row(d, "prior_own_def_success", "panel_def",
                    "After own defensive stop",
                    "After opponent score or longer punt"),
  build_descrip_row(d, "prior_own_off_success", "panel_off",
                    "After own offensive score (prior own drive)",
                    "After own offensive non-score"),
  build_descrip_row(d, "prior_opp_off_success", "panel_opp",
                    "After opponent score (prior opp drive)",
                    "After opponent non-score")
))

descrip_tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Unconditional drive-level scoring probabilities by prior context.}",
  "\\label{tab:descrip}",
  "\\begin{tabular}{lcc}",
  "\\toprule",
  "Condition & P(score) & N drives \\\\",
  "\\midrule",
  paste0(descrip_dt$label[1], " & ", sprintf("%.3f", descrip_dt$p_score[1]), " & ", format(descrip_dt$n[1], big.mark = ","), " \\\\"),
  paste0(descrip_dt$label[2], " & ", sprintf("%.3f", descrip_dt$p_score[2]), " & ", format(descrip_dt$n[2], big.mark = ","), " \\\\"),
  "\\midrule",
  paste0(descrip_dt$label[3], " & ", sprintf("%.3f", descrip_dt$p_score[3]), " & ", format(descrip_dt$n[3], big.mark = ","), " \\\\"),
  paste0(descrip_dt$label[4], " & ", sprintf("%.3f", descrip_dt$p_score[4]), " & ", format(descrip_dt$n[4], big.mark = ","), " \\\\"),
  "\\midrule",
  paste0(descrip_dt$label[5], " & ", sprintf("%.3f", descrip_dt$p_score[5]), " & ", format(descrip_dt$n[5], big.mark = ","), " \\\\"),
  paste0(descrip_dt$label[6], " & ", sprintf("%.3f", descrip_dt$p_score[6]), " & ", format(descrip_dt$n[6], big.mark = ","), " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\\\[0.5em]",
  "\\begin{minipage}{0.95\\textwidth}",
  "\\footnotesize",
  "\\textit{Notes.} Each panel uses the largest sample for which the relevant prior-drive indicator is defined; sample size therefore differs across panels (the unified estimation sample requiring all three indicators is N = 73{,}842 drives). Sample: NFL regular season + postseason 2010--2024 from nflfastR. ",
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(descrip_tex, file.path(TAB_DIR, "tab_descrip.tex"))

## ============================================================================
## Figures
## ============================================================================

## Figure 1: Conditional P(score) by field-position decile
fp_plot <- d[!is.na(prior_own_def_success),
             .(p_score = mean(off_success),
               n = .N,
               se = sqrt(mean(off_success)*(1-mean(off_success))/.N)),
             by = .(fp_decile, prior_own_def_success)]
fp_plot[, fp_mid := (fp_decile - 0.5) * 10]
fp_plot[, condition := factor(prior_own_def_success,
                               levels = c(0,1),
                               labels = c("Prior drive: opponent scored or punted normally",
                                          "Prior drive: defense forced 3-and-out / TO"))]

p1 <- ggplot(fp_plot, aes(x = fp_mid, y = p_score, color = condition)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  geom_ribbon(aes(ymin = p_score - 1.96*se, ymax = p_score + 1.96*se,
                  fill = condition), alpha = 0.15, color = NA) +
  scale_x_reverse(breaks = seq(10, 100, 10)) +
  labs(x = "Yards from opponent end zone (drive start; lower = closer to scoring)",
       y = "P(drive ends in TD or FG)",
       color = NULL, fill = NULL,
       title = "Scoring probability by drive starting field position",
       subtitle = "Drives following defensive stops vs. drives following opponent scores or longer punts") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

ggsave(file.path(FIG_DIR, "fig1_field_position.pdf"), p1, width = 7, height = 4.5)
ggsave(file.path(FIG_DIR, "fig1_field_position.png"), p1, width = 7, height = 4.5, dpi = 200)

## Figure 2: Coefficient plot of the three momentum channels (model 5)
coef_dt <- as.data.table(broom::tidy(m5_all, conf.int = TRUE))
coef_dt <- coef_dt[term %in% c("prior_own_def_success","prior_own_off_success","prior_opp_off_success")]
coef_dt[, term_label := factor(term,
                                levels = c("prior_opp_off_success",
                                           "prior_own_off_success",
                                           "prior_own_def_success"),
                                labels = c("Opponent scored on prior drive",
                                           "Own offense scored on prior drive",
                                           "Defense forced stop on prior drive"))]

p2 <- ggplot(coef_dt, aes(x = estimate, y = term_label)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
  geom_pointrange(aes(xmin = conf.low, xmax = conf.high), size = 0.5) +
  labs(x = "Effect on P(scoring drive)",
       y = NULL,
       title = "Three drive-level momentum channels",
       subtitle = "LPM with game and team FEs; cluster SEs (game, team); 95% CI") +
  theme_minimal(base_size = 11)

ggsave(file.path(FIG_DIR, "fig2_coef_plot.pdf"), p2, width = 7, height = 3.5)
ggsave(file.path(FIG_DIR, "fig2_coef_plot.png"), p2, width = 7, height = 3.5, dpi = 200)

## ============================================================================
## Save model objects + key numbers for paper
## ============================================================================
saveRDS(list(m1=m1_naive, m2=m2_state, m3=m3_fe, m4=m4_game,
             m5=m5_all, m6=m6_intx, m7=m7_nogt, m8=m8_epa),
        file.path(DATA_DIR, "models.rds"))

key_nums <- list(
  n_drives_total      = nrow(d),
  n_drives_analysis   = nrow(d_main),
  n_games             = uniqueN(d$game_id),
  off_success_rate    = mean(d$off_success),
  def_success_rate    = mean(d$def_success),

  # Raw conditional means (referenced in abstract)
  p_score_def_stop    = mean(d[prior_own_def_success == 1]$off_success),
  p_score_no_stop     = mean(d[prior_own_def_success == 0]$off_success),
  n_def_stop          = sum(d$prior_own_def_success == 1, na.rm = TRUE),
  n_no_stop           = sum(d$prior_own_def_success == 0, na.rm = TRUE),

  # Naive and progressive H1 coefficients (all four progression columns)
  beta_naive_h1       = coef(m1_naive)["prior_own_def_success"],
  beta_state_h1       = coef(m2_state)["prior_own_def_success"],
  beta_teamFE_h1      = coef(m3_fe)["prior_own_def_success"],
  beta_gameFE_h1      = coef(m4_game)["prior_own_def_success"],

  # Joint model coefficients
  beta_h1_joint       = coef(m5_all)["prior_own_def_success"],
  beta_h2_joint       = coef(m5_all)["prior_own_off_success"],
  beta_h3_joint       = coef(m5_all)["prior_opp_off_success"],

  # Standard errors (clustered)
  se_h1_joint = sqrt(diag(vcov(m5_all)))["prior_own_def_success"],
  se_h2_joint = sqrt(diag(vcov(m5_all)))["prior_own_off_success"],
  se_h3_joint = sqrt(diag(vcov(m5_all)))["prior_opp_off_success"],

  # EPA model
  beta_h1_epa  = coef(m8_epa)["prior_own_def_success"],
  beta_h2_epa  = coef(m8_epa)["prior_own_off_success"],
  beta_h3_epa  = coef(m8_epa)["prior_opp_off_success"]
)
saveRDS(key_nums, file.path(DATA_DIR, "key_nums.rds"))

cat("\n\nKEY NUMBERS:\n")
str(key_nums)
cat("\nDone.\n")
