## 05_mirage_analyses.R
## NEW 2026-10-05 (v3, "The Momentum Mirage"): six analyses added on top of the
## v2 pipeline. Reads data/drives_v2.rds plus the nflfastR play-by-play (for
## swing-moment flags, first-scrimmage-play EP, and kickoff types).
##
##   A1. Gelbach (2016) decomposition of the raw H1 gap
##   A2. Equivalence tests (TOST) for H1-H3 at +/-1 and +/-2 pp
##   A3. The equivalence bound in yards of field position and expected points
##   A4. Momentum-swing moments (takeaways, 4th-down stops, long TDs, non-offensive TDs)
##   A5. Cost of chasing momentum: onside and short kickoffs vs standard kickoffs
##   A6. Game fixed effects with few drives per game: simulated no-momentum null
##
## Outputs (all *_v3; v1/v2 outputs untouched):
##   data/key_nums_v3.rds, data/key_nums_v3.json
##   tables/tab_decomp_v3.tex, tab_tost_v3.tex, tab_swing_v3.tex,
##   tables/tab_kickoff_v3.tex, tab_nickell_v3.tex
##   figures/fig3_decomp_v3.{pdf,png}

suppressPackageStartupMessages({
  library(data.table)
  library(fixest)
  library(ggplot2)
  library(jsonlite)
  library(here)
})

set.seed(20261005)
DATA_DIR <- here::here("data")
TAB_DIR  <- here::here("tables")
FIG_DIR  <- here::here("figures")

N_BOOT <- 200   # cluster bootstrap reps for the decomposition
N_SIM <- 100    # simulation reps for A6

## ---------------------------------------------------------------------------
## Data: v2 drives, same sample rules as 02_analyze_momentum.R
## ---------------------------------------------------------------------------
d_all <- readRDS(file.path(DATA_DIR, "drives_v2.rds"))
setDT(d_all)
d <- d_all[!is.na(off_success) & !is.na(yardline_100) & !is.na(score_differential)]
d[, def_success := fifelse(is.na(def_success), 0L, def_success)]

ctrl <- "yardline_100 + score_differential + qtr + half_seconds_remaining + home_off"
fe_head <- "posteam + defteam + season"
cl <- ~game_id + posteam

d_main <- d[!is.na(prior_own_def_success) & !is.na(prior_own_off_success) &
            !is.na(prior_opp_off_success) & !is.na(half_seconds_remaining)]
cat("Joint estimation sample:", nrow(d_main), "\n")

f_head <- as.formula(paste("off_success ~ prior_own_def_success + prior_own_off_success +",
                           "prior_opp_off_success +", ctrl, "|", fe_head))
m_head <- feols(f_head, data = d_main, cluster = cl)

## ---------------------------------------------------------------------------
## Play-by-play: drive-level flags (same filter and drive keys as 01_build_drives.R)
## ---------------------------------------------------------------------------
cat("Loading 2010-2024 pbp...\n")
pbp <- nflreadr::load_pbp(2010:2024)
setDT(pbp)
pbp <- pbp[!is.na(posteam) & !is.na(defteam) & !is.na(fixed_drive)]
setorder(pbp, game_id, fixed_drive, play_id)

scrim <- pbp[!is.na(play_type) & !(play_type %in% c("kickoff", "no_play")) &
               kickoff_attempt == 0 & !is.na(yardline_100)]
ep_scrim <- scrim[, .(ep_scrim = ep[1]), by = .(game_id, fixed_drive)]

drv_flags <- pbp[, .(
  long_td = as.integer(any(touchdown == 1 & td_team == posteam & kickoff_attempt == 0 &
                             punt_attempt == 0 & yards_gained >= 40, na.rm = TRUE)),
  kr_td   = as.integer(any(kickoff_attempt == 1 & return_touchdown == 1 &
                             td_team == posteam, na.rm = TRUE)),
  ## scrimmage snaps only (v2's `plays` counts every pbp row, incl. kickoffs,
  ## timeouts and end-of-quarter rows, so few punts met its "3 or fewer" rule)
  n_scrim = sum(play_type %in% c("pass", "run", "qb_kneel", "qb_spike"), na.rm = TRUE)
), by = .(game_id, fixed_drive)]
setorder(drv_flags, game_id, fixed_drive)
drv_flags[, `:=`(prev_long_td = shift(long_td), prev_kr_td = shift(kr_td),
                 prev_n_scrim = shift(n_scrim)), by = game_id]

d <- merge(d, drv_flags[, .(game_id, fixed_drive, long_td, prev_long_td, prev_kr_td,
                            prev_n_scrim)],
           by = c("game_id", "fixed_drive"), all.x = TRUE)
## NEW v3: H1 with three-and-outs counted on scrimmage snaps
d[, h1_scrim := fifelse(!is.na(prev_posteam) & prev_posteam != posteam,
                        as.integer(prev_drive_result %in% c("Turnover", "Turnover on downs", "Safety") |
                                     (prev_drive_result == "Punt" & prev_n_scrim <= 3)),
                        NA_integer_)]
d_main <- merge(d_main, d[, .(game_id, fixed_drive, h1_scrim)],
                by = c("game_id", "fixed_drive"), all.x = TRUE)
d <- merge(d, ep_scrim, by = c("game_id", "fixed_drive"), all.x = TRUE)
d_main <- merge(d_main, ep_scrim, by = c("game_id", "fixed_drive"), all.x = TRUE)

key <- list()
key$n_joint <- nobs(m_head)
key$h1 <- unname(coef(m_head)["prior_own_def_success"])
key$h2 <- unname(coef(m_head)["prior_own_off_success"])
key$h3 <- unname(coef(m_head)["prior_opp_off_success"])
key$se_h1 <- unname(se(m_head)["prior_own_def_success"])
key$se_h2 <- unname(se(m_head)["prior_own_off_success"])
key$se_h3 <- unname(se(m_head)["prior_opp_off_success"])
key$fp_slope <- unname(coef(m_head)["yardline_100"])

## ===========================================================================
## A1. Gelbach (2016) decomposition of the raw H1 gap
## ===========================================================================
## beta_base - beta_full = sum_g delta_g, delta_g = coef of D in OLS of
## (X_g %*% b_g) on D, i.e. the stop-vs-no-stop mean difference of each group's
## fitted contribution. Exact with a common sample.
gelbach <- function(dt, full_rhs, groups, fe) {
  f_full <- as.formula(paste("off_success ~ prior_own_def_success +", full_rhs, "|", fe))
  mf <- feols(f_full, data = dt, cluster = cl)
  b <- coef(mf)
  s <- dt[obs(mf)]
  D <- s$prior_own_def_success
  gap <- mean(s$off_success[D == 1]) - mean(s$off_success[D == 0])
  dif <- function(h) mean(h[D == 1]) - mean(h[D == 0])
  out <- sapply(groups, function(vars) {
    h <- as.matrix(s[, ..vars]) %*% b[vars]
    dif(as.numeric(h))
  })
  fx <- fixef(mf)
  for (nm in names(fx)) {
    h <- fx[[nm]][as.character(s[[nm]])]
    out[paste0("FE: ", nm)] <- dif(as.numeric(h))
  }
  list(gap = gap, full = unname(b["prior_own_def_success"]), contrib = out, n = nobs(mf))
}

grp_h1 <- list("Field position" = "yardline_100",
               "Score state"    = "score_differential",
               "Time remaining" = c("qtr", "half_seconds_remaining"),
               "Home"           = "home_off")
grp_joint <- c(grp_h1, list("Other momentum indicators" =
                              c("prior_own_off_success", "prior_opp_off_success")))

d_h1 <- d[!is.na(prior_own_def_success) & !is.na(half_seconds_remaining)]
g_h1 <- gelbach(d_h1, ctrl, grp_h1, fe_head)
g_jt <- gelbach(d_main, paste("prior_own_off_success + prior_opp_off_success +", ctrl),
                grp_joint, fe_head)
cat("Gelbach H1-only: gap", round(g_h1$gap, 4), "full", round(g_h1$full, 4),
    "sum check", round(g_h1$gap - sum(g_h1$contrib) - g_h1$full, 8), "\n")
print(round(g_h1$contrib, 4))
cat("Gelbach joint: gap", round(g_jt$gap, 4), "full", round(g_jt$full, 4),
    "sum check", round(g_jt$gap - sum(g_jt$contrib) - g_jt$full, 8), "\n")
print(round(g_jt$contrib, 4))

## Game-cluster bootstrap SEs for the H1-only decomposition
games <- unique(d_h1$game_id)
setkey(d_h1, game_id)
boot <- replicate(N_BOOT, {
  gs <- sample(games, length(games), replace = TRUE)
  bd <- d_h1[J(gs), allow.cartesian = TRUE]
  r <- gelbach(bd, ctrl, grp_h1, fe_head)
  c(gap = r$gap, full = r$full, r$contrib)
})
boot_se <- apply(boot, 1, sd)

decomp_tab <- data.table(
  component = c("Raw gap (stop minus no stop)", names(g_h1$contrib),
                "Residual (H1 coefficient, full spec)"),
  estimate = c(g_h1$gap, g_h1$contrib, g_h1$full),
  se = c(boot_se["gap"], boot_se[names(g_h1$contrib)], boot_se["full"]))
decomp_tab[, share := estimate / g_h1$gap]
print(decomp_tab)

key$decomp_n <- g_h1$n
key$decomp_gap <- g_h1$gap
key$decomp_full <- g_h1$full
for (nm in names(g_h1$contrib)) key[[paste0("decomp_", nm)]] <- unname(g_h1$contrib[nm])
key$decomp_share_fp <- unname(g_h1$contrib["Field position"] / g_h1$gap)
key$decomp_share_resid <- g_h1$full / g_h1$gap
key$decomp_joint_gap <- g_jt$gap
key$decomp_joint_fp <- unname(g_jt$contrib["Field position"])
key$decomp_joint_share_fp <- unname(g_jt$contrib["Field position"] / g_jt$gap)
key$decomp_joint_full <- g_jt$full

pretty_comp <- function(x) {
  x <- sub("FE: posteam", "Possessing-team FE", x)
  x <- sub("FE: defteam", "Opponent FE", x)
  sub("FE: season", "Season FE", x)
}
fmt <- function(x, k = 1) formatC(x, format = "f", digits = k)
rows <- decomp_tab[, sprintf("%s & %s & (%s) & %s \\\\", pretty_comp(component),
                            fmt(100 * estimate), fmt(100 * se), fmt(100 * share, 0))]
writeLines(c(
  "\\begin{table}[htbp]", "\\centering",
  paste0("\\caption{\\label{tab:decomp} Where the defensive-stop gap goes: Gelbach (2016) ",
         "decomposition. Percentage points. The raw gap is the difference in scoring ",
         "rates after a defensive stop versus any other opponent outcome (N = ",
         format(g_h1$n, big.mark = ","), " drives with the H1 indicator and time ",
         "remaining). Each row is the part of the gap explained by that block of the ",
         "headline specification (field position, score differential, quarter and ",
         "half-seconds remaining, home, and possessing-team, opponent, and season fixed ",
         "effects). Bootstrap standard errors (", N_BOOT, " game-cluster resamples) in ",
         "parentheses. Shares are of the raw gap.}"),
  "\\begin{tabular}{lrrr}", "\\toprule",
  "Component & Estimate (pp) & (s.e.) & Share (\\%) \\\\", "\\midrule",
  rows[1], "\\midrule", rows[2:(length(rows) - 1)], "\\midrule", rows[length(rows)],
  "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  file.path(TAB_DIR, "tab_decomp_v3.tex"))

plot_dt <- decomp_tab[2:(.N - 1)]
plot_dt[, component := pretty_comp(component)]
plot_dt <- rbind(plot_dt, decomp_tab[.N][, component := "Residual (\"momentum\")"])
plot_dt[, component := factor(component, levels = rev(component))]
p3 <- ggplot(plot_dt, aes(x = 100 * estimate, y = component)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_col(width = 0.6, fill = "grey35") +
  geom_errorbarh(aes(xmin = 100 * (estimate - 1.96 * se), xmax = 100 * (estimate + 1.96 * se)),
                 height = 0.2) +
  labs(x = "Percentage points of the raw defensive-stop gap", y = NULL,
       title = sprintf("Decomposing the %.1f-pp gap after a defensive stop", 100 * g_h1$gap),
       subtitle = "Gelbach decomposition; 95% game-cluster bootstrap intervals") +
  theme_minimal(base_size = 11)
ggsave(file.path(FIG_DIR, "fig3_decomp_v3.pdf"), p3, width = 7, height = 3.5)
ggsave(file.path(FIG_DIR, "fig3_decomp_v3.png"), p3, width = 7, height = 3.5, dpi = 200)

## ===========================================================================
## A2. TOST equivalence tests (alpha = .05: 90% CI inside +/- bound)
## A3. Bound in yards and expected points
## ===========================================================================
z <- qnorm(0.95)
ep_fit <- feols(ep_scrim ~ yardline_100, data = d_main)
ep_per_yard <- -unname(coef(ep_fit)["yardline_100"])   # EP gained per yard closer
yards_per_pp <- 1 / abs(100 * key$fp_slope)           # yards that move P(score) by 1 pp
key$ep_per_yard <- ep_per_yard
key$yards_per_pp <- yards_per_pp

tost <- rbindlist(lapply(1:3, function(j) {
  b <- c(key$h1, key$h2, key$h3)[j]; s <- c(key$se_h1, key$se_h2, key$se_h3)[j]
  p_at <- function(bound) max(pnorm((b - (-bound)) / s, lower.tail = FALSE),
                              pnorm((b - bound) / s))
  minb <- abs(b) + z * s
  data.table(channel = c("H1", "H2", "H3")[j], est = b, se = s,
             ci90_lo = b - z * s, ci90_hi = b + z * s,
             p_tost_2 = p_at(0.02), p_tost_1 = p_at(0.01),
             min_bound = minb,
             min_bound_yards = 100 * minb * yards_per_pp,
             min_bound_ep = 100 * minb * yards_per_pp * ep_per_yard)
}))
print(tost)
for (j in 1:3) {
  ch <- tolower(tost$channel[j])
  key[[paste0("tost_", ch, "_p2")]] <- tost$p_tost_2[j]
  key[[paste0("tost_", ch, "_p1")]] <- tost$p_tost_1[j]
  key[[paste0("tost_", ch, "_minbound")]] <- tost$min_bound[j]
  key[[paste0("tost_", ch, "_yards")]] <- tost$min_bound_yards[j]
  key[[paste0("tost_", ch, "_ep")]] <- tost$min_bound_ep[j]
  key[[paste0("ci90_", ch, "_lo")]] <- tost$ci90_lo[j]
  key[[paste0("ci90_", ch, "_hi")]] <- tost$ci90_hi[j]
}
pf <- function(p) ifelse(p < 0.001, "$<$0.001", fmt(p, 3))
tost_rows <- tost[, sprintf("%s & %s (%s) & [%s, %s] & %s & %s & %s & %s & %s \\\\",
                            channel, fmt(100 * est), fmt(100 * se), fmt(100 * ci90_lo),
                            fmt(100 * ci90_hi), pf(p_tost_2), pf(p_tost_1),
                            fmt(100 * min_bound, 2), fmt(min_bound_yards, 2),
                            fmt(min_bound_ep, 2))]
writeLines(c(
  "\\begin{table}[htbp]", "\\centering",
  paste0("\\caption{\\label{tab:tost} Equivalence tests (two one-sided tests, $\\alpha=0.05$). ",
         "Estimates from the headline specification (column 1 of Table \\ref{tab:allchan}). ",
         "$p_{\\pm b}$ is the TOST $p$-value against the equivalence bound $\\pm b$ ",
         "percentage points; $p<0.05$ means the effect is shown to lie inside the bound. ",
         "The smallest bound is the tightest $\\pm b$ the 90\\% interval rules out. ",
         "Yards convert that bound with the model's own field-position slope (",
         fmt(abs(100 * key$fp_slope), 2), " pp per yard); expected points (EP) convert ",
         "the yards with the slope of \\texttt{nflfastR} EP at the first scrimmage play (",
         fmt(ep_per_yard, 3), " EP per yard).}"),
  "\\begin{tabular}{lccccccc}", "\\toprule",
  " & Est. (s.e.), pp & 90\\% CI & $p_{\\pm 2}$ & $p_{\\pm 1}$ & Smallest bound (pp) & Yards & EP \\\\",
  "\\midrule", tost_rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  file.path(TAB_DIR, "tab_tost_v3.tex"))

## NEW v3: H1 re-estimated with three-and-outs counted on scrimmage snaps
key$n_stop_v2 <- d[prior_own_def_success == 1, .N]
key$n_punt_3out_v2 <- d[prior_own_def_success == 1 & prev_drive_result == "Punt", .N]
key$n_stop_scrim <- d[h1_scrim == 1, .N]
key$n_punt_3out_scrim <- d[h1_scrim == 1 & prev_drive_result == "Punt", .N]
key$raw_gap_scrim <- d[h1_scrim == 1, mean(off_success)] - d[h1_scrim == 0, mean(off_success)]
key$p_score_stop_scrim <- d[h1_scrim == 1, mean(off_success)]
key$p_score_nostop_scrim <- d[h1_scrim == 0, mean(off_success)]
m_head_c <- feols(as.formula(paste("off_success ~ h1_scrim + prior_own_off_success +",
                                   "prior_opp_off_success +", ctrl, "|", fe_head)),
                  data = d_main, cluster = cl)
print(summary(m_head_c))
key$h1c <- unname(coef(m_head_c)["h1_scrim"]); key$se_h1c <- unname(se(m_head_c)["h1_scrim"])
key$h2c <- unname(coef(m_head_c)["prior_own_off_success"]); key$se_h2c <- unname(se(m_head_c)["prior_own_off_success"])
key$h3c <- unname(coef(m_head_c)["prior_opp_off_success"]); key$se_h3c <- unname(se(m_head_c)["prior_opp_off_success"])
key$h1c_ci90_lo <- key$h1c - z * key$se_h1c
key$h1c_ci90_hi <- key$h1c + z * key$se_h1c
key$tost_h1c_minbound <- abs(key$h1c) + z * key$se_h1c
key$tost_h1c_p1 <- max(pnorm((key$h1c + 0.01) / key$se_h1c, lower.tail = FALSE),
                       pnorm((key$h1c - 0.01) / key$se_h1c))
key$tost_h1c_p2 <- max(pnorm((key$h1c + 0.02) / key$se_h1c, lower.tail = FALSE),
                       pnorm((key$h1c - 0.02) / key$se_h1c))
key$tost_h1c_yards <- 100 * key$tost_h1c_minbound * yards_per_pp
key$tost_h1c_ep <- key$tost_h1c_yards * ep_per_yard
cat(sprintf("Corrected H1: stops %d (punts %d) vs v2 %d (punts %d); est %.4f (%.4f)\n",
            key$n_stop_scrim, key$n_punt_3out_scrim, key$n_stop_v2, key$n_punt_3out_v2,
            key$h1c, key$se_h1c))
tost_c <- rbindlist(lapply(c("h1c", "h2c", "h3c"), function(ch) {
  b <- key[[ch]]; s <- key[[paste0("se_", ch)]]
  p_at <- function(bound) max(pnorm((b + bound) / s, lower.tail = FALSE), pnorm((b - bound) / s))
  minb <- abs(b) + z * s
  data.table(channel = ch, est = b, se = s, lo = b - z * s, hi = b + z * s,
             p2 = p_at(0.02), p1 = p_at(0.01), minb = minb,
             yards = 100 * minb * yards_per_pp, ep = 100 * minb * yards_per_pp * ep_per_yard)
}))
for (j in 2:3) {
  ch <- tost_c$channel[j]
  key[[paste0(ch, "_ci90_lo")]] <- tost_c$lo[j]; key[[paste0(ch, "_ci90_hi")]] <- tost_c$hi[j]
  key[[paste0("tost_", ch, "_minbound")]] <- tost_c$minb[j]
  key[[paste0("tost_", ch, "_p1")]] <- tost_c$p1[j]; key[[paste0("tost_", ch, "_p2")]] <- tost_c$p2[j]
  key[[paste0("tost_", ch, "_yards")]] <- tost_c$yards[j]; key[[paste0("tost_", ch, "_ep")]] <- tost_c$ep[j]
}
tost_c_row <- tost_c[, sprintf("%s & %s (%s) & [%s, %s] & %s & %s & %s & %s & %s \\\\",
                               c("H1 (snap-count stops)", "H2 (same model)", "H3 (same model)"),
                               fmt(100 * est), fmt(100 * se), fmt(100 * lo), fmt(100 * hi),
                               pf(p2), pf(p1), fmt(100 * minb, 2), fmt(yards, 2), fmt(ep, 2))]
tl <- readLines(file.path(TAB_DIR, "tab_tost_v3.tex"))
i <- grep("^\\\\bottomrule", tl)
tl <- append(tl, c("\\midrule", "\\multicolumn{8}{l}{\\emph{Stops counted on scrimmage snaps (three-and-out punts included)}} \\\\", tost_c_row), after = i - 1)
writeLines(tl, file.path(TAB_DIR, "tab_tost_v3.tex"))

## ===========================================================================
## A4. Momentum-swing moments
## ===========================================================================
## Classify the drive immediately before each drive (game order, before the
## end-of-half drops). Opponent possession categories are mutually exclusive;
## reference = an opponent drive that ended in a punt after 4+ plays, a missed
## field goal, or another non-scoring, non-stop outcome.
d[, prev_cat := fcase(
  is.na(prev_posteam), NA_character_,
  prev_posteam == posteam & prev_drive_result == "Opp touchdown", "nonoff_td_against",
  prev_posteam == posteam, "same_team_other",
  prev_drive_result == "Turnover", "takeaway",
  prev_drive_result == "Turnover on downs", "fourth_down_stop",
  prev_drive_result == "Safety" | (prev_drive_result == "Punt" & prev_n_scrim <= 3), "three_and_out",
  prev_drive_result == "Touchdown" & prev_kr_td == 1, "nonoff_td_against",
  prev_drive_result == "Touchdown" & prev_long_td == 1, "long_td_against",
  prev_drive_result == "Touchdown", "short_td_against",
  prev_drive_result == "Field goal", "fg_against",
  prev_drive_result == "Opp touchdown", "own_def_td",   # our defense scored; rare (we receive only after onside etc.)
  default = "ref")]
print(d[, .N, by = prev_cat][order(-N)])

d_sw <- d[!is.na(prev_cat) & !(prev_cat %in% c("same_team_other", "own_def_td")) &
            !is.na(half_seconds_remaining)]
m_sw <- feols(as.formula(paste("off_success ~ i(prev_cat, ref = 'ref') +", ctrl, "|", fe_head)),
              data = d_sw, cluster = cl)
print(summary(m_sw))
cf <- coef(m_sw); V <- vcov(m_sw)
nm <- function(x) paste0("prev_cat::", x)
contrast <- function(a, b) {
  v <- setNames(rep(0, length(cf)), names(cf)); v[nm(a)] <- 1; if (!is.null(b)) v[nm(b)] <- -1
  c(est = sum(v * cf), se = sqrt(drop(t(v) %*% V %*% v)))
}
sw_rows <- list(
  takeaway = contrast("takeaway", NULL),
  fourth_down_stop = contrast("fourth_down_stop", NULL),
  three_and_out = contrast("three_and_out", NULL),
  fg_against = contrast("fg_against", NULL),
  short_td_against = contrast("short_td_against", NULL),
  long_td_against = contrast("long_td_against", NULL),
  nonoff_td_against = contrast("nonoff_td_against", NULL),
  long_vs_short_td = contrast("long_td_against", "short_td_against"),
  nonoff_vs_short_td = contrast("nonoff_td_against", "short_td_against"),
  takeaway_vs_3out = contrast("takeaway", "three_and_out"),
  fourth_vs_3out = contrast("fourth_down_stop", "three_and_out"))
sw_n <- d_sw[, .N, by = prev_cat]

## Beneficiary side: the team that made the swing play, on ITS next offensive drive.
## own_long_td: own previous offensive drive was a TD of >= 40 yards on the scoring play
## (vs. any other own TD); own_def_td: the opponent drive just before ended in our
## defensive/return TD (the opponent then had a drive in between).
d <- d[order(game_id, fixed_drive)]
d[, own_last_long := {
  res <- rep(NA_integer_, .N); last <- list()
  for (i in seq_len(.N)) {
    tm <- posteam[i]
    if (!is.null(last[[tm]])) res[i] <- last[[tm]]
    last[[tm]] <- long_td[i]
  }
  res
}, by = game_id]
d[, own_def_td_since := {
  ## did this team's defense/returners score a non-offensive TD since its last possession?
  res <- rep(0L, .N); flag <- list()
  for (i in seq_len(.N)) {
    tm <- posteam[i]; op <- defteam[i]
    res[i] <- if (!is.null(flag[[tm]])) flag[[tm]] else 0L
    flag[[tm]] <- 0L
    if (!is.na(drive_result[i]) && drive_result[i] == "Opp touchdown") flag[[op]] <- 1L
  }
  res
}, by = game_id]
d_ben <- d[!is.na(prior_own_off_result) & !is.na(half_seconds_remaining)]
d_ben[, own_td := as.integer(prior_own_off_result == "Touchdown")]
d_ben[, own_long := as.integer(own_td == 1 & own_last_long == 1)]
m_ben <- feols(as.formula(paste("off_success ~ own_td + own_long + own_def_td_since +",
                                "prior_opp_off_success +", ctrl, "|", fe_head)),
               data = d_ben[!is.na(prior_opp_off_success)], cluster = cl)
print(summary(m_ben))
ben_rows <- list(
  own_long_vs_other_td = c(est = unname(coef(m_ben)["own_long"]), se = unname(se(m_ben)["own_long"])),
  own_def_td = c(est = unname(coef(m_ben)["own_def_td_since"]),
                 se = unname(se(m_ben)["own_def_td_since"])))
key$n_swing <- nobs(m_sw)
key$n_benef <- nobs(m_ben)
key$n_own_long <- sum(d_ben[!is.na(prior_opp_off_success)]$own_long, na.rm = TRUE)
key$n_own_def_td <- sum(d_ben[!is.na(prior_opp_off_success)]$own_def_td_since)
for (k in names(sw_rows)) {
  key[[paste0("sw_", k)]] <- unname(sw_rows[[k]]["est"])
  key[[paste0("sw_", k, "_se")]] <- unname(sw_rows[[k]]["se"])
}
for (k in names(ben_rows)) {
  key[[paste0("ben_", k)]] <- unname(ben_rows[[k]]["est"])
  key[[paste0("ben_", k, "_se")]] <- unname(ben_rows[[k]]["se"])
}
for (k in sw_n$prev_cat) key[[paste0("n_prev_", k)]] <- sw_n[prev_cat == k]$N

lab <- c(takeaway = "Takeaway (INT or fumble lost)",
         fourth_down_stop = "Fourth-down stop",
         three_and_out = "Punt after $\\leq$3 scrimmage snaps, or safety",
         fg_against = "Opponent field goal",
         short_td_against = "Opponent offensive TD, scoring play $<$40 yds",
         long_td_against = "Opponent offensive TD, scoring play $\\geq$40 yds",
         nonoff_td_against = "Opponent defensive or return TD",
         long_vs_short_td = "\\quad Long TD minus short TD",
         nonoff_vs_short_td = "\\quad Defensive/return TD minus short TD",
         takeaway_vs_3out = "\\quad Takeaway minus three-and-out",
         fourth_vs_3out = "\\quad Fourth-down stop minus three-and-out",
         own_long_vs_other_td = "Own prior drive was a $\\geq$40-yd TD (vs.\\ other own TD)",
         own_def_td = "Own defense or returners just scored a TD")
nfor <- function(k) { x <- sw_n[prev_cat == k]$N; if (length(x)) format(x, big.mark = ",") else "" }
mkrow <- function(k, r, n = "") sprintf("%s & %s & (%s) & %s \\\\", lab[k], fmt(100 * r["est"]),
                                       fmt(100 * r["se"]), n)
writeLines(c(
  "\\begin{table}[htbp]", "\\centering",
  paste0("\\caption{\\label{tab:swing} Scoring after momentum-swing moments. Percentage-point ",
         "change in the probability that the next drive ends in a score, from linear ",
         "probability models with the headline controls (field position at the first ",
         "scrimmage play, score differential, quarter, half-seconds remaining, home) and ",
         "possessing-team, opponent, and season fixed effects. Panel A: categories of the ",
         "opponent possession immediately before the drive; reference category is an ",
         "opponent drive ending in a punt after four or more plays, a missed field goal, or ",
         "another non-scoring outcome (N = ", format(nobs(m_sw), big.mark = ","), "). ",
         "Panel B: the team that made the big play, on its next offensive drive, also ",
         "controlling for whether its previous drive was a touchdown and whether the opponent ",
         "just scored (N = ", format(nobs(m_ben), big.mark = ","), "). Standard errors ",
         "two-way clustered by game and possessing team.}"),
  "\\begin{tabular}{lrrr}", "\\toprule",
  "Preceding event & Est. (pp) & (s.e.) & N \\\\", "\\midrule",
  "\\multicolumn{4}{l}{\\emph{A. Drive right after the event}} \\\\",
  mapply(function(k) mkrow(k, sw_rows[[k]], nfor(k)), names(sw_rows)[1:7]),
  "\\multicolumn{4}{l}{\\emph{Contrasts}} \\\\",
  mapply(function(k) mkrow(k, sw_rows[[k]]), names(sw_rows)[8:11]),
  "\\midrule", "\\multicolumn{4}{l}{\\emph{B. The play-maker's next offensive drive}} \\\\",
  mkrow("own_long_vs_other_td", ben_rows$own_long_vs_other_td, format(key$n_own_long, big.mark = ",")),
  mkrow("own_def_td", ben_rows$own_def_td, format(key$n_own_def_td, big.mark = ",")),
  "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  file.path(TAB_DIR, "tab_swing_v3.tex"))

## ===========================================================================
## A5. Cost of chasing momentum: kickoff type after scores
## ===========================================================================
## nflfastR labels onside kicks in the play description ("kicks onside"). It
## does NOT label squib kicks (zero descriptions contain "squib"), so short
## non-onside kicks are identified by where the kick lands: at or beyond the
## receiving team's 15 (landing = kick spot distance to goal - kick_distance).
## That class mixes squibs and pooch kicks; it is a proxy, not a squib label.
## Sample: 2010-2023 (the 2024 dynamic kickoff changed kick rules), regulation,
## post-score kickoffs (not half-opening, not safety free kicks), and only
## "discretionary" game states: excluding the final 5 minutes of Q4 and the
## final 2 minutes of Q2, where onside/short kicks are clock-driven.
ko <- pbp[kickoff_attempt == 1 & season <= 2023 & qtr <= 4 & !is.na(kick_distance) &
            !is.na(epa) & yardline_100 != 20]   # spot 20 = safety free kick (receiver view)
## half-opening kickoffs: the first kickoff of each half in each game
ko[, half_id := fifelse(qtr <= 2, 1L, 2L)]
setorder(ko, game_id, play_id)
ko[, first_of_half := seq_len(.N) == 1L, by = .(game_id, half_id)]
ko <- ko[first_of_half == FALSE]
ko[, onside := grepl("onside", desc, ignore.case = TRUE)]
ko[, landing := (100 - yardline_100) - kick_distance]
ko[, ktype := fifelse(onside, "onside", fifelse(landing >= 15, "short", "standard"))]
ko[, kick_epa := -epa]                    # kicking team's perspective
ko[, kick_margin := -score_differential]  # kicking team's lead
ko[, discretionary := !(qtr == 4 & game_seconds_remaining <= 300) &
       !(qtr == 2 & half_seconds_remaining <= 120)]
print(ko[, .N, by = .(discretionary, ktype)][order(discretionary, ktype)])
kd <- ko[discretionary == TRUE]
kd[, ktype := relevel(factor(ktype), ref = "standard")]
m_ko <- feols(kick_epa ~ ktype + kick_margin + qtr | season, data = kd, cluster = ~game_id)
print(summary(m_ko))
ko_sum <- kd[, .(n = .N, mean_epa = mean(kick_epa),
                 recov = mean(own_kickoff_recovery == 1, na.rm = TRUE)), by = ktype]
print(ko_sum)
key$ko_n_total <- nrow(kd)
for (t in c("standard", "short", "onside")) {
  key[[paste0("ko_n_", t)]] <- ko_sum[ktype == t]$n
  key[[paste0("ko_mean_epa_", t)]] <- ko_sum[ktype == t]$mean_epa
  key[[paste0("ko_recov_", t)]] <- ko_sum[ktype == t]$recov
}
key$ko_cost_short <- unname(coef(m_ko)["ktypeshort"])
key$ko_cost_short_se <- unname(se(m_ko)["ktypeshort"])
key$ko_cost_onside <- unname(coef(m_ko)["ktypeonside"])
key$ko_cost_onside_se <- unname(se(m_ko)["ktypeonside"])
## Clock-driven states, for contrast
kc <- ko[discretionary == FALSE]
kc[, ktype := relevel(factor(ktype), ref = "standard")]
m_ko_c <- feols(kick_epa ~ ktype + kick_margin + qtr | season, data = kc, cluster = ~game_id)
key$ko_cost_onside_late <- unname(coef(m_ko_c)["ktypeonside"])
key$ko_cost_onside_late_se <- unname(se(m_ko_c)["ktypeonside"])
key$ko_n_onside_late <- kc[ktype == "onside", .N]
## Break-even: momentum effect needed on the kicking team's next drive to repay
## the cost, in pp of P(score), using mean points per scoring drive in the sample.
pts_per_score <- d[off_success == 1, mean(fifelse(drive_result == "Touchdown", 7, 3))]
key$pts_per_score <- pts_per_score
key$ko_breakeven_short_pp <- -key$ko_cost_short / pts_per_score
key$ko_breakeven_onside_pp <- -key$ko_cost_onside / pts_per_score

ko_rows <- c(
  sprintf("Standard (lands inside the 15 or deeper) & %s & %s & --- & --- \\\\",
          format(key$ko_n_standard, big.mark = ","), fmt(key$ko_mean_epa_standard, 2)),
  sprintf("Short, non-onside (lands at the 15 or beyond) & %s & %s & %s & (%s) \\\\",
          format(key$ko_n_short, big.mark = ","), fmt(key$ko_mean_epa_short, 2),
          fmt(key$ko_cost_short, 2), fmt(key$ko_cost_short_se, 2)),
  sprintf("Onside & %s & %s & %s & (%s) \\\\",
          format(key$ko_n_onside, big.mark = ","), fmt(key$ko_mean_epa_onside, 2),
          fmt(key$ko_cost_onside, 2), fmt(key$ko_cost_onside_se, 2)))
writeLines(c(
  "\\begin{table}[htbp]", "\\centering",
  paste0("\\caption{\\label{tab:kickoff} The price of chasing momentum on the kickoff. ",
         "Kicking team's expected points added (EPA, \\texttt{nflfastR}) on post-score ",
         "kickoffs, 2010--2023 regulation, excluding the final five minutes of the fourth ",
         "quarter and final two minutes of the second quarter. Onside kicks are labeled in the ",
         "play-by-play; squib kicks are not, so short kicks are identified by landing spot and ",
         "include pooch kicks. Difference vs.\\ standard: OLS with kicking-team score margin, ",
         "quarter, and season fixed effects; standard errors clustered by game.}"),
  "\\begin{tabular}{lrrrr}", "\\toprule",
  "Kick type & N & Mean EPA & Diff.\\ vs.\\ standard & (s.e.) \\\\", "\\midrule",
  ko_rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  file.path(TAB_DIR, "tab_kickoff_v3.tex"))

## ===========================================================================
## A6. Game fixed effects with ~20 drives per game: simulated null
## ===========================================================================
## Simulate drive outcomes with NO momentum: each drive scores with probability
## p = p0 + u, where p0 is the fitted value from the headline controls and FE
## (no momentum terms) and u is a team-game shock. The shocks' variance and the
## covariance between the two teams in a game are calibrated by method of
## moments to the observed team-game mean residuals. Stops are drawn among
## non-scoring drives at the observed rate. Lags and the lead placebo are built
## exactly as for the observed data, and both specifications are estimated.
dp <- d[!is.na(half_seconds_remaining)][order(game_id, fixed_drive)]
build_lags_fast <- function(x) {
  x[, prev_pt := shift(posteam), by = game_id]
  x[, l_def := fifelse(!is.na(prev_pt) & prev_pt != posteam, shift(def_s), NA_integer_), by = game_id]
  x[, h_out := fifelse(posteam == home_team, off_s, NA_integer_)]
  x[, a_out := fifelse(posteam != home_team, off_s, NA_integer_)]
  x[, last_h := nafill(shift(h_out), "locf"), by = game_id]
  x[, last_a := nafill(shift(a_out), "locf"), by = game_id]
  x[, l_own := fifelse(posteam == home_team, last_h, last_a)]
  x[, l_opp := fifelse(posteam == home_team, last_a, last_h)]
  x[, next_pt := shift(posteam, type = "lead"), by = game_id]
  x[, l_lead := fifelse(!is.na(next_pt) & next_pt != posteam,
                        shift(off_s, type = "lead"), NA_integer_), by = game_id]
  x
}
fit_both <- function(x) {
  s <- x[!is.na(l_def) & !is.na(l_own) & !is.na(l_opp) & !is.na(l_lead)]
  f1 <- as.formula(paste("off_s ~ l_def + l_own + l_opp + l_lead +", ctrl, "|", fe_head))
  f2 <- as.formula(paste("off_s ~ l_def + l_own + l_opp + l_lead +",
                         "yardline_100 + score_differential + qtr + half_seconds_remaining",
                         "| game_id + posteam"))
  c(head = coef(feols(f1, data = s))[c("l_def", "l_own", "l_opp", "l_lead")],
    game = coef(feols(f2, data = s))[c("l_def", "l_own", "l_opp", "l_lead")])
}
dp[, `:=`(off_s = off_success, def_s = def_success)]
real <- fit_both(build_lags_fast(copy(dp)))

m0 <- feols(as.formula(paste("off_success ~", ctrl, "|", fe_head)), data = dp)
dp[, p0 := fitted(m0)]
dp[, home_side := as.integer(posteam == home_team)]
tg <- dp[, .(m = mean(off_success - p0), v = mean(p0 * (1 - p0)) / .N), by = .(game_id, home_side)]
sig2_u <- max(tg[, mean(m^2) - mean(v)], 0)
tgw <- dcast(tg, game_id ~ home_side, value.var = "m")
cov_u <- tgw[, mean(`0` * `1`, na.rm = TRUE)]
q_stop <- dp[off_success == 0, mean(def_success)]
key$sim_sd_u <- sqrt(sig2_u)
key$sim_rho_u <- cov_u / sig2_u
cat("Null sim: sd(u) =", round(sqrt(sig2_u), 4), " rho =", round(cov_u / sig2_u, 3), "\n")
Sig <- matrix(c(sig2_u, cov_u, cov_u, sig2_u), 2)
L <- chol(Sig)
gids <- unique(dp$game_id)
sim <- replicate(N_SIM, {
  x <- copy(dp)
  u <- matrix(rnorm(2 * length(gids)), ncol = 2) %*% L
  ug <- data.table(game_id = gids, u0 = u[, 1], u1 = u[, 2])
  x <- merge(x, ug, by = "game_id")[order(game_id, fixed_drive)]
  x[, p := pmin(pmax(p0 + fifelse(home_side == 1L, u1, u0), 0.001), 0.999)]
  x[, off_s := rbinom(.N, 1, p)]
  x[, def_s := fifelse(off_s == 1L, 0L, rbinom(.N, 1, q_stop))]
  fit_both(build_lags_fast(x))
})
sim_mean <- rowMeans(sim); sim_sd <- apply(sim, 1, sd)
print(round(rbind(real = real, sim_mean = sim_mean, sim_sd = sim_sd), 4))
drives_per_game <- dp[, .N, by = game_id][, mean(N)]
key$drives_per_game <- drives_per_game
key$nickell_approx <- -1 / (drives_per_game - 1)
key$n_sim <- N_SIM
for (nmx in names(sim_mean)) {
  k2 <- gsub("\\.", "_", nmx)
  key[[paste0("sim_", k2)]] <- unname(sim_mean[nmx])
  key[[paste0("sim_sd_", k2)]] <- unname(sim_sd[nmx])
  key[[paste0("real_", k2)]] <- unname(real[nmx])
  key[[paste0("adj_", k2)]] <- unname(real[nmx] - sim_mean[nmx])
}
nk_lab <- c(l_def = "H1: prior own defensive stop", l_own = "H2: prior own score",
            l_opp = "H3: prior opponent score", l_lead = "Placebo: opponent's next drive scores")
nk_rows <- sapply(c("l_def", "l_own", "l_opp", "l_lead"), function(v) {
  h <- paste0("head.", v); g <- paste0("game.", v)
  sprintf("%s & %s & %s (%s) & %s & %s (%s) \\\\", nk_lab[v],
          fmt(100 * real[h]), fmt(100 * sim_mean[h]), fmt(100 * sim_sd[h]),
          fmt(100 * real[g]), fmt(100 * sim_mean[g]), fmt(100 * sim_sd[g]))
})
writeLines(c(
  "\\begin{table}[htbp]", "\\centering",
  paste0("\\caption{\\label{tab:nickell} Game fixed effects manufacture negative momentum. ",
         "Coefficients (pp) on lagged and lead drive outcomes in the observed data and in ",
         N_SIM, " simulated seasons with no momentum: each drive scores with its fitted ",
         "probability from the headline controls plus a team-game shock (s.d.\\ ",
         fmt(sqrt(sig2_u), 3), ", correlation between the two teams in a game ",
         fmt(cov_u / sig2_u, 2), ", both calibrated to the data). ",
         "Each model includes all four indicators and the controls; the game-FE ",
         "specification uses game and possessing-team fixed effects. Sample: drives with ",
         "all four indicators defined. Simulation means, standard deviations in parentheses. ",
         "Mean drives per game: ", fmt(drives_per_game, 1), "; the textbook within-group ",
         "bias for one other observation in the same group is $-1/(T-1) = ",
         fmt(100 * key$nickell_approx), "$ pp.}"),
  "\\begin{tabular}{lrrrr}", "\\toprule",
  " & \\multicolumn{2}{c}{Headline FE (team, opp., season)} & \\multicolumn{2}{c}{Game + team FE} \\\\",
  "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
  "Regressor & Observed & Null sim. & Observed & Null sim. \\\\", "\\midrule",
  nk_rows, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
  file.path(TAB_DIR, "tab_nickell_v3.tex"))

## ---------------------------------------------------------------------------
saveRDS(key, file.path(DATA_DIR, "key_nums_v3.rds"))
write_json(key, file.path(DATA_DIR, "key_nums_v3.json"), auto_unbox = TRUE, digits = 8,
           pretty = TRUE)
cat("\nKEY NUMBERS v3:\n"); str(key)
cat("\nDone.\n")
