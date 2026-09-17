## 01_build_drives.R
## Build drive-level dataset from nflfastR play-by-play 2010-2024
## Output: data/drives.rds with one row per offensive drive

suppressPackageStartupMessages({
  library(nflfastR)
  library(data.table)
  library(dplyr)
  library(here)
})

set.seed(42)

DATA_DIR <- here::here("data")
dir.create(DATA_DIR, recursive = TRUE, showWarnings = FALSE)

cat("Loading 2010-2024 pbp...\n")
seasons <- 2010:2024
pbp <- load_pbp(seasons)
setDT(pbp)
cat("  rows:", nrow(pbp), "  games:", uniqueN(pbp$game_id), "\n")

## Keep regular-season + playoffs, both team possessions valid
pbp <- pbp[!is.na(posteam) & !is.na(defteam) & !is.na(fixed_drive)]

## Aggregate to drive level using fixed_drive (within game) -------------------
drive_keys <- c("game_id", "fixed_drive")

## First play of drive for state vars, last play for outcome
first_play <- pbp[order(game_id, fixed_drive, play_id),
                  .SD[1], by = drive_keys,
                  .SDcols = c("posteam","defteam","season","week","season_type",
                              "home_team","away_team","qtr","game_seconds_remaining",
                              "half_seconds_remaining","yardline_100",
                              "score_differential","wp","ep")]

last_play  <- pbp[order(game_id, fixed_drive, play_id),
                  .SD[.N], by = drive_keys,
                  .SDcols = c("fixed_drive_result","drive_play_count",
                              "drive_time_of_possession","play_id",
                              "drive_first_downs")]

## Plays per drive count and yards (sum yards_gained on rush+pass)
play_cnt <- pbp[, .(plays = .N,
                    drive_yards = sum(yards_gained, na.rm = TRUE)),
                by = drive_keys]

drives <- merge(first_play, last_play, by = drive_keys)
drives <- merge(drives, play_cnt, by = drive_keys)

## Outcome categories ---------------------------------------------------------
drives[, drive_result := fixed_drive_result]

## Offensive success (binary): TD or FG
drives[, off_success := as.integer(drive_result %in% c("Touchdown","Field goal"))]

## Defensive success: forced punt with <=3 plays (three-and-out),
## turnover (INT, fumble lost), turnover on downs, safety
drives[, def_success := as.integer(
  (drive_result == "Punt" & plays <= 3) |
  drive_result %in% c("Turnover","Turnover on downs","Safety")
)]

## EPA per drive (sum of play EPA on this drive, posteam perspective)
epa_drv <- pbp[, .(drive_epa = sum(epa, na.rm = TRUE),
                   pass_attempts = sum(pass == 1, na.rm = TRUE),
                   rush_attempts = sum(rush == 1, na.rm = TRUE)),
               by = drive_keys]
drives <- merge(drives, epa_drv, by = drive_keys)

## Sequence within game ------------------------------------------------------
setorder(drives, game_id, fixed_drive)

## Build lagged outcomes for the SAME OFFENSIVE TEAM ------------------------
## prior_own_off: outcome of this offense's previous offensive drive
## prior_own_def: outcome when this team was on defense most recently before now
## prior_opp_off: outcome of opponent's most recent offensive drive

build_lags <- function(d) {
  ## Two perspectives: offense (this team has the ball) and defense (other team has ball)
  ## For each game, walk drives in order:
  d <- as.data.table(d)
  setorder(d, game_id, fixed_drive)

  ## For each drive, the team possessing the ball (posteam) had its defense on the field
  ## during the immediately preceding drive (if any, and only if posteam differs).
  d[, prev_posteam := shift(posteam, type = "lag"), by = game_id]
  d[, prev_off_success := shift(off_success, type = "lag"), by = game_id]
  d[, prev_def_success := shift(def_success, type = "lag"), by = game_id]
  d[, prev_drive_result := shift(drive_result, type = "lag"), by = game_id]
  d[, prev_drive_epa := shift(drive_epa, type = "lag"), by = game_id]
  d[, prev_yardline_100 := shift(yardline_100, type = "lag"), by = game_id]

  ## prior_own_def: did THIS team's defense just get a stop?
  ## i.e., the previous drive was the opponent on offense, and that drive's def_success
  ## from THIS team's perspective = previous drive's def_success indicator.
  d[, prior_own_def_success := ifelse(!is.na(prev_posteam) & prev_posteam != posteam,
                                       prev_def_success, NA_integer_)]
  d[, prior_own_def_disaster := ifelse(!is.na(prev_posteam) & prev_posteam != posteam,
                                        as.integer(prev_drive_result %in%
                                                   c("Touchdown","Field goal")), NA_integer_)]
  d[, prior_opp_drive_epa := ifelse(!is.na(prev_posteam) & prev_posteam != posteam,
                                     prev_drive_epa, NA_real_)]

  ## prior_own_off: previous OFFENSIVE drive by THIS team (look back further if needed)
  ## Walk through each game and for each row, find the most recent earlier drive
  ## with posteam == current posteam.
  d[, drive_idx := seq_len(.N), by = game_id]
  prior_own_off <- d[, {
    rs <- rep(NA_integer_, .N)
    re <- rep(NA_real_, .N)
    rd <- rep(NA_character_, .N)
    last_off <- list()
    for (i in seq_len(.N)) {
      tm <- posteam[i]
      if (!is.null(last_off[[tm]])) {
        rs[i] <- off_success[last_off[[tm]]]
        re[i] <- drive_epa[last_off[[tm]]]
        rd[i] <- drive_result[last_off[[tm]]]
      }
      last_off[[tm]] <- i
    }
    list(drive_idx = drive_idx,
         prior_own_off_success = rs,
         prior_own_off_epa     = re,
         prior_own_off_result  = rd)
  }, by = game_id]
  d <- merge(d, prior_own_off, by = c("game_id","drive_idx"))

  ## prior_opp_off: most recent offensive drive by the OPPONENT
  prior_opp_off <- d[, {
    rs <- rep(NA_integer_, .N)
    re <- rep(NA_real_, .N)
    last_off <- list()
    for (i in seq_len(.N)) {
      opp <- defteam[i]
      if (!is.null(last_off[[opp]])) {
        rs[i] <- off_success[last_off[[opp]]]
        re[i] <- drive_epa[last_off[[opp]]]
      }
      last_off[[posteam[i]]] <- i
    }
    list(drive_idx = drive_idx,
         prior_opp_off_success = rs,
         prior_opp_off_epa     = re)
  }, by = game_id]
  d <- merge(d, prior_opp_off, by = c("game_id","drive_idx"))

  d
}

cat("Building lag structure...\n")
drives <- build_lags(drives)

## Game-state controls --------------------------------------------------------
drives[, garbage_time := as.integer(abs(score_differential) >= 21 & qtr >= 4)]
drives[, half := ifelse(qtr <= 2, 1L, 2L)]
drives[, fp_decile := cut(yardline_100, breaks = seq(0, 100, 10),
                           include.lowest = TRUE, labels = FALSE)]
drives[, home_off := as.integer(posteam == home_team)]

## Drop drives that are end-of-half kneel-downs or end-of-game
drives <- drives[!drive_result %in% c("End of half","End of game") | is.na(drive_result)]

## Save -----------------------------------------------------------------------
saveRDS(drives, file.path(DATA_DIR, "drives.rds"))
fwrite(drives, file.path(DATA_DIR, "drives.csv"))

cat("Saved drives:", nrow(drives), "rows\n")
cat("  off_success rate:", round(mean(drives$off_success), 3), "\n")
cat("  def_success rate:", round(mean(drives$def_success), 3), "\n")
cat("  with prior_own_def_success:", sum(!is.na(drives$prior_own_def_success)), "\n")
cat("  with prior_own_off_success:", sum(!is.na(drives$prior_own_off_success)), "\n")
cat("  with prior_opp_off_success:", sum(!is.na(drives$prior_opp_off_success)), "\n")
