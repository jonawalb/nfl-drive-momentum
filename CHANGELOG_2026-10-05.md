# Changelog 2026-10-05 (portfolio-audit fix pass)

- `code/01_build_drives.R`: drive-start `yardline_100` now taken from the offense's first scrimmage play (kickoffs and no-play rows skipped). v1 used the drive's first row, i.e. the kickoff spot (35) for 91.9% of post-score drives. Output `data/drives_v2.{rds,csv}`.
- `code/02_analyze_momentum.R`: headline spec uses possessing-team + opponent + season FE (+ home indicator); game-FE spec kept as robustness with a lead placebo (opponent's next drive). Outputs `*_v2`.
- `code/03_robustness.R`: R1–R7 re-run on v2 data with headline FE; output `data/robustness_v2.rds`.
- `paper_v2.tex` / `.pdf`: joint H1/H2/H3 +2.0/−4.2/−32.7 → −0.5/+2.3/−0.6 pp; robustness numbers replaced; new R9 (game FE −1.0/−3.8/−6.6, lead placebo −10.4 vs −3.5; −7.7 vs −0.4 holding next start fixed); H1 now absorbed by field position; v1 Discussion 7.2–7.4 commented out with [TK]; subtitle truncated; Berger & Pope "NCAA football" → "NCAA basketball"; "in JQAS" dropped for non-JQAS works.
- README counts and `run_all.sh` updated to v2. `Sloan_SSAC27_Abstract_NFL_v2.{md,docx}` added.
- Not regenerated: `figures/interactive/` (still v1 numbers).
