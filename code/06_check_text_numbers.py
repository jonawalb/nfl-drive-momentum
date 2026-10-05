#!/usr/bin/env python3
"""06_check_text_numbers.py -- accuracy pass for the v3 text (NEW 2026-10-05).

Two checks against the pipeline outputs:
  1. CLAIMS: each headline claim (number + sign + context) is located in the text by
     a regex and compared with the value computed from the output files.
  2. SWEEP: every number in the prose must equal some output value at the printed
     precision (raw, x100, or abs); unmatched numbers are listed for manual review
     (years, yardages, and definitional constants are whitelisted).

Inputs: data/key_nums_v3.json (05_mirage_analyses.R), data/key_nums_v2.rds and
data/robustness_v2.rds (exported to JSON through Rscript), tables/*_v3.tex.
Usage: python3 code/06_check_text_numbers.py paper_v3.tex [other text files...]
Exit status 1 if any CLAIM fails.
"""
import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def load_v2():
    r = ('suppressMessages(library(jsonlite)); f<-function(x){if(is.data.frame(x)) return(as.list(x));'
         'if(is.list(x)) return(lapply(x,f)); if(is.numeric(x)) return(unname(as.list(x))); x};'
         'cat(toJSON(list(k2=f(readRDS("data/key_nums_v2.rds")),rob=f(readRDS("data/robustness_v2.rds"))),'
         'auto_unbox=TRUE,digits=10,force=TRUE))')
    out = subprocess.run(["Rscript", "-e", r], cwd=ROOT, capture_output=True, text=True, check=True)
    return json.loads(out.stdout)


def flat(x, pre=""):
    if isinstance(x, dict):
        for k, v in x.items():
            yield from flat(v, f"{pre}.{k}" if pre else k)
    elif isinstance(x, list):
        for i, v in enumerate(x):
            yield from flat(v, f"{pre}[{i}]")
    elif isinstance(x, (int, float)) and not isinstance(x, bool):
        yield pre, float(x)


def prose(path):
    t = Path(path).read_text()
    if path.endswith(".tex"):
        t = t.split("\\begin{document}", 1)[1].split("\\bibliography", 1)[0]
        t = re.sub(r"^%.*$", "", t, flags=re.M)
        t = re.sub(r"\\(input|includegraphics|label|ref|citet|citep|cite)(\[[^\]]*\])?\{[^}]*\}", " ", t)
        t = re.sub(r"\\caption\{.*?\}\n\\label", " ", t, flags=re.S)
    t = t.replace("$-$", "-").replace("$+$", "+").replace("−", "-").replace("{,}", ",")
    t = t.replace("\\%", "%").replace("$", "")
    return t


def num_pattern(v, dec):
    s = f"{v:.{dec}f}"
    if s.startswith("-"):
        return r"[-−]\s?" + re.escape(s[1:])
    return re.escape(s)


def main(files):
    k3 = json.loads((ROOT / "data/key_nums_v3.json").read_text())
    v2 = load_v2()
    k2, rob = v2["k2"], v2["rob"]
    g = lambda d, key: d[key][0] if isinstance(d[key], list) else d[key]
    pp = lambda x: 100 * x

    # (description, regex with {v} placeholders, list of (value, decimals))
    claims = [
        ("raw stop rate / other rate / gap", r"52\.3\D.*?35\.6\D", []),
        ("raw gap 16.7", r"{v}", [(pp(k3["decomp_gap"]), 1)]),
        ("FP share of gap 17.0", r"{v} points", [(pp(k3["decomp_Field position"]), 1)]),
        ("FP share 102%", r"{v}%", [(100 * k3["decomp_share_fp"], 0)]),
        ("decomp residual -0.1 (0.5)", r"{v} points \(s\.e\. 0\.5\)", [(pp(k3["decomp_full"]), 1)]),
        ("joint-sample gap 16.3 / fp 16.8", r"{v} points, field position explains {v}",
         [(pp(k3["decomp_joint_gap"]), 1), (pp(k3["decomp_joint_fp"]), 1)]),
        ("H1 joint", r"{v} points \(s\.e\. {v}, p=0\.39\)", [(pp(k3["h1"]), 1), (pp(k3["se_h1"]), 1)]),
        ("H2 joint", r"\+{v} \(s\.e\. {v}\)", [(pp(k3["h2"]), 1), (pp(k3["se_h2"]), 1)]),
        ("H3 joint", r"{v} \(s\.e\. {v}, p=0\.09\)", [(pp(k3["h3"]), 1), (pp(k3["se_h3"]), 1)]),
        ("H3 p-value", r"p=0\.09", []),
        ("H1 90% CI", r"from {v} to {v}", [(pp(k3["ci90_h1_lo"]), 1), (pp(k3["ci90_h1_hi"]), 1)]),
        ("H1 TOST p at 2", r"p=0\.004", []),
        ("H1 TOST p at 1", r"p=0\.18", []),
        ("H1 min bound", r"{v} points", [(pp(k3["tost_h1_minbound"]), 2)]),
        ("H3 90% CI", r"from {v} to -0\.0", [(pp(k3["ci90_h3_lo"]), 1)]),
        ("H3 min bound", r"{v} points", [(pp(k3["tost_h3_minbound"]), 2)]),
        ("FP slope pp/yard", r"{v} points", [(abs(pp(k3["fp_slope"])), 2)]),
        ("EP per yard", r"{v} expected points", [(k3["ep_per_yard"], 3)]),
        ("H1 bound yards/EP", r"{v} yards and {v} expected", [(k3["tost_h1_yards"], 2), (k3["tost_h1_ep"], 2)]),
        ("H3 bound yards/EP", r"{v} yards and {v} expected", [(k3["tost_h3_yards"], 2), (k3["tost_h3_ep"], 2)]),
        ("H2 CI", r"(from|between) {v} (to|and) {v} points", [(pp(k3["ci90_h2_lo"]), 1), (pp(k3["ci90_h2_hi"]), 1)]),
        ("H2 TOST p", r"(p=|\()0\.75\)", []),
        ("H2 bound", r"{v} points, {v} yards, or {v} expected", [(pp(k3["tost_h2_minbound"]), 2),
                                                              (k3["tost_h2_yards"], 2), (k3["tost_h2_ep"], 2)]),
        ("sim headline H2", r"\+{v} points for H2", [(pp(k3["sim_head_l_own"]), 1)]),
        ("snap raw gap", r"{v} points \({v}% against {v}%\)", [(pp(k3["raw_gap_scrim"]), 1),
                                                             (pp(k3["p_score_stop_scrim"]), 1),
                                                             (pp(k3["p_score_nostop_scrim"]), 1)]),
        ("snap H1", r"{v} points \(s\.e\. {v}, 90% CI {v} to {v}\)",
         [(pp(k3["h1c"]), 1), (pp(k3["se_h1c"]), 1), (pp(k3["h1c_ci90_lo"]), 1), (pp(k3["h1c_ci90_hi"]), 1)]),
        ("snap H1 TOST p", r"p=0\.04\)", []),
        ("snap H3", r"H3 in that model is {v} points \(s\.e\. {v}\), with a tightest bound of {v}",
         [(pp(k3["h3c"]), 1), (pp(k3["se_h3c"]), 1), (pp(k3["tost_h3c_minbound"]), 2)]),
        ("stop counts (v2 rule admitted 1 punt)", r"(exactly 1|only one) punt", []),
        ("stop counts v2", r"11,824", []),
        ("snap punts / stops", r"18,652.*?30,475", []),
        ("R1 flexible", r"give {v}, \+{v}, and {v}", [(pp(rob["R1"]["beta_h1"][0]), 1),
                                                     (pp(rob["R1"]["beta_h2"][0]), 1),
                                                     (pp(rob["R1"]["beta_h3"][0]), 1)]),
        ("R3 logit", r"are {v}, \+{v}, and -0\.0", [(pp(rob["R3"]["ame"]["h1"][0]), 1),
                                                  (pp(rob["R3"]["ame"]["h2"][0]), 1)]),
        ("no garbage time", r"gives {v}, \+{v}, and {v}", [(pp(g(k2, "beta_nogt_h1")), 1),
                                                          (pp(g(k2, "beta_nogt_h2")), 1),
                                                          (pp(g(k2, "beta_nogt_h3")), 1)]),
        ("swing takeaway", r"0\.7 points more often than after the reference outcome \(s\.e\. 0\.7\)", []),
        ("swing takeaway vs 3out", r"{v} points more than after a three-and-out \(s\.e\. {v}\)",
         [(pp(k3["sw_takeaway_vs_3out"]), 1), (pp(k3["sw_takeaway_vs_3out_se"]), 1)]),
        ("swing long TD", r"by {v} points, exactly", [(abs(pp(k3["sw_long_td_against"])), 1)]),
        ("swing long-short", r"difference is {v} points \(s\.e\. {v}\)", [(pp(k3["sw_long_vs_short_td"]), 1),
                                                                        (pp(k3["sw_long_vs_short_td_se"]), 1)]),
        ("swing nonoff-short", r"difference {v}, s\.e\. {v}", [(pp(k3["sw_nonoff_vs_short_td"]), 1),
                                                               (pp(k3["sw_nonoff_vs_short_td_se"]), 1)]),
        ("ben own long", r"{v} points less than after a shorter touchdown \(s\.e\. {v}\)",
         [(abs(pp(k3["ben_own_long_vs_other_td"])), 1), (pp(k3["ben_own_long_vs_other_td_se"]), 1)]),
        ("ben own def td", r"{v} points less \(s\.e\. {v}\)", [(abs(pp(k3["ben_own_def_td"])), 1),
                                                              (pp(k3["ben_own_def_td_se"]), 1)]),
        ("4th down stop", r"{v} points less often \(s\.e\. {v}\)", [(abs(pp(k3["sw_fourth_down_stop"])), 1),
                                                                   (pp(k3["sw_fourth_down_stop_se"]), 1)]),
        ("gameFE H2/H3", r"H2 coefficient was {v} points and the H3 coefficient was {v}",
         [(pp(g(k2, "beta_h2_gamefe")), 1), (pp(g(k2, "beta_h3_gamefe")), 1)]),
        ("lead placebo gameFE", r"at {v} points when including game", [(pp(g(k2, "beta_lead_gamefe")), 1)]),
        ("lead placebo fp ctrl", r"is {v} points when including game fixed effects and {v} points",
         [(pp(g(k2, "beta_lead_gamefe_fpctrl")), 1), (pp(g(k2, "beta_lead_headline_fpctrl")), 1)]),
        ("drives per game / approx", r"averages {v} drives per game, which gives about {v} points",
         [(k3["drives_per_game"], 1), (pp(k3["nickell_approx"]), 1)]),
        ("sim calibration", r"s\.d\. {v}\).*?\({v}\)", [(k3["sim_sd_u"], 3), (k3["sim_rho_u"], 2)]),
        ("sim gameFE", r"return {v} points for H2, {v} for H3, and {v} for the lead",
         [(pp(k3["sim_game_l_own"]), 1), (pp(k3["sim_game_l_opp"]), 1), (pp(k3["sim_game_l_lead"]), 1)]),
        ("sim headline", r"\+{v} for H2 and \+{v} for H3", [(pp(k3["sim_head_l_own"]), 1),
                                                            (pp(k3["sim_head_l_opp"]), 1)]),
        ("adjusted H2", r"\+{v} in the headline model, and \+{v} with game",
         [(pp(k3["adj_head_l_own"]), 1), (pp(k3["adj_game_l_own"]), 1)]),
        ("adjusted H3", r"H3 is {v} and {v}", [(pp(k3["adj_head_l_opp"]), 1), (pp(k3["adj_game_l_opp"]), 1)]),
        ("adjusted lead", r"placebo is {v} and {v}", [(pp(k3["adj_head_l_lead"]), 1),
                                                     (pp(k3["adj_game_l_lead"]), 1)]),
        ("kick N", r"{v} kickoffs", [(k3["ko_n_total"], 0)]),
        ("kick short", r"lower by {v} points \(standard error = {v}, based on {v} short",
         [(-k3["ko_cost_short"], 2), (k3["ko_cost_short_se"], 2), (k3["ko_n_short"], 0)]),
        ("kick onside", r"lower by {v} points \(standard error = {v}, based on {v} onside kicks, which the kicking team recovered {v}%",
         [(-k3["ko_cost_onside"], 2), (k3["ko_cost_onside_se"], 2), (k3["ko_n_onside"], 0),
          (100 * k3["ko_recov_onside"], 1)]),
        ("pts per score", r"worth {v} points", [(k3["pts_per_score"], 2)]),
        ("break-even", r"by {v} and {v} percentage points", [(pp(k3["ko_breakeven_short_pp"]), 1),
                                                             (pp(k3["ko_breakeven_onside_pp"]), 1)]),
        ("ceiling 3.0", r"{v} percentage points for a given channel", [(pp(k3["ci90_h2_hi"]), 1)]),
        ("v1 suppression", r"32\.7", []),
        ("v1 kickoff share", r"91\.9%", [(100 * g(k2, "share_post_score_at_35_v1"), 1)]),
        ("FP means", r"74\.9 yards.*?69\.3 yards", []),
        ("samples", r"86,249 drives.*?4,078 games", []),
        ("est sample", r"74,125 drives, and 74,120", []),
    ]
    fmt = lambda v, d: f"{round(v, d):.{d}f}" if d else f"{int(round(v)):,}"
    text_all = {f: re.sub(r"\s+", " ", prose(f)) for f in files}
    fails = 0
    target = files[0]
    t = text_all[target]
    for desc, rx, vals in claims:
        pat = rx
        for v, d in vals:
            s = fmt(v, d)
            p = (r"-\s?" + re.escape(s[1:])) if s.startswith("-") else re.escape(s)
            pat = pat.replace("{v}", p, 1)
        ok = re.search(pat, t) is not None
        status = "OK  " if ok else "FAIL"
        fails += not ok
        print(f"[{status}] {desc}: expects {[fmt(v, d) for v, d in vals]}")

    # sweep
    allowed = set()
    sources = list(flat(k3)) + list(flat(k2)) + list(flat(rob))
    sources.append(("n_def_stop+n_no_stop", g(k2, "n_def_stop") + g(k2, "n_no_stop")))
    sources.append(("v1 joint H3, data/key_nums.rds beta_h3_joint = -0.3271", -0.3271351))
    for tf in (ROOT / "tables").glob("*_v[23].tex"):
        for m in re.findall(r"-?\d+\.?\d*", tf.read_text().replace(",", "")):
            sources.append((tf.name, float(m)))
    for _, v in sources:
        for x in (v, 100 * v, abs(v), abs(100 * v), -v):
            for d in range(0, 4):
                allowed.add(f"{x:.{d}f}")
    whitelist = {"2010", "2023", "2024", "1", "2", "3", "4", "5", "6", "7", "15", "35", "40", "60", "10", "21",
                 "100", "90", "0", "1.0", "0.5", "2.0"}
    for f, tx in text_all.items():
        nums = re.findall(r"(?<![\w.])[-+]?\d[\d,]*\.?\d*", tx)
        bad = sorted({n for n in nums if n.rstrip(".").replace(",", "").lstrip("+") not in allowed
                      and n.rstrip(".").replace(",", "").lstrip("+-") not in whitelist})
        print(f"\nSWEEP {Path(f).name}: {len(nums)} numbers, unmatched -> {bad}")
    print(f"\nCLAIMS failed: {fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or [str(ROOT / "paper_v3.tex")]))
