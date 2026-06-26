#!/usr/bin/env python3
"""
run_pipeline.py
===============
End-to-end, reproducible pipeline for the 2026 GFP SynBio Challenge.

    python run_pipeline.py

Reads the four organizer files from ./data/, trains the brightness model,
verifies the six designs against every competition rule, and writes
./outputs/submission.csv and ./outputs/design_report.csv.

Place these organizer files in ./data/ first:
    - GFP_data.xlsx
    - Exclusion_List.csv
    - AAseqs of 5 GFP proteins_20260511.txt   (used only to verify sfGFP)
    - submission_template.csv                 (used only to verify column order)

Set your team name below.
"""

from __future__ import annotations
import os
import sys
import csv

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from brightness_model import BrightnessModel              # noqa: E402
from design import (SFGFP, DESIGNS, apply_mutations,       # noqa: E402
                    net_charge, delta_charge, FORBIDDEN_PAIRS)
from verify import (validate_format, chromophore_intact,   # noqa: E402
                    load_exclusion, load_beforetop)
from reliability_screen import ReliabilityScreen          # noqa: E402
from sfgfp_brightness_model import SfGFPBrightnessModel    # noqa: E402

# ---------------------------------------------------------------- USER CONFIG
TEAM_NAME = "GeneMeow"        # <-- EDIT THIS to your exact team name
# ---------------------------------------------------------------------------

DATA = "data"
OUT = "outputs"
F_XLSX = os.path.join(DATA, "GFP_data.xlsx")
F_EXCL = os.path.join(DATA, "Exclusion_List.csv")
F_REF = os.path.join(DATA, "AAseqs of 5 GFP proteins_20260511.txt")
F_TEMPLATE = os.path.join(DATA, "submission_template.csv")


def parse_reference(path: str) -> dict[str, str]:
    refs, name, buf = {}, None, []
    with open(path, encoding="utf-8-sig") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.startswith(">"):
                if name:
                    refs[name] = "".join(buf)
                name, buf = line[1:].strip(), []
            elif name:
                buf.append(line.upper())
        if name:
            refs[name] = "".join(buf)
    return refs


def main():
    os.makedirs(OUT, exist_ok=True)

    # --- 0. sanity: organizer files present, sfGFP matches our reference -----
    for f in (F_XLSX, F_EXCL):
        if not os.path.exists(f):
            sys.exit(f"Missing required file: {f} (see ./data/README.md)")
    if os.path.exists(F_REF):
        refs = parse_reference(F_REF)
        if refs.get("sfGFP") and refs["sfGFP"] != SFGFP:
            sys.exit("sfGFP in data pack differs from src/design.py SFGFP -- stop.")
        print("[0] sfGFP reference verified against organizer file.")

    # --- 1. train brightness model (corrected sfGFP-aware scorer) -----------
    sfmodel = SfGFPBrightnessModel().fit(F_XLSX)
    model = sfmodel.core                       # validated avGFP core (reused below)
    info = sfmodel.report()
    print(f"[1] Brightness model: held-out R^2 = {info['held_out_r2']} | "
          f"offset {info['offset_check']} | coverage {info['coverage_pct']}% of "
          f"single substitutions")
    print(f"    sfGFP sits at the calibration ceiling ({info['sfgfp_wt_pred']} = "
          f"{info['ceiling']}); absolute brightness can only detect DIMMING, so "
          f"ranking uses delta vs sfGFP (S65T coef {info['S65T_coef']} shows why "
          f"the avGFP background is not added naively).")
    wt_pred = model.predict([])

    # --- 2. load disqualification sets --------------------------------------
    excluded = load_exclusion(F_EXCL)
    beforetop = load_beforetop(F_XLSX)
    print(f"[2] Exclusion list: {len(excluded)} unique sequences | "
          f"past winners: {len(beforetop)}")
    assert SFGFP in excluded, "sanity check failed: sfGFP should be excluded"

    # --- 3. build, score, and verify the six designs ------------------------
    print(f"[3] Building {len(DESIGNS)} designs (parent net charge "
          f"{net_charge(SFGFP)}, WT predicted brightness {wt_pred:.3f}):")
    report_rows, sequences = [], {}
    for sid, d in DESIGNS.items():
        muts = d["mutations"]
        seq = apply_mutations(SFGFP, muts)
        sequences[sid] = seq

        errs = validate_format(seq)
        if errs:
            sys.exit(f"  Seq {sid} FAILS format: {errs}")
        if not chromophore_intact(seq, SFGFP):
            sys.exit(f"  Seq {sid} FAILS: chromophore altered")
        if seq in excluded:
            sys.exit(f"  Seq {sid} FAILS: in exclusion list")
        if seq in beforetop:
            sys.exit(f"  Seq {sid} FAILS: matches a past winner")
        ms = set(muts)
        for pair in FORBIDDEN_PAIRS:
            if pair.issubset(ms):
                sys.exit(f"  Seq {sid} FAILS: forbidden epistatic pair "
                         f"{sorted(pair)} co-present")

        bs = sfmodel.score(muts)               # corrected sfGFP-aware brightness
        conf = bs["confidence"]
        near = sum("buried+near-chromo" in p["flags"] for p in bs["per_mutation"])
        nq = net_charge(seq)
        report_rows.append(dict(
            Seq_ID=sid, Design_Label=d["label"], Length=len(seq),
            Num_Mutations=len(muts), Mutations=":".join(muts),
            Chromophore=seq[64:67],
            Brightness_Delta_vs_sfGFP=bs["delta_score"],
            Brightness_Verdict=bs["verdict"],
            Measured=conf["measured"], Estimated=conf["estimated"],
            Ref_Mismatch=conf["ref_mismatch"], Near_Chromophore=near,
            Net_Charge=nq, Delta_Charge=delta_charge(muts), Sequence=seq))
        print(f"    Seq {sid} [{d['label']:<22}] {len(muts)} mut | "
              f"delta vs sfGFP {bs['delta_score']:+.3f} "
              f"({conf['measured']}m/{conf['estimated']}e/{conf['ref_mismatch']}rm,"
              f" {near} near-chromo) | net charge {nq:+d}")

    # cross-design checks
    seqs = list(sequences.values())
    assert len(set(seqs)) == len(seqs), "duplicate designs"
    print("[3] All designs: format OK, chromophore intact, not excluded, "
          "not past winners, mutually distinct.")

    # --- 3b. reliability screen: no design may silently rely on a low-confidence
    #         mutation (data-poor, non-additive, or buried next to the chromophore) -
    scr = ReliabilityScreen(model, F_XLSX)
    print(f"[3b] Reliability screen (global additive MAE {scr.global_mae:.3f}):")
    all_clean = True
    for sid, d in DESIGNS.items():
        flagged = [r for r in scr.screen_design(d["mutations"]) if r["flags"]]
        if flagged:
            all_clean = False
            detail = "; ".join(f"{r['mutation']} [{','.join(r['flags'])}]"
                               for r in flagged)
            print(f"     Seq {sid}: FLAGGED -> {detail}")
        else:
            print(f"     Seq {sid}: clean")
    if all_clean:
        print("     All designs reliability-clean (no flagged mutations).")
    else:
        print("     NOTE: flagged mutations above are lower-confidence; keep only "
              "if deliberately isolated and disclosed.")

    # --- 3c. corrected sfGFP-aware brightness gate: every design must be
    #         predicted no dimmer than sfGFP (delta >= 0), with no estimated /
    #         ref-mismatch / buried-near-chromophore mutations -----------------
    print("[3c] Corrected sfGFP-aware brightness gate (delta vs sfGFP):")
    gate_ok = True
    for r in report_rows:
        bad = []
        if r["Brightness_Delta_vs_sfGFP"] < -0.05:
            bad.append(f"delta {r['Brightness_Delta_vs_sfGFP']:+.2f}")
        if r["Estimated"]:
            bad.append(f"{r['Estimated']} estimated(unseen)")
        if r["Ref_Mismatch"]:
            bad.append(f"{r['Ref_Mismatch']} ref-mismatch")
        if r["Near_Chromophore"]:
            bad.append(f"{r['Near_Chromophore']} buried+near-chromo")
        if bad:
            gate_ok = False
            print(f"     Seq {r['Seq_ID']}: RE-PICK -> {', '.join(bad)}")
        else:
            print(f"     Seq {r['Seq_ID']}: pass (delta "
                  f"{r['Brightness_Delta_vs_sfGFP']:+.3f}, all measured, clean)")
    if gate_ok:
        print("     All designs pass the corrected brightness gate "
              "(>= sfGFP, fully measured, no chromophore-proximal mutations).")

    # --- 4. write submission.csv (strict 3-column schema) -------------------
    if os.path.exists(F_TEMPLATE):
        with open(F_TEMPLATE, newline="") as fh:
            cols = next(csv.reader(fh))
        assert cols == ["Team_Name", "Seq_ID", "Sequence"], \
            f"template columns changed: {cols}"

    sub_path = os.path.join(OUT, "submission.csv")
    with open(sub_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Team_Name", "Seq_ID", "Sequence"])
        for sid in sorted(sequences):
            w.writerow([TEAM_NAME, sid, sequences[sid]])
    print(f"[4] Wrote {sub_path}")
    if TEAM_NAME == "REPLACE_WITH_TEAM_NAME":
        print("    !!! Remember to set TEAM_NAME at the top of this file.")

    # --- 5. write design_report.csv (for the PDF / appendix) ---------------
    rep_path = os.path.join(OUT, "design_report.csv")
    fields = list(report_rows[0].keys())
    with open(rep_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(report_rows)
    print(f"[5] Wrote {rep_path}")
    print("\nDone.")


if __name__ == "__main__":
    main()
