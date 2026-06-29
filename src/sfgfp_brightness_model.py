#!/usr/bin/env python3
"""
sfgfp_brightness_model.py
=========================
Final, corrected brightness model for the 2026 GFP challenge.

Predicts the *initial* fluorescence brightness (F_initial) of an sfGFP variant,
trained on the avGFP deep-mutational-scanning data (Sarkisyan et al., 2016).
This module is built to avoid the specific mistakes seen across other pipelines.

WHAT IT GETS RIGHT
------------------
1. Numbering. The avGFP DMS numbers the protein WITHOUT the initial methionine
   (dataset_position = protein_position - 1). This is verified by alignment at
   fit time (the reconstructed avGFP reference must match the bundled avGFP at
   that offset, ~222/233 vs ~10/233 at offset 0) and asserted, so a silent
   off-by-one (which flips the sign of several mutations) cannot slip through.

2. avGFP -> sfGFP transfer done as a DELTA, not an absolute baseline.
   sfGFP differs from avGFP at 11 positions, and several of those defining
   mutations are epistatic with misleading additive coefficients - most starkly
   S65T, whose avGFP-DMS coefficient is NEGATIVE (-0.349) even though sfGFP is in
   reality brighter than avGFP. We therefore do NOT read absolute brightness off
   the avGFP intercept plus sfGFP's background. Instead we score the *change* a
   design makes relative to its parent (sum of the new mutations' effects), which
   is the quantity the competition score actually depends on. The net of the 11
   background coefficients (+0.62) puts sfGFP at the model's calibration ceiling,
   which is the honest reason no design can be predicted *brighter* than sfGFP -
   only dimmer. So the absolute scale is used to detect/quantify DIMMING; the raw
   delta is used for ranking.

3. Reference-residue check. A mutation is scored directly only if the residue it
   mutates *from* matches both the parent (sfGFP) and the avGFP DMS reference at
   that position. Mutations at positions where sfGFP differs from avGFP (or whose
   'from' residue disagrees) are flagged - their avGFP coefficient is for the
   wrong starting residue. (Using such a coefficient here would be a silent error.)

4. Coverage-aware, never silently zero. Only ~40% of single substitutions were
   measured in the DMS. Unseen mutations are not set to 0; they receive a
   position-level fallback (mean measured effect at that residue position) and
   are flagged 'estimated'. The model thus reports what it knows vs guesses.

5. Honest confidence on every call: counts of measured / estimated /
   reference-mismatched mutations, plus structural flags (burial, chromophore
   proximity) where the additive assumption is weakest.

WHAT IT DOES NOT DO
-------------------
It predicts initial brightness only. It says nothing about fluorescence retention
after the 72 C heat challenge - the axis that actually decides the competition and
for which no validated predictor exists here. Use net charge / supercharging
(anti-aggregation) and an explicit stability tool for that, separately.

Training data: avGFP only (51,715 variants). amacGFP (8,286) is a different family
and is deliberately NOT mixed in - cross-family mixing dilutes the additive signal.
"""
from __future__ import annotations
import os
import numpy as np

_HERE = os.path.dirname(os.path.abspath(__file__))
import sys
sys.path.insert(0, _HERE)
from brightness_model import BrightnessModel, _to_dataset_token  # validated core

_REF_FASTA = os.path.join(_HERE, "..", "data", "reference_sequences.fasta")
_FEATURES = os.path.join(_HERE, "..", "structural_validation", "sfgfp_site_features.csv")

MIN_SUPPORT = 15
RSASA_BURIED = 0.15
CHROMO_A = 8.0


def _load_fasta(path):
    seqs, name, buf = {}, None, []
    for line in open(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name:
                seqs[name] = "".join(buf)
            name, buf = line[1:].strip(), []
        else:
            buf.append(line.upper())
    if name:
        seqs[name] = "".join(buf)
    return seqs


class SfGFPBrightnessModel:
    def __init__(self, ref_fasta: str = _REF_FASTA, features_csv: str = _FEATURES):
        refs = _load_fasta(ref_fasta)
        self.avgfp = refs["avGFP"]
        self.sfgfp = refs["sfGFP"]
        assert len(self.avgfp) == len(self.sfgfp), "reference length mismatch"
        self.feat = {}
        if os.path.exists(features_csv):
            import csv
            self.feat = {int(r["position"]): (float(r["rSASA"]),
                         float(r["dist_to_chromophore_A"]))
                         for r in csv.DictReader(open(features_csv))}

    # ---- fitting -----------------------------------------------------------
    def fit(self, xlsx_path: str):
        self.core = BrightnessModel(gfp_type="avGFP").fit(xlsx_path)
        self._xlsx = xlsx_path
        rows, _ = BrightnessModel._load_rows(xlsx_path, "avGFP")
        # per-token support and per-position fallback (mean measured effect)
        self.support, pos_eff = {}, {}
        for s, _ in rows:
            for t in s:
                self.support[t] = self.support.get(t, 0) + 1
        for tok, c in self.core.coef.items():
            p = int(tok[1:-1])
            pos_eff.setdefault(p, []).append(c)
        self.pos_mean = {p: float(np.mean(v)) for p, v in pos_eff.items()}
        self.global_mean = float(np.mean(list(self.core.coef.values())))
        # verify the numbering offset against the bundled avGFP reference
        self.offset_match, self.offset_total = self._verify_offset()
        frac = self.offset_match / max(self.offset_total, 1)
        assert frac > 0.8, (f"numbering offset check FAILED ({self.offset_match}/"
                            f"{self.offset_total}); refusing to run with a bad offset")
        # sfGFP background (avGFP -> sfGFP) net coefficient and ceiling
        self.background_raw = sum(
            self.core.coef.get(self._token(self.avgfp[p], p + 1, self.sfgfp[p]), 0.0)
            for p in range(len(self.avgfp)) if self.avgfp[p] != self.sfgfp[p])
        self.ceiling = float(self.core.iso.predict([self.core.intercept_ + 50.0])[0])
        self.avgfp_wt = float(self.core.iso.predict([self.core.intercept_])[0])
        self.sfgfp_wt = float(self.core.iso.predict(
            [self.core.intercept_ + self.background_raw])[0])
        return self

    def _token(self, from_aa, protein_pos, to_aa):
        return f"{from_aa}{protein_pos - 1}{to_aa}"

    def _verify_offset(self):
        ok = tot = 0
        for tok in self.core.coef:
            d = int(tok[1:-1])             # dataset index
            if 0 <= d < len(self.avgfp):   # dataset d  <->  Met-incl protein d+1
                tot += 1
                ok += (self.avgfp[d] == tok[0])
        return ok, tot

    # ---- scoring -----------------------------------------------------------
    def _effect(self, mut: str, parent: str):
        from_aa, to_aa = mut[0], mut[-1]
        pos = int(mut[1:-1])
        if not (1 <= pos <= len(parent)):
            raise ValueError(f"{mut}: position out of range")
        if parent[pos - 1] != from_aa:
            raise ValueError(f"{mut}: parent has {parent[pos-1]} at {pos}, not {from_aa}")
        flags = []
        if self.avgfp[pos - 1] != from_aa:
            flags.append("ref-mismatch(sfGFP!=avGFP here)")
        tok = self._token(from_aa, pos, to_aa)
        coef = self.core.coef.get(tok)
        measured = coef is not None
        if not measured:
            coef = self.pos_mean.get(pos - 1, self.global_mean)
            flags.append("estimated(unseen)")
        elif self.support.get(tok, 0) < MIN_SUPPORT:
            flags.append("low-support")
        rs, d = self.feat.get(pos, (0.5, 99.0))
        if rs < RSASA_BURIED and d < CHROMO_A:
            flags.append("buried+near-chromo")
        return dict(mutation=mut, effect=float(coef), measured=measured, flags=flags)

    def score(self, mutations, parent: str = "sfGFP") -> dict:
        """Score a design. delta_score>=0 => predicted at least as bright as parent.
        abs_brightness is calibrated (avGFP log scale) but saturates at the ceiling
        where sfGFP already sits; use it to gauge DIMMING, use delta_score to rank."""
        pseq = self.sfgfp if parent == "sfGFP" else self.avgfp
        praw = self.background_raw if parent == "sfGFP" else 0.0
        per = [self._effect(m, pseq) for m in mutations]
        delta = float(sum(p["effect"] for p in per))
        abs_b = float(self.core.iso.predict(
            [self.core.intercept_ + praw + delta])[0])
        n_meas = sum(p["measured"] for p in per)
        n_est = sum(not p["measured"] for p in per)
        n_ref = sum("ref-mismatch(sfGFP!=avGFP here)" in p["flags"] for p in per)
        if delta >= -1e-9:
            verdict = ">= WT (predicted no dimmer than sfGFP)"
        else:
            verdict = f"below WT (predicted dimmer by ~{-delta:.2f} in additive units)"
        return dict(
            parent=parent, delta_score=round(delta, 3), verdict=verdict,
            abs_brightness=round(abs_b, 3),
            confidence=dict(n_mutations=len(per), measured=n_meas,
                            estimated=n_est, ref_mismatch=n_ref),
            at_ceiling=abs_b >= self.ceiling - 0.01,
            per_mutation=per)

    # ---- self-test ---------------------------------------------------------
    def report(self):
        rep = self.core.train_and_report(self._xlsx)
        s65t = self.core.coef.get(self._token(self.avgfp[64], 65, self.sfgfp[64]))
        return dict(
            held_out_r2=round(rep["r2_ridge_isotonic"], 3),
            offset_check=f"{self.offset_match}/{self.offset_total} match at -1 offset",
            coverage_pct=round(100 * len(self.core.coef) / (len(self.avgfp) * 19)),
            avgfp_wt_pred=round(self.avgfp_wt, 3),
            sfgfp_wt_pred=round(self.sfgfp_wt, 3),
            ceiling=round(self.ceiling, 3),
            sfgfp_at_ceiling=self.sfgfp_wt >= self.ceiling - 0.01,
            S65T_coef=round(s65t, 3) if s65t is not None else None)


def _main():
    from design import DESIGNS
    xlsx = os.path.join(_HERE, "..", "data", "GFP_data.xlsx")
    m = SfGFPBrightnessModel().fit(xlsx)
    r = m.report()
    print("=== model report ===")
    for k, v in r.items():
        print(f"  {k}: {v}")
    print("\n=== our six designs (delta_score >= 0 means predicted >= sfGFP) ===")
    for sid, d in DESIGNS.items():
        s = m.score(d["mutations"])
        print(f"  #{sid} {d['label']:22s} delta {s['delta_score']:+.3f} | "
              f"{s['verdict']} | conf {s['confidence']} | at_ceiling={s['at_ceiling']}")
    print("\n=== sanity checks ===")
    chk = m.score(["S205Y"])
    print(f"  S205Y (chromophore proton-wire): delta {chk['delta_score']:+.3f} -> "
          f"{chk['verdict']}; flags {chk['per_mutation'][0]['flags']}")
    uns = m.score(["D190W"])  # likely unseen
    pm = uns['per_mutation'][0]
    print(f"  D190W: measured={pm['measured']} (fallback used={'estimated(unseen)' in pm['flags']}), "
          f"delta {uns['delta_score']:+.3f}")
    try:
        m.score(["A30E"])  # wrong 'from' residue for sfGFP (sfGFP has R30)
    except ValueError as e:
        print(f"  invalid mutation correctly rejected: {e}")
    rm = m.score(["R30E"])  # valid on sfGFP but sfGFP!=avGFP at 30 -> ref-mismatch
    print(f"  R30E: flags {rm['per_mutation'][0]['flags']} (sfGFP R30 vs avGFP S30)")


if __name__ == "__main__":
    _main()
