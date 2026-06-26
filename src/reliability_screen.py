#!/usr/bin/env python3
"""
reliability_screen.py
=====================
Systematic per-mutation reliability screen for the additive brightness model.

The additive model assumes each mutation's effect is independent. That assumption
is least trustworthy where a mutation is (a) poorly supported by the DMS, (b)
demonstrably non-additive in the multi-mutant data, or (c) buried in the core right
next to the chromophore (packing + direct chromophore effects that do not compose
additively, and where the avGFP->sfGFP transfer is most fragile).

This module flags any candidate mutation on four independent signals so that the
pipeline never silently relies on a low-confidence mutation (the failure mode that
otherwise gets caught only by eye).

Signals
-------
data-driven (direct):
  * support        - how many DMS variants contain the (offset-corrected) token
  * additivity     - bias and mean-abs-residual of the additive model over
                     multi-mutant variants containing the token, vs the global
                     baseline. Directly measures whether the model works for it.
structural (proxy, available even when data are sparse):
  * burial         - relative SASA on the sfGFP model
  * chromophore    - distance to the chromophore (residue 66)

Flag rules (a mutation is flagged = lower confidence, handle with care):
  UNSEEN              coefficient never observed in the DMS
  low-support        0 < support < MIN_SUPPORT
  non-additive       >=5 multi-mutant obs and (MAE > MAE_FACTOR*global or |bias|>BIAS)
  buried+near-chromo rSASA < RSASA_BURIED AND distance < CHROMO_A
                     (both required: a surface residue a little inside the cutoff,
                      e.g. K166E at 8 A, is not flagged; a buried one, e.g. L44M, is)

A flag means "the brightness prediction here is less reliable", not "forbidden".
Structurally flagged but evidence-backed mutations may still be used deliberately
if isolated and disclosed; chromophore residues (e.g. S65T) are genuinely off-limits.
"""
from __future__ import annotations
import csv
import os
import numpy as np
from brightness_model import BrightnessModel, _to_dataset_token

MIN_SUPPORT = 15
MAE_FACTOR = 1.5
BIAS = 0.30
RSASA_BURIED = 0.15
CHROMO_A = 8.0

_HERE = os.path.dirname(os.path.abspath(__file__))
_FEATURES = os.path.join(_HERE, "..", "structural_validation", "sfgfp_site_features.csv")


class ReliabilityScreen:
    def __init__(self, model: BrightnessModel, xlsx_path: str,
                 features_csv: str = _FEATURES):
        self.model = model
        rows, _ = BrightnessModel._load_rows(xlsx_path, model.gfp_type)
        raw = np.array([model.intercept_ + sum(model.coef.get(t, 0.0) for t in s)
                        for s, _ in rows])
        pred = model.iso.predict(raw)
        meas = np.array([b for _, b in rows])
        nmut = np.array([len(s) for s, _ in rows])
        self.global_mae = float(np.abs(pred - meas)[nmut >= 2].mean())
        self.support, self._rsum, self._rabs, self._rc = {}, {}, {}, {}
        for i, (s, _) in enumerate(rows):
            for t in s:
                self.support[t] = self.support.get(t, 0) + 1
                if nmut[i] >= 2:
                    r = pred[i] - meas[i]
                    self._rsum[t] = self._rsum.get(t, 0.0) + r
                    self._rabs[t] = self._rabs.get(t, 0.0) + abs(r)
                    self._rc[t] = self._rc.get(t, 0) + 1
        self.feat = {int(r["position"]): (float(r["rSASA"]),
                     float(r["dist_to_chromophore_A"]))
                     for r in csv.DictReader(open(features_csv))}

    def evaluate(self, mut: str) -> dict:
        tok = _to_dataset_token(mut)
        pos = int(mut[1:-1])
        rs, d = self.feat.get(pos, (0.5, 99.0))
        coef = self.model.coef.get(tok)
        n = self.support.get(tok, 0)
        nm = self._rc.get(tok, 0)
        mae = self._rabs.get(tok, 0.0) / nm if nm else float("nan")
        bias = self._rsum.get(tok, 0.0) / nm if nm else float("nan")
        flags = []
        if coef is None:
            flags.append("UNSEEN")
        elif n < MIN_SUPPORT:
            flags.append("low-support")
        if coef is not None and nm >= 5 and (mae > MAE_FACTOR * self.global_mae
                                             or abs(bias) > BIAS):
            flags.append("non-additive")
        if rs < RSASA_BURIED and d < CHROMO_A:
            flags.append("buried+near-chromo")
        return dict(mutation=mut, coef=coef, support=n, add_mae=mae, bias=bias,
                    rSASA=rs, dist_chromo=d, flags=flags)

    def screen_design(self, mutations: list[str]) -> list[dict]:
        return [self.evaluate(m) for m in mutations]


def main():
    import sys
    sys.path.insert(0, _HERE)
    from design import DESIGNS
    xlsx = os.path.join(_HERE, "..", "data", "GFP_data.xlsx")
    model = BrightnessModel("avGFP").fit(xlsx)
    scr = ReliabilityScreen(model, xlsx)
    print(f"global additive-model MAE = {scr.global_mae:.3f}  "
          f"(non-additive flag if mut-MAE > {MAE_FACTOR}x = {MAE_FACTOR*scr.global_mae:.3f})\n")
    all_clean = True
    for sid, d in DESIGNS.items():
        flagged = [r for r in scr.screen_design(d["mutations"]) if r["flags"]]
        tag = "CLEAN" if not flagged else "; ".join(
            f"{r['mutation']} [{','.join(r['flags'])}]" for r in flagged)
        if flagged:
            all_clean = False
        print(f"#{sid} {d['label']:24s} {tag}")
    print(f"\nall designs reliability-clean: {all_clean}")


if __name__ == "__main__":
    main()
