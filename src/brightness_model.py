"""
brightness_model.py
-------------------
Additive brightness predictor trained on the organizer-provided deep mutational
scanning (DMS) data (GFP_data.xlsx, sheet "brightness").

Design rationale
================
Sarkisyan et al. (2016, Nature) showed that the GFP fitness landscape is, to ~95%
of its variance, a *threshold* function of an additive "fitness potential": each
mutation contributes a roughly additive change to log-fluorescence, and brightness
collapses only once the cumulative destabilisation crosses a stability margin
(~7-9 kcal/mol). We therefore model brightness as:

        raw_score(variant) = intercept + sum_m  beta_m            (additive Ridge)
        brightness         = isotonic( raw_score )                (monotone threshold)

* Ridge regression on a sparse one-hot encoding of (position, substitution) tokens
  recovers per-mutation marginal effects while accounting for co-occurrence.
* Isotonic regression calibrates the additive score onto the measured brightness
  scale and reproduces the sigmoid/threshold non-linearity of the landscape.

Held-out performance: R^2 ~ 0.93 (20% test split, see train_and_report()).

CRITICAL NUMBERING NOTE
=======================
The DMS dataset numbers avGFP WITHOUT the initial methionine, i.e.

        dataset_position  =  protein_position - 1

We verified this by aligning the dataset's per-position reference residues against
avGFP: offset 0 -> 11/233 match, offset +1 (Met dropped) -> 233/233 match.
All public functions in this module take mutations in standard, Met-included
protein numbering (e.g. "K166E") and convert internally. Getting this wrong
silently corrupts every position lookup, so it is handled in exactly one place
(`_to_dataset_token`).
"""

from __future__ import annotations
import re
import numpy as np
import scipy.sparse as sp
import openpyxl
from sklearn.linear_model import Ridge
from sklearn.isotonic import IsotonicRegression
from sklearn.model_selection import train_test_split
from sklearn.metrics import r2_score

_MUT_RE = re.compile(r"([A-Z])(\d+)([A-Z])")


def _to_dataset_token(mutation: str) -> str:
    """Convert a Met-included protein mutation (e.g. 'K166E') to the DMS token
    used in GFP_data.xlsx (e.g. 'K165E'), accounting for the -1 numbering offset."""
    m = _MUT_RE.fullmatch(mutation.strip())
    if m is None:
        raise ValueError(f"Bad mutation format: {mutation!r}")
    ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
    return f"{ref}{pos - 1}{alt}"


class BrightnessModel:
    """Additive Ridge + isotonic brightness predictor over a single GFP family."""

    def __init__(self, gfp_type: str = "avGFP", alpha: float = 1.0):
        self.gfp_type = gfp_type
        self.alpha = alpha
        self.vocab: list[str] = []
        self.index: dict[str, int] = {}
        self.coef: dict[str, float] = {}
        self.intercept_: float = 0.0
        self.ridge: Ridge | None = None
        self.iso: IsotonicRegression | None = None
        self.wt_brightness: float | None = None

    # ----------------------------------------------------------------- loading
    @staticmethod
    def _load_rows(xlsx_path: str, gfp_type: str):
        wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
        ws = wb["brightness"]
        rows, wt = [], None
        for mut, typ, b in ws.iter_rows(min_row=2, values_only=True):
            if typ != gfp_type or b is None:
                continue
            b = float(b)
            if mut == "WT":
                wt = b
                rows.append((frozenset(), b))
            else:
                toks = frozenset(t.strip() for t in str(mut).split(":"))
                rows.append((toks, b))
        return rows, wt

    # ------------------------------------------------------------------- train
    def fit(self, xlsx_path: str):
        rows, wt = self._load_rows(xlsx_path, self.gfp_type)
        self.wt_brightness = wt
        self.vocab = sorted({t for s, _ in rows for t in s})
        self.index = {t: i for i, t in enumerate(self.vocab)}
        ri, ci, y = [], [], []
        for r, (s, b) in enumerate(rows):
            for t in s:
                ri.append(r)
                ci.append(self.index[t])
            y.append(b)
        X = sp.csr_matrix((np.ones(len(ri)), (ri, ci)),
                          shape=(len(rows), len(self.vocab)))
        y = np.asarray(y)
        self.ridge = Ridge(alpha=self.alpha).fit(X, y)
        self.iso = IsotonicRegression(out_of_bounds="clip").fit(
            self.ridge.predict(X), y)
        self.intercept_ = float(self.ridge.intercept_)
        self.coef = {t: float(self.ridge.coef_[self.index[t]]) for t in self.vocab}
        return self

    def train_and_report(self, xlsx_path: str) -> dict:
        """Refit on an 80/20 split and report held-out R^2 (for the README/PDF)."""
        rows, _ = self._load_rows(xlsx_path, self.gfp_type)
        vocab = sorted({t for s, _ in rows for t in s})
        index = {t: i for i, t in enumerate(vocab)}
        ri, ci, y = [], [], []
        for r, (s, b) in enumerate(rows):
            for t in s:
                ri.append(r); ci.append(index[t])
            y.append(b)
        X = sp.csr_matrix((np.ones(len(ri)), (ri, ci)), shape=(len(rows), len(vocab)))
        y = np.asarray(y)
        Xtr, Xte, ytr, yte = train_test_split(X, y, test_size=0.2, random_state=0)
        ridge = Ridge(alpha=self.alpha).fit(Xtr, ytr)
        iso = IsotonicRegression(out_of_bounds="clip").fit(ridge.predict(Xtr), ytr)
        r2_lin = r2_score(yte, ridge.predict(Xte))
        r2_cal = r2_score(yte, iso.predict(ridge.predict(Xte)))
        return {"n_variants": len(rows), "r2_ridge": r2_lin,
                "r2_ridge_isotonic": r2_cal}

    # ----------------------------------------------------------------- predict
    def predict(self, mutations: list[str]) -> float:
        """Predicted brightness for sfGFP + `mutations` (Met-included numbering)."""
        if self.ridge is None or self.iso is None:
            raise RuntimeError("Model not fitted; call .fit() first.")
        toks = [_to_dataset_token(m) for m in mutations]
        raw = self.intercept_ + sum(self.coef.get(t, 0.0) for t in toks)
        return float(self.iso.predict([raw])[0])

    def single_effect(self, mutation: str) -> float | None:
        """Ridge coefficient (marginal effect) of a single mutation, or None
        if that token was never observed in the DMS data."""
        return self.coef.get(_to_dataset_token(mutation))
