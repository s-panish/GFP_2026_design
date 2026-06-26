"""
design.py
---------
Mutation pools and the six final candidate designs.

Strategy (see the PDF "Algorithm design introduction" for the full argument)
==========================================================================
The competition's combined score simplifies algebraically to

        Relative_Brightness x Thermal_Retention
            = (F_initial / F_initial_WT) x (F_final / F_initial)
            = F_final / F_initial_WT

i.e. *absolute post-heat brightness relative to WT initial brightness*. Because
GFP fluorescence is gated by folding stability (Sarkisyan 2016), both factors are
stability read-outs at two temperatures. We therefore:

1. Anchor on the sfGFP core (robust folding -> high F_initial, large stability
   margin). sfGFP itself, avGFP, and all 20 past-winner sequences are in the
   exclusion list, so we mutate sfGFP.

2. Brightness mutations: substitutions that are positive on BOTH the DMS single-
   mutant brightness AND the Ridge marginal coefficient (dual-validated). This
   guards against the negative epistasis that makes many individually-bright
   mutations dim in combination.

3. Stability mutations: negative surface supercharging (K/R/N/Q -> E/D), the
   Thermal Green Protein recipe (Close et al. 2015), restricted to positions where
   the charge change is brightness-neutral in the DMS data. Brightness-neutral in
   a stability-gated assay == folding-safe == "free" in the stability budget; the
   added negative surface charge raises Tm and suppresses irreversible aggregation
   on partial unfolding, which is what the 72 C retention score measures.

4. Six designs span a brightness<->stability frontier indexed by net charge
   (sfGFP = -6; TGP = -10; supercharged GFPs tolerate <= -30 while remaining
   fluorescent, Lawrence et al. 2007 / Der et al. 2013).

Every mutation below is in standard Met-included protein numbering.
"""

from __future__ import annotations

# sfGFP reference (Pedelacq et al. 2006); chromophore at 65-67 = TYG.
SFGFP = (
    "MSKGEELFTGVVPILVELDGDVNGHKFSVRGEGEGDATNGKLTLKFICTTGKLPVPWPTLVTTLTYG"
    "VQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGTYKTRAEVKFEGDTLVNRIELKGIDFKEDG"
    "NILGHKLEYNFNSHNVYITADKQKNGIKANFKIRHNVEDGSVQLADHYQQNTPIGDGPVLLPDNHYLS"
    "TQSVLSKDPNEKRDHMVLLEFVTAAGITHGMDELYK"
)

CHARGE = {"K": +1, "R": +1, "D": -1, "E": -1}

# --- Brightness drivers: positive Ridge coef AND bright single-mutant AND
#     reliability-clean (surface, away from chromophore, supported, additive).
#     Note: the highest-coefficient DMS mutations (D180E, K162T, ...) were rejected
#     because their coefficient far exceeds their single-mutant effect (co-occurrence
#     inflation). L44M was dropped after the structural check flagged it (buried,
#     ~7.4 A from the chromophore); Y237N replaces it as a clean surface driver. ---
BRIGHTNESS_POOL = ["D19E", "R73L", "H231Y", "Y237N"]
# Negative-supercharging, brightness-neutral/positive (thermal/aggregation resistance):
SUPERCHARGE_POOL = ["N198D", "K214E", "K166E", "N212D", "K101E", "K156E"]
# Never co-present (measured negative epistasis, -0.18 in the DMS double mutants):
FORBIDDEN_PAIRS = [frozenset({"N198D", "K214E"})]

# --- The six submitted designs (all reliability-screen-clean) ---
# Every engineered position is surface-exposed and reliability-clean. Two designs are
# dedicated single-objective champions: #5 maximises initial brightness (adds the
# screened surface driver Y237N); #4 maximises thermal retention (the DMS-validated
# supercharging maximum at net charge -15 on a bright core). #6 is an independent
# insurance design with no shared brightness drivers.
DESIGNS = {
    1: dict(label="conservative_surface",
            mutations=["D19E", "R73L", "H231Y"]),
    2: dict(label="balanced_core",
            mutations=["D19E", "R73L", "H231Y", "K166E", "N212D"]),
    3: dict(label="MAIN_balanced",
            mutations=["D19E", "R73L", "H231Y", "N198D", "N212D", "K166E"]),
    4: dict(label="BEST_thermal",
            mutations=["D19E", "R73L", "H231Y", "N198D", "N212D", "K166E",
                       "K156E", "K101E"]),
    5: dict(label="BEST_brightness",
            mutations=["D19E", "R73L", "H231Y", "Y237N"]),
    6: dict(label="INSURANCE_supercharge",
            mutations=["K166E", "N212D", "K101E", "K156E", "N198D"]),
}


def apply_mutations(parent: str, mutations: list[str]) -> str:
    """Apply mutations (Met-included numbering) to a parent sequence, verifying
    the reference residue at each position before substituting."""
    import re
    seq = list(parent)
    for mut in mutations:
        m = re.fullmatch(r"([A-Z])(\d+)([A-Z])", mut.strip())
        if m is None:
            raise ValueError(f"Bad mutation: {mut}")
        ref, pos, alt = m.group(1), int(m.group(2)), m.group(3)
        if not (1 <= pos <= len(seq)):
            raise ValueError(f"{mut}: position out of range")
        if seq[pos - 1] != ref:
            raise ValueError(
                f"{mut}: expected {ref} at {pos}, found {seq[pos - 1]}")
        seq[pos - 1] = alt
    return "".join(seq)


def net_charge(seq: str) -> int:
    return sum(CHARGE.get(a, 0) for a in seq)


def delta_charge(mutations: list[str]) -> int:
    return sum(CHARGE.get(m[-1], 0) - CHARGE.get(m[0], 0) for m in mutations)
