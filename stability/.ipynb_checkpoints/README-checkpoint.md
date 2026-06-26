# Stability axis — supercharging rationale and optional ThermoMPNN validation

The brightness axis of this project is fully data-driven and reproducible
(`src/brightness_model.py`, R^2 ~ 0.93 on held-out DMS variants). This note
documents the **stability axis**, which is what the 72 C retention score rewards.

## Why the combined score is a stability problem

The competition's per-sequence score is

```
Relative_Brightness x Thermal_Retention
   = (F_initial / F_initial_WT) x (F_final / F_initial)
   = F_final / F_initial_WT
```

so the ranking quantity is the **absolute post-heat brightness**. In GFP,
fluorescence is gated by folding stability through a threshold function
(Sarkisyan et al. 2016): brightness collapses once cumulative destabilisation
crosses ~7-9 kcal/mol. Hence both `F_initial` (how completely the barrel folds in
the cell-free system) and the 72 C retention (how well it resists unfolding) are
stability read-outs at two temperatures.

## The supercharging recipe (and why it is folding-safe here)

Thermal Green Protein (Close et al. 2015) was made by replacing exposed K/R/N/Q
residues on the beta-barrel with glutamate, taking the net charge to -10 and
greatly improving thermal stability and aggregation resistance with negligible
change to the chromophore (spectra and quantum yield essentially unchanged).
GFPs tolerate much more extreme supercharging (down to about -30) while remaining
fluorescent (Lawrence et al. 2007; Der et al. 2013).

We select each surface charge mutation (K/R/N/Q -> E/D) under one hard constraint:
**it must be brightness-neutral or positive in the DMS data.** In a
stability-gated assay, brightness-neutral == the mutation does not destabilise the
native fold == it is "free" in the stability budget. The added negative surface
charge then contributes thermal/aggregation resistance "on top" of an undisturbed
fold. This is why our six designs raise net charge from -6 (sfGFP) to between -7
and -17 without any predicted loss of initial brightness.

## Optional: quantitative ddG validation with ThermoMPNN (run in Colab)

ThermoMPNN (Dieckhaus et al. 2024, PNAS) predicts per-mutation folding ddG from a
structure. Use it to (a) confirm the supercharging mutations are ddG-neutral or
stabilising and (b) keep each design's cumulative ddG net-stabilising, well inside
the Sarkisyan threshold. It needs a PDB structure and downloadable weights, so it
is easiest to run on Google Colab (PDB / model hosts are reachable there but were
blocked in the offline environment used to build this repo).

```bash
# In a Colab cell (GPU runtime recommended):
!git clone https://github.com/Kuhlman-Lab/ThermoMPNN.git
%cd ThermoMPNN
!pip install -r requirements.txt

# Fetch the sfGFP structure (recommended PDB from the data pack):
!wget https://files.rcsb.org/download/2B3P.pdb

# Run the single-mutant ddG scan over 2B3P, then look up our positions
# (Met-included numbering -> match the chain numbering in 2B3P), and verify:
#   - every supercharging mutation has ddG <= ~0 (neutral or stabilising)
#   - sum of ddG per design stays well below the ~7-9 kcal/mol fluorescence
#     threshold of Sarkisyan et al. (2016)
```

Record the resulting per-mutation ddG table here when run, and cite it in the
"Algorithm design introduction" PDF as the structural confirmation of the
DMS-based folding-safety argument.
