# Structural validation (ColabFold) — run it on WT + our six designs

The third, independent line of evidence: a structural QC check that the designed
substitutions preserve the GFP beta-barrel fold. It does **not** predict brightness or
heat retention (use the brightness gate and the supercharging rationale for those).

## Input (already prepared)

`colabfold_input_WT_plus_6designs.fasta` — a single multi-FASTA with `WT_sfGFP` plus the
six designs (`Seq1_…` … `Seq6_…`).

## Run ColabFold (AlphaFold2-batch)

Official batch notebook:
<https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/AlphaFold2_batch.ipynb>

Upload the FASTA and set:

| Setting | Value |
|---|---|
| msa_mode | `mmseqs2_uniref_env` |
| model_type | `alphafold2_ptm` |
| num_models | 5 |
| num_recycles | 3 |
| use_templates | off |
| use_amber (relaxation) | off |
| rank_by | pLDDT |
| zip_results | on |

Scripted alternative:

```bash
!pip -q install "colabfold[alphafold-minus-jax]"
!colabfold_batch --num-models 5 --num-recycle 3 --rank plddt \
    --msa-mode mmseqs2_uniref_env --model-type alphafold2_ptm \
    colabfold_input_WT_plus_6designs.fasta colabfold_out/
```

## Process the results

Upload the ColabFold results zip and run:

```bash
python process_colabfold_results_final6_memory_safe.py
```

It is memory-safe (streams/extracts incrementally) and writes:

- `stage4_colabfold_metrics_final6.csv` — per design: mean pLDDT, pTM, chromophore
  pLDDT (min/mean over 65–67), core pLDDT and **Cα core RMSD to WT over residues 1–220**.
- `final6_structural_summary.csv` — compact summary joined with the brightness gate.
- figures: mean pLDDT, min chromophore pLDDT, core RMSD to WT, brightness delta, and the
  raw ColabFold pLDDT / PAE / MSA-coverage plots.
- `GeneMeow_v3_figures_for_PDF.zip` and a minimal submission package.

## What to expect

All designs are point mutants of WT at surface positions, so AF2 should return
high-confidence GFP barrels (mean pLDDT ~96, high chromophore-region confidence,
sub-Ångström core RMSD to WT). High confidence is necessary but **not** sufficient: it
confirms the fold is intact, not that the protein is bright or thermostable. Drop the
resulting metrics CSV and figures into `docs/` to populate the PDF's structural section.
