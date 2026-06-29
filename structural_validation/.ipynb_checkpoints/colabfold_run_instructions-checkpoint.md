# Structural validation with ColabFold

This step was performed in a **Google Colab GPU runtime** using the official ColabFold AlphaFold2_batch notebook.

ColabFold is not maintained as a runnable local script inside this repository. The repository only stores:

- the FASTA inputs used for the run;
- the downloaded ColabFold result files;
- the downstream parser used to summarize the saved results.

This avoids committing an exported Colab notebook as a `.py` file, because such exports often depend on the current Google Colab runtime and can break when run later as standalone scripts.

## Official ColabFold resources

```text
Official ColabFold GitHub:
https://github.com/sokrypton/ColabFold

Official AlphaFold2_batch notebook:
https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/batch/AlphaFold2_batch.ipynb

Notebook source on GitHub:
https://github.com/sokrypton/ColabFold/blob/main/batch/AlphaFold2_batch.ipynb
```

## Input FASTA files

The repository contains two equivalent ColabFold input formats.

### Option A — individual FASTA files

These files were used in the Google Drive batch layout:

```text
input/WT_sfGFP.fasta
input/Seq1_conservative_surface_D19E_R73L_H231Y.fasta
input/Seq2_balanced_core_D19E_R73L_H231Y_K166E_N212D.fasta
input/Seq3_MAIN_balanced_D19E_R73L_H231Y_N198D_N212D_K166E.fasta
input/Seq4_BEST_thermal_D19E_R73L_H231Y_N198D_N212D_K166E_K156E_K101E.fasta
input/Seq5_BEST_brightness_D19E_R73L_H231Y_Y237N.fasta
input/Seq6_INSURANCE_supercharge_K166E_N212D_K101E_K156E_N198D.fasta
```

In Google Drive, copy them to:

```text
/content/drive/MyDrive/colab_fold/input_fasta/
```

The result directory used in the original run was:

```text
/content/drive/MyDrive/colab_fold/result/
```

### Option B — one multi-FASTA file

For manual upload to the official batch notebook, use:

```text
structural_validation/colabfold_input_WT_plus_6designs.fasta
```

This file contains `WT_sfGFP` and the six final GeneMeow designs.

## Recommended AlphaFold2_batch settings

Use the official notebook and set:

```text
msa_mode      = MMseqs2 (UniRef+Environmental)
model_type    = AlphaFold2-ptm
num_models    = 5
num_recycles  = 3
use_templates = False
num_relax     = 0
rank_by       = pLDDT
zip_results   = True
```

The exact option names can change in the official ColabFold notebook, so the official notebook interface should be treated as the source of truth.

## Saved outputs in this repository

The already downloaded ColabFold outputs used in this project are stored in:

```text
models/colabfold_rank001_final6/
figures/colabfold_raw_plots_final6/
```

The summarized structural metrics are stored in:

```text
outputs/colabfold_validation_metrics.csv
outputs/final6_colabfold_brightness_metrics.csv
```

## Downstream processing of saved ColabFold results

After downloading the ColabFold result archive from Google Colab, the project-specific parser can be run from the `structural_validation/` directory:

```bash
python process_colabfold_results_final6_memory_safe.py
```

This parser is memory-safe and summarizes saved ColabFold outputs into per-design metrics such as mean pLDDT, pTM, chromophore-region confidence, and core Cα RMSD to WT.

## Interpretation

ColabFold was used only as a structural sanity check. High pLDDT and preservation of the GFP barrel are necessary but not sufficient for experimental brightness or thermal stability. Brightness and post-heat fluorescence retention still require experimental measurement.
