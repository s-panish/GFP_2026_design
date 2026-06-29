# ColabFold official workflow reference

This repository does not maintain a custom runnable ColabFold script. The structural prediction step should be run through the official ColabFold notebook.

Official resources:

```text
ColabFold GitHub:
https://github.com/sokrypton/ColabFold

AlphaFold2_batch notebook in Google Colab:
https://colab.research.google.com/github/sokrypton/ColabFold/blob/main/batch/AlphaFold2_batch.ipynb

Notebook source on GitHub:
https://github.com/sokrypton/ColabFold/blob/main/batch/AlphaFold2_batch.ipynb
```

Use the FASTA files in `input/` or the multi-FASTA file in `structural_validation/colabfold_input_WT_plus_6designs.fasta`.

The saved results from the original Google Colab run are already included in `models/colabfold_rank001_final6/`.
