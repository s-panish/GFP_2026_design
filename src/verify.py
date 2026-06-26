"""
verify.py
---------
Hard gate before submission. Reproduces every disqualification rule from the
competition brief and the data pack.
"""

from __future__ import annotations
import re
import csv
import openpyxl

ALLOWED_AA = set("ACDEFGHIKLMNPQRSTVWY")
CHROMOPHORE = (65, 66, 67)  # TYG in sfGFP


def validate_format(seq: str) -> list[str]:
    """Return a list of rule violations; empty list == valid."""
    errors = []
    if not (220 <= len(seq) <= 250):
        errors.append(f"length {len(seq)} not in 220-250")
    if not seq.startswith("M"):
        errors.append("does not start with M")
    bad = sorted(set(seq) - ALLOWED_AA)
    if bad:
        errors.append("non-standard characters: " + "".join(bad))
    if "*" in seq:
        errors.append("contains stop '*'")
    return errors


def chromophore_intact(seq: str, parent: str) -> bool:
    return all(seq[p - 1] == parent[p - 1] for p in CHROMOPHORE)


def load_exclusion(path: str) -> set[str]:
    """Robustly load excluded amino-acid sequences from Exclusion_List.csv."""
    excluded = set()
    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader, None)
        for row in reader:
            for cell in row:
                s = re.sub(r"\s+", "", str(cell)).upper()
                if len(s) >= 50 and set(s) <= ALLOWED_AA:
                    excluded.add(s)
    return excluded


def load_beforetop(xlsx_path: str) -> set[str]:
    wb = openpyxl.load_workbook(xlsx_path, read_only=True, data_only=True)
    ws = wb["beforetopseqs"]
    out = set()
    for _yr, seq in ws.iter_rows(min_row=2, values_only=True):
        if seq:
            out.add(re.sub(r"\s+", "", str(seq)).upper())
    return out
