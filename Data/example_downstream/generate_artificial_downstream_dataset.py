#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd

np.random.seed(42)

out_dir = Path("Data/example_downstream")
out_dir.mkdir(parents=True, exist_ok=True)

datasets = ["DatasetA", "DatasetB", "DatasetC"]
kits = ["KitA", "KitB", "KitC"]

samples_per_dataset = 40
n_genes = 200

metadata_records = []

for dataset, kit in zip(datasets, kits):
    for i in range(1, samples_per_dataset + 1):
        metadata_records.append({
            "Sample": f"{dataset}_sample{i:03d}",
            "Dataset": dataset,
            "Kit": kit
        })

metadata = pd.DataFrame(metadata_records)
genes = [f"G{i:03d}" for i in range(1, n_genes + 1)]

# Systematic kit-specific under-detection.
# These Sample x Gene records are omitted, not set to zero.
kit_dropout_genes = {
    "KitA": set(genes[50:80]),    # G051-G080 missing in KitA
    "KitB": set(genes[80:110]),   # G081-G110 missing in KitB
    "KitC": set(genes[110:140])   # G111-G140 missing in KitC
}

def draw_cwaf(a, b):
    value = float(np.random.beta(a, b))
    return round(max(value, 0.001), 6)

records = []
systematic_missing = 0
random_missing = 0

for _, row in metadata.iterrows():
    sample = row["Sample"]
    dataset = row["Dataset"]
    kit = row["Kit"]

    for gene in genes:
        # Systematic kit-specific missingness: no record is written.
        if gene in kit_dropout_genes.get(kit, set()):
            systematic_missing += 1
            continue

        # Shared genes, detected in almost all samples.
        if gene in genes[:50]:
            cwaf = draw_cwaf(5, 2)

        # Genes affected by kit-specific dropout in another kit.
        elif gene in genes[50:140]:
            if np.random.rand() < 0.02:
                random_missing += 1
                continue
            cwaf = draw_cwaf(4, 3)

        # Partially detected noisy genes.
        elif gene in genes[140:170]:
            if np.random.rand() < 0.35:
                random_missing += 1
                continue
            cwaf = draw_cwaf(3, 4)

        # Low-detection background genes.
        else:
            if np.random.rand() < 0.75:
                random_missing += 1
                continue
            cwaf = draw_cwaf(2, 5)

        records.append({
            "Sample": sample,
            "Gene": gene,
            "CADD_weighted_avg_AF": cwaf,
            "Dataset": dataset,
            "Kit": kit
        })

raw_ft_long = pd.DataFrame(records)

metadata.to_csv(out_dir / "metadata.tsv", sep="\t", index=False)
raw_ft_long.to_csv(out_dir / "raw_ft_long.tsv", sep="\t", index=False)

print(f"Wrote: {out_dir / 'metadata.tsv'}")
print(f"Wrote: {out_dir / 'raw_ft_long.tsv'}")
print(f"Samples: {metadata.shape[0]}")
print(f"Genes: {len(genes)}")
print(f"Rows written: {raw_ft_long.shape[0]}")
print(f"Systematic missing records: {systematic_missing}")
print(f"Random missing records: {random_missing}")
print(f"Total possible Sample x Gene pairs: {metadata.shape[0] * len(genes)}")