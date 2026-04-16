#!/usr/bin/env python3

"""
calculate_features.py

Compute per-patient, per-gene features by integrating original and imputed VCFs.
Features:
- CADD-weighted allele frequency (AF)
- CADD_PHRED statistics (avg, max, sum)
- HFI (High-Function Impact) and LoH metrics
- Average AF

Usage:
  ./cal_features_multi.py \
    --original path/to/original.vcf.gz \
    --imputed  path/to/imputed.vcf.gz \
    --out_dir  /path/to/output_dir
"""

import argparse  # Parsing command-line arguments
from cyvcf2 import VCF  # pyright: ignore[reportMissingImports] # Fast multisample VCF parser
import pandas as pd  # Data handling
import numpy as np   # Numeric operations
from sklearn.impute import KNNImputer  # For CADD imputation
import re  # Regex parsing


# Path to loss-of-function variants list (gnomAD LoF)
default_lof_file = 'Data/gnomad_lofs/gnomad.v2.1.1.all_lofs.hg38.txt'

# Encoding maps for categorical annotations
variant_type_map = {
    'synonymous_SNV': 0, 'unknown': 1, '.': 2,
    'nonframeshift_insertion': 3, 'nonframeshift_deletion': 4,
    'startloss': 5, 'nonsynonymous_SNV': 6,
    'frameshift_insertion': 7, 'frameshift_deletion': 8,
    'stoploss': 9, 'stopgain': 10
}
clnsig_map = {
    '.': 0, 'not_provided': 1, 'no_classification_for_the_single_variant': 2,
    'other': 3, 'Uncertain_significance': 4,
    'protective': 5, 'drug_response': 6,
    'confers_sensitivity': 7, 'Affects': 8,
    'association': 9, 'Benign': 10,
    'Likely_benign': 11, 'Conflicting_classifications_of_pathogenicity': 12,
    'risk_factor': 13, 'Likely_pathogenic': 14, 'Pathogenic': 15
}

GENE_SEP_REGEX = re.compile(r'[;,|/&]')

def parse_gene_field(raw):
    if raw is None:
        return []
    s = str(raw)
    s = re.sub(r'\\x3b|\x3b', ';', s)  # decode literal \x3b to ';'
    parts = [p.strip() for p in GENE_SEP_REGEX.split(s) if p.strip()]
    parts = [re.sub(r'\(.*\)$', '', p) for p in parts]  # drop GENE(123kb)
    return parts

BAD_PREFIX = re.compile(r'^(LOC|LINC|MIR|RP[SL])', re.IGNORECASE)

def pick_primary_gene(genes):
    if not genes:
        return 'NA'
    # prefer names without  „LOC/LINC/MIR/RP*”
    def key(g):
        penalty = 1 if BAD_PREFIX.match(g) else 0
        return (penalty, g)
    return sorted(genes, key=key)[0]


def encode_variant_type(vt):
    """Map variant_type string to numeric code."""
    return variant_type_map.get(vt, 1)


def encode_clnsig(val):
    """Encode CLNSIG to numeric code, splitting on '|' or '/'."""
    if pd.isna(val) or val == '.':
        return 0
    parts = re.split(r'[|/]', str(val))
    codes = [clnsig_map.get(p, -1) for p in parts]
    codes = [c for c in codes if c != -1]
    return max(codes) if codes else 0


def custom_aggregate(group):
    """Aggregate features per gene for one patient."""
    af_vals = group['AF_val'].values
    ph_vals = group['CADD_PHRED'].values
    total_ph = ph_vals.sum()
    # weighted allele frequency (ensure scalar float)
    # weighted allele frequency (ensure scalar float)
    weighted_af = np.average(af_vals, weights=ph_vals/total_ph).item() if total_ph > 0 else 0.0

    # High-Function Impact flag
    hfi_flag = int(
        ((group['variant_type'] == 'nonsynonymous_SNV') &
         ((group['MetaSVM_pred'] == 'D') |
          group['CLNSIG'].str.contains('Pathogenic'))).any()
    )
    # Loss of Heterozygosity proxy
    loh_flag = int(
        (group['variant_type'].isin([
            'frameshift_insertion','frameshift_deletion','stopgain','stoploss'
        ]) |
         (group['LoF'] == 1) |
         group['CLNSIG'].str.contains('Pathogenic')).any()
    )

    # ensure aggregated values are native floats
    # average CADD_PHRED as Python float
    avg_ph = ph_vals.mean().item()
    # max CADD_PHRED as Python float
    max_ph = ph_vals.max().item()
    # sum of CADD_PHRED as Python float
    sum_ph = total_ph  # already float

    avg_af = af_vals.mean().item()

    return pd.Series({
        'CADD_weighted_avg_AF': weighted_af,
        'avg_AF': avg_af,
        'avg_CADD_PHRED': avg_ph,
        'max_CADD_PHRED': max_ph,
        'sum_CADD_PHRED': sum_ph,
        'HFI': hfi_flag,
        'LoH': loh_flag
    })


def load_lof(path):
    """Load LoF variant keys (chrom_pos_ref_alt) into a set."""
    lof_set = set()
    with open(path) as f:
        next(f)
        for line in f:
            cols = line.strip().split('\t')
            lof_set.add('_'.join(cols[0:4]))
    return lof_set


def process_sample(sample, orig_path, imp_ds, variant_info, out_dir):
    """Process one sample: extract AF, CADD, aggregate and write features."""


    vcf_orig = VCF(args.original, regions=args.regions) if args.regions else VCF(args.original)
    vcf_orig.set_threads(4)  # adjust to cores
    
    idx = vcf_orig.samples.index(sample)
    records = []
        # Define output path and header for feature file
    sample_out = f"{out_dir}/{sample}.feature.txt"
    header = (
    "Gene\tCADD_weighted_avg_AF\tavg_CADD_PHRED\tmax_CADD_PHRED\t"
    "sum_CADD_PHRED\tHFI\tLoH\n"
)
    

    for var in vcf_orig:
        if not var.is_snp:
            continue

        chrom, pos, ref, alt = var.CHROM, var.POS, var.REF, var.ALT[0]
        key = f"{chrom}_{pos}_{ref}_{alt}"
        info = variant_info.get(key)
        gene = info['gene']
        ph = info['CADD_PHRED']
        vt = info['variant_type']
        metasvm = info['MetaSVM_pred']
        cln = info['CLNSIG']
        lof_flag = info['LoF']


        # APPROACH: 0/0  --> AF=0,  ./. --> DS
        af = None
        ad = var.format('AD')
        gt = var.genotypes[idx][:2]  # two first alleles

        if ad is not None and len(ad[idx]) > 1:
            counts = ad[idx]
            if counts[1] > 0:
                af = counts[1] / counts.sum()

        if af is None:
            if gt == (0, 0):
                af = 0.0  # homozugous reference → AF = 0
            else:
                ds_arr = imp_ds.get(key)
                if ds_arr is not None:
                    af = ds_arr[idx] / 2

        if af is None or af < 0.2:
            continue


        records.append({
            'Gene': gene,
            'AF_val': af,
            'variant_type': vt,
            'CADD_PHRED': ph,
            'MetaSVM_pred': metasvm,
            'CLNSIG': cln,
            'LoF': lof_flag
        })

    df = pd.DataFrame(records)
    # If no variants for this sample, write only header and return
    if df.empty:
        with open(sample_out, 'w') as f:
            f.write(header)
        return

    # Aggregate features per gene
    result = df.groupby('Gene', group_keys=False).apply(
    lambda g: custom_aggregate(g.drop(columns='Gene'))
    ).reset_index()



    result.to_csv(f"{out_dir}/{sample}.feature.txt", sep='\t', index=False)

    vcf_orig.close()



def parse_args():
    parser = argparse.ArgumentParser(description='Compute per-sample gene features')
    parser.add_argument('--original', required=True)
    parser.add_argument('--imputed', required=True)
    parser.add_argument('--out_dir', required=True)
    parser.add_argument('--regions', help='Optional region, e.g. chr1 or chr1:1-1000000')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    lof_set = load_lof(default_lof_file)

    vcf_orig = VCF(args.original, regions=args.region) if args.regions else VCF(args.original)
    vcf_orig.set_threads(4)
    
    vcf_imp = VCF(args.imputed, regions=args.region) if args.regions else VCF(args.imputed)
    vcf_imp.set_threads(4)

    # Preload DS for all variants into a dictionary
    imp_ds = {}
    for iv in vcf_imp:
        key = f"{iv.CHROM}_{iv.POS}_{iv.REF}_{iv.ALT[0]}"
        imp_ds[key] = iv.format('DS')

    # Collect annotations from INFO

    # Collect annotations from INFO
    variant_info = {}
    for var in vcf_orig:
        chrom, pos, ref, alt = var.CHROM, var.POS, var.REF, var.ALT[0]
        key = f"{chrom}_{pos}_{ref}_{alt}"
        info = var.INFO
        af = var.INFO.get('AF')
        vt = info.get('ExonicFunc.refGene', '.')
        metasvm = info.get('MetaSVM_pred', '.')
        cln = info.get('CLNSIG', '.')
        genes = parse_gene_field(info.get('Gene.refGene', 'NA'))
        gene = pick_primary_gene(genes)

        try:
            ph = float(info.get('CADD1.7_PHRED', '.'))
        except (TypeError, ValueError):
            ph = np.nan

        variant_info[key] = {
            'CADD_PHRED': ph,
            'variant_type_enc': encode_variant_type(vt),
            'MetaSVM_pred_enc': {'D': 1, 'T': 0}.get(metasvm, -1),
            'CLNSIG_enc': encode_clnsig(cln),
            'LoF': int(key in lof_set),
            'variant_type': vt,
            'MetaSVM_pred': metasvm,
            'CLNSIG': cln,
            'gene': gene,
            'global_AF': float(af) if af is not None else 0
        }


    # CADD imputation
    df_variants = pd.DataFrame([
        {
            'key': key,
            **{
                k: v for k, v in val.items()
                if k in ['CADD_PHRED', 'variant_type_enc', 'MetaSVM_pred_enc', 'CLNSIG_enc', 'LoF', 'global_AF']
            }
        }
        for key, val in variant_info.items()
    ])

    imputer = KNNImputer(n_neighbors=5)
    imputed = imputer.fit_transform(df_variants[['CADD_PHRED', 'variant_type_enc', 'MetaSVM_pred_enc', 'CLNSIG_enc', 'LoF', 'global_AF']])

    df_variants['CADD_PHRED'] = imputed[:, 0]

    for idx, row in df_variants.iterrows():
        variant_info[row['key']]['CADD_PHRED'] = row['CADD_PHRED']

    # Process each sample
    for sample in vcf_orig.samples:
        process_sample(sample, args.original, imp_ds, variant_info, args.out_dir)

    print('Done feature calculation')


