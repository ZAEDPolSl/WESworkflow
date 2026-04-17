# functionspython/run_parc.py
import numpy as np
import pandas as pd
import parc

def parc_labels_from_long_tsv(path, knn=30, jac=0.15, res=1.0, seed=42, n_jobs=-1):
    df = pd.read_csv(path, sep="\t", header=0)
    if "Gene" in df.columns:
        df = df.rename(columns={"Gene": "gene"})
    elif "gene" not in df.columns:
        raise ValueError("Missing required column: 'gene' or 'Gene'")
    need = {"Sample", "gene", "CADD_weighted_avg_AF"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"Missing columns: {missing}")

    wide = (df.groupby(["Sample","gene"], as_index=False)["CADD_weighted_avg_AF"].mean()
              .pivot(index="Sample", columns="gene", values="CADD_weighted_avg_AF")
              .sort_index().fillna(0.0))
    X = wide.to_numpy(dtype=np.float32)
    P = parc.PARC(X, knn=int(knn), jac_std_global=jac, resolution_parameter=res,
                  random_seed=seed)
    P.run_PARC()
    return list(wide.index), list(map(int, P.labels))