#pyBigwig environment
import os
import numpy as np
import pandas as pd
import pyBigWig

# -------------------------------------------------------------------------- 从predicted bigwig files计算predicted counts ---------------------------------------------------------------
# =========================
# 0. 函数
# =========================
def get_cts(peaks_df, bw, width):
    """
    Fetches values from a bigwig bw, given a df with minimally
    chr, start and summit columns. Summit is relative to start.
    Retrieves values of specified width centered at summit.

    "cts" = per base counts across a region
    """
    vals = []
    for i, r in peaks_df.iterrows():
        vals.append(np.nan_to_num(bw.values(r['chr'], 
                                            r['start'] + r['summit'] - width//2,
                                            r['start'] + r['summit'] + width//2)))
        
    return np.array(vals)

# =========================
# 1. cell type lists
# =========================

hs_immune_cell_list=["C10.1", "C10.2", "C10.3", "C10.5", "C10.6", "C10.7", "C10.8", "C10.9", "C10.10", "C10.11", "C10.14", "C3.1", "C3.3"]
TS_immune_cell_list=["C24.1", "C24.3", "C24.4", "C24.8", "C24.10", "C24.11", "C24.12", "C24.14", "C24.15", "C24.16", "C24.17", "C24.18", "C24.21"]

# =========================
# 2. input/output paths
# =========================

hs_peak_df = pd.read_csv("~/chrombpnet/Hs_TS_immune_chrombpnet/human_immune_peak_valid.modified10col.bed",sep="\t",header=None)
hs_peak_df.columns = [
    "chr", "start", "end", "name", "score", "strand",
    "signal", "pvalue", "qvalue", "summit"
]

TS_peak_df=pd.read_csv("~/chrombpnet/Hs_TS_immune_chrombpnet/Treeshrew_immune_peak_valid.modified10col.bed",sep="\t",header=None)
TS_peak_df.columns = [
    "chr", "start", "end", "name", "score", "strand",
    "signal", "pvalue", "qvalue", "summit"
]

for df in [hs_peak_df, TS_peak_df]:
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["summit"] = df["summit"].astype(int)


#human_model_predict
hs_bigwig_dir="~/chrombpnet/Hs_TS_immune_chrombpnet/human_model_predict/pre_bw"

#treeshrew_model_predict
TS_bigwig_dir="~/chrombpnet/Hs_TS_immune_chrombpnet/treeshrew_model_predict/pre_bw"

outdir="~/chrombpnet/Hs_TS_immune_chrombpnet"

# =========================
# 3. shared parameters
# =========================

width=502
folds=range(5)

# =========================
# 4. human model loop
# =========================

#统计人的model预测
#hs_pre_hs_cts_dict={}
hs_pre_hs_total_counts_dict = {}
#hs_pre_hs_logcounts_dict = {}

#hs_pre_ts_cts_dict={}
hs_pre_ts_total_counts_dict = {}
#hs_pre_ts_logcounts_dict = {}

for hs_immune_cell in hs_immune_cell_list:
    for peak in ["human_immune_peak", "treeshrew_immune_peak"]:
        fold_total_counts = []
        for fold in folds:
            bw_path = os.path.join(hs_bigwig_dir, f"fold_{fold}", hs_immune_cell, peak, f"{hs_immune_cell}_fold_{fold}_{peak}_prebw_chrombpnet.bw")
            with pyBigWig.open(bw_path) as bw:
                if peak == "human_immune_peak":
                    cts = get_cts(hs_peak_df, bw, width=width)
                    #hs_pre_hs_cts_dict[hs_immune_cell] = cts
                    #hs_pre_hs_total_counts_dict[hs_immune_cell] = cts.sum(axis=1)
                    #hs_pre_hs_logcounts_dict[hs_immune_cell] = np.log1p(cts.sum(axis=1))
                else:
                    cts = get_cts(TS_peak_df, bw, width=width)
                    #hs_pre_ts_cts_dict[hs_immune_cell] = cts
                    #hs_pre_ts_total_counts_dict[hs_immune_cell] = cts.sum(axis=1)
                    #hs_pre_ts_logcounts_dict[hs_immune_cell] = np.log1p(cts.sum(axis=1))
                
                total_counts = cts.sum(axis=1)
                fold_total_counts.append(total_counts)

        mean_total_counts = np.mean(fold_total_counts, axis=0)
        if peak == "human_immune_peak":
            hs_pre_hs_total_counts_dict[hs_immune_cell] = mean_total_counts
        else:
            hs_pre_ts_total_counts_dict[hs_immune_cell] = mean_total_counts

#hs_pre_hs
hs_pre_hs_total_df = hs_peak_df[["chr", "start", "end"]].copy()

for celltype, total_counts in hs_pre_hs_total_counts_dict.items():
    hs_pre_hs_total_df[f"{celltype}_pre_counts"] = total_counts
    hs_pre_hs_total_df[f"{celltype}_pre_logcounts"] = np.log1p(total_counts)

hs_pre_hs_total_df.to_csv(os.path.join(outdir,"hs_predicted_hs_counts_immune_celltypes.csv"), index=False)

#hs_pre_ts
hs_pre_ts_total_df = TS_peak_df[["chr", "start", "end"]].copy()

for celltype, total_counts in hs_pre_ts_total_counts_dict.items():
    hs_pre_ts_total_df[f"{celltype}_pre_counts"] = total_counts
    hs_pre_ts_total_df[f"{celltype}_pre_logcounts"] = np.log1p(total_counts)

hs_pre_ts_total_df.to_csv(os.path.join(outdir,"hs_predicted_ts_counts_immune_celltypes.csv"), index=False)

# =========================
# 6. tree shrew model loop
# =========================

#统计树鼩的model预测
#TS_pre_hs_cts_dict={}
TS_pre_hs_total_counts_dict = {}
#TS_pre_hs_logcounts_dict = {}

#TS_pre_ts_cts_dict={}
TS_pre_ts_total_counts_dict = {}
#TS_pre_ts_logcounts_dict = {}

for TS_immune_cell in TS_immune_cell_list:
    for peak in ["human_immune_peak", "treeshrew_immune_peak"]:
        fold_total_counts = []
        for fold in folds:
            bw_path = os.path.join(TS_bigwig_dir, f"fold_{fold}", TS_immune_cell, peak,f"{TS_immune_cell}_fold_{fold}_{peak}_prebw_chrombpnet.bw")
            with pyBigWig.open(bw_path) as bw:
                if peak == "human_immune_peak":
                    cts = get_cts(hs_peak_df, bw, width=width)
                    #TS_pre_hs_cts_dict[TS_immune_cell] = cts
                    #TS_pre_hs_total_counts_dict[TS_immune_cell] = cts.sum(axis=1)
                    #TS_pre_hs_logcounts_dict[TS_immune_cell] = np.log1p(cts.sum(axis=1))
                else:
                    cts = get_cts(TS_peak_df, bw, width=width)
                    #TS_pre_ts_cts_dict[TS_immune_cell] = cts
                    #TS_pre_ts_total_counts_dict[TS_immune_cell] = cts.sum(axis=1)
                    #TS_pre_ts_logcounts_dict[TS_immune_cell] = np.log1p(cts.sum(axis=1))
                total_counts = cts.sum(axis=1)
                fold_total_counts.append(total_counts)

        mean_total_counts = np.mean(fold_total_counts, axis=0)
        if peak == "human_immune_peak":
            TS_pre_hs_total_counts_dict[TS_immune_cell] = mean_total_counts
        else:
            TS_pre_ts_total_counts_dict[TS_immune_cell] = mean_total_counts

#TS_pre_hs
TS_pre_hs_total_df = hs_peak_df[["chr", "start", "end"]].copy()

for celltype, total_counts in TS_pre_hs_total_counts_dict.items():
    TS_pre_hs_total_df[f"{celltype}_pre_counts"] = total_counts
    TS_pre_hs_total_df[f"{celltype}_pre_logcounts"] = np.log1p(total_counts)

TS_pre_hs_total_df.to_csv(os.path.join(outdir,"TS_predicted_hs_counts_immune_celltypes.csv"), index=False)

#TS_pre_ts
TS_pre_ts_total_df = TS_peak_df[["chr", "start", "end"]].copy()

for celltype, total_counts in TS_pre_ts_total_counts_dict.items():
    TS_pre_ts_total_df[f"{celltype}_pre_counts"] = total_counts
    TS_pre_ts_total_df[f"{celltype}_pre_logcounts"] = np.log1p(total_counts)

TS_pre_ts_total_df.to_csv(os.path.join(outdir,"TS_predicted_ts_counts_immune_celltypes.csv"), index=False)