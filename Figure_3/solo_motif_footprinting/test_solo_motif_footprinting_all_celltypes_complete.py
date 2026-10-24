
# SETUP ------------------------------------------------------------------------

import sys
import os
import re
import numpy as np
import pandas as pd
import pickle
import gc
import logging
import datetime

from pyfaidx import Fasta

import torch
import tangermeme
import tangermeme.io
from tangermeme.ersatz import substitute
from tangermeme.predict import predict

import matplotlib
matplotlib.use("Agg")   # 后台运行时不要弹窗
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import bpnetlite
from bpnetlite import BPNet

sys.path.append("~/HDMA/code/03-chrombpnet")
import tangermeme_utils as tanutils
from tangermeme_utils.constants import *
from tangermeme_utils.eval import summarize_across_folds, mean_across_sequences
from tangermeme_utils.wrappers import ChromBPNetWrapper, reshape_folds


print(torch.__version__)
print(tangermeme.__version__)
print(bpnetlite.__version__)
print(np.__version__)

torch.manual_seed(100)

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["ps.fonttype"] = 42


# CONFIGS ----------------------------------------------------------------------
# 每个字典代表一个要跑的细胞大类中的一个 cluster

BASE_DIR = "~/chrombpnet/chrombpnet_use_snakefile/Treeshrew"
REFS_DIR = "~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta"

# 组合图设置：只影响最后的跨细胞类型 footprint 图
PLOT_WINDOW = 200          # 显示 motif 中心前后 200 bp
COMBINED_NCOLS = 4         # 每页放 4 列子图
YLIM_PADDING_FRAC = 0.05   # ylim 上下各增加 5% 空白

genome_fa_path = os.path.join(REFS_DIR, "genome.fa")
common_motif_dict={
'NR5A2':'TGACCTTGA',
'RFX2':'GTTGCCATGGCAAC',
'MEF2D':'ACTATAAATAGA',
'OTX2':'GGGATTA',
'TCF7L2':'AGATCAAAG',
'GATA3':'TCTTATCT',
'ETV6':'GCGGAAGTG',
'P73':'ACATGTCTGGACATGT',
'KLF3':'GACCACGCCC',
'HNF1B':'GTTAATGATTAAC',
'NGN2':'CAGATGG',
'PAX2':'TCAGTGAAGCGTGACG',
'NKX2-8':'CCACTTGA',
'TWIST1':'CCAGATGT',
'MYOD1':'GCACCTGTC',
'FOSB+JUN':'GATGACGTCAT',
'POU4F1':'TGAATAATTAAT',
'LHX2':'CAATTA',
'EBF1':'TCCCCAGGGGA',
'ASCL1':'GCAGCTGC',
'SIX1':'GTAACCTGA',
'SP2':'GGGGCGGGG',
'PAX6':'TTCACGCATGAGTT',
'SPI1':'AAAGAGGAAGTGG',
'NFYC':'CCAATCA',
'TEAD1':'ACATTCCAG',
'POU6F2':'AGCTCATTA',
'RFX1':'CGTTGCCATGGCAACG',
'NFYA':'CCAATCA',
'CTCF':'GCCACCAGGGGGCGC',
'SP1':'GGGGCGGGG',
'NFIA':'TGCCAA',
'FOXB1':'TATGTAAATAT',
'SOX9':'CCATTGTT',
'NFYB':'TCATTGGCC',
'KLF12':'GGGGCGGGG',
'ZNF76':'ACCCACAATGCACCGCG',
'NRF1':'CTGCGCCTGCGC'
}

CONFIGS = [
    {
        "cell_class": "TS_C1Adrenal_cells",
        "cluster": "C1.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C2PNS",
        "cluster": "C2.5",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C3Cardiomyocytes",
        "cluster": "C3.2",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C4NAG",
        "cluster": "C4.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C5CNS_neurons",
        "cluster": "C5.16",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C6Definitive_erythroid_cells",
        "cluster": "C6.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C7Endothelial_cells",
        "cluster": "C7.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C8Epithelial_cells",
        "cluster": "C8.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C9Gastrointestinal_tract",
        "cluster": "C9.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C10Hepatocytes",
        "cluster": "C10.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C11Intermediate_neuronal_progenitors",
        "cluster": "C11.2",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C12Kidney_cells",
        "cluster": "C12.4",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C13Lung_cells",
        "cluster": "C13.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C14Mesoderm",
        "cluster": "C14.7",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C15Muscle_cells",
        "cluster": "C15.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C16Neural_crest_PNS_glia",
        "cluster": "C16.5",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C17Neural_crest_PNS_neurons",
        "cluster": "C17.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C18Neuroectoderm_and_glia",
        "cluster": "C18.20",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C19Olfactory_sensory_neurons",
        "cluster": "C19.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C20Oligodendrocytes",
        "cluster": "C20.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C21Otic_cells",
        "cluster": "C21.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C22Primitive_erythroid_cells",
        "cluster": "C22.1",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C23Retinal_cells",
        "cluster": "C23.6",
        "motifs": common_motif_dict
    },
    {
        "cell_class": "TS_C24white_blood_cells",
        "cluster": "C24.3",
        "motifs": common_motif_dict
    },
]


# LOGGING ----------------------------------------------------------------------

logger = logging.getLogger()
logger.setLevel(logging.WARNING)

ch1 = logging.StreamHandler()
ch1.setLevel(logging.WARNING)

formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)

ch1.setFormatter(formatter)
logger.addHandler(ch1)


# LOAD GENOME ------------------------------------------------------------------

logger.critical("Loading genome...")
genome_fa = Fasta(genome_fa_path)


#
def is_acgt_region(row, genome_fa):
    chrom = row[0]
    start = int(row[1])
    end = int(row[2])

    if chrom not in genome_fa:
        return False

    if start < 0 or end > len(genome_fa[chrom]):
        return False

    seq = genome_fa[chrom][start:end].seq.upper()

    return set(seq).issubset({"A", "C", "G", "T"})


# MAIN LOOP --------------------------------------------------------------------
# 外层循环：每次处理一个细胞大类中的一个 cluster

# 保存所有细胞类型的 footprint 结果：
# all_motif_delta["MEF2D"] 中会保存 MEF2D 在所有细胞类型中的结果
all_motif_delta = {motif_name: [] for motif_name in common_motif_dict}

for config in CONFIGS:

    cell_class = config["cell_class"]
    cluster = config["cluster"]
    input_motifs_solo = config["motifs"]


    PROJ_IN = os.path.join(BASE_DIR, cell_class)
    PROJ_OUT = os.path.join(BASE_DIR, cell_class)

    out = os.path.join(
        PROJ_OUT,
        f"output/predict/{cluster}/in_silico_test_solo_completeResults"
    )

    os.makedirs(out, exist_ok=True)

    logger.critical("========================================================")
    logger.critical(f"Processing cell class: {cell_class}")
    logger.critical(f"Processing cluster: {cluster}")
    logger.critical(f"Project input directory: {PROJ_IN}")
    logger.critical(f"Project output directory: {PROJ_OUT}")
    logger.critical(f"References directory: {REFS_DIR}")
    logger.critical(f"Output directory: {out}")
    logger.critical(f"GPU available: {torch.cuda.is_available()}")
    logger.critical("========================================================")


    # -------------------------------------------------------------------------
    # 对于一个细胞类型执行一次生成 random 100 negative bed
    # -------------------------------------------------------------------------

    negatives_dir = os.path.join(PROJ_IN, "data/nonpeaks")
    output_dir = os.path.join(PROJ_IN, "data/nonpeaks")
    folds = [0, 1, 2, 3, 4]

    print("------------" + cluster + "------------")

    nonpeaks = []

    for fold in folds:

        nonpeaks_bed = os.path.join(
            negatives_dir,
            f"fold_{fold}/{cluster}_negatives.bed"
        )

        nonpeaks_fold = pd.read_csv(nonpeaks_bed, sep="\t", header=None)

        # 随机打乱 index，只检查到找到 20 条 clean nonpeaks 为止
        np.random.seed(123 + fold)
        candidate_indices = np.random.permutation(
            nonpeaks_fold.index.to_numpy()
        )

        selected_indices = []

        for idx in candidate_indices:

            row = nonpeaks_fold.loc[idx]

            if is_acgt_region(row, genome_fa):
                selected_indices.append(idx)

            if len(selected_indices) == 20:
                break

        if len(selected_indices) < 20:
            raise ValueError(
                f"{cluster}, fold {fold} only found "
                f"{len(selected_indices)} valid nonpeaks"
            )

        # 注意这里用 loc，不用 iloc
        nonpeaks_subset = nonpeaks_fold.loc[selected_indices]

        nonpeaks.append(nonpeaks_subset)

        logger.critical(
            f"{cluster}, fold {fold}: selected "
            f"{len(nonpeaks_subset)} valid nonpeaks"
        )

    nonpeaks = pd.concat(nonpeaks, ignore_index=True)

    rand100_bed = os.path.join(
        output_dir,
        f"{cluster}__output_negatives_rand100.bed"
    )

    nonpeaks.to_csv(
        rand100_bed,
        sep="\t",
        header=False,
        index=False
    )

    logger.critical(f"Saved random negatives: {rand100_bed}")


    # -------------------------------------------------------------------------
    # RUN SOLO MOTIF ISM
    # -------------------------------------------------------------------------

    logger.critical(f"--- Processing Cluster: {cluster} ---")

    base_models = []

    logger.critical("\tLoading base BPNet models for all folds...")

    for fold_idx in range(5):

        model_path = os.path.join(
            PROJ_IN,
            f"output/models/fold_{fold_idx}/{cluster}/models/chrombpnet_nobias.h5"
        )

        logger.critical(f"\tLoading model: {model_path}")

        model = BPNet.from_chrombpnet(filename=model_path)
        base_models.append(model)


    # Load background sequences ONCE per cluster
    logger.critical("\tLoading background sequences...")

    nonpeaks_bed = os.path.join(
        PROJ_IN,
        f"data/nonpeaks/{cluster}__output_negatives_rand100.bed"
    )

    X_nonpeaks = tangermeme.io.extract_loci(
        nonpeaks_bed,
        genome_fa_path,
        ignore=CHAR_IGNORE
    )

    X_nonpeaks = X_nonpeaks.float()


    path_solo_pkl = os.path.join(out, "solo_predictions.pkl")
    path_solo_df = os.path.join(out, "summarized_effects_solo.tsv")


    # solo effects -------------------------------------------------------------

    logger.critical("\tRunning solo motif ISM...")
    logger.critical(input_motifs_solo)

    solo_dfs, solo_preds = [], {}
    models = [ChromBPNetWrapper(m) for m in base_models]

    with torch.no_grad():#关闭梯度计算，仅进行模型预测，不需要反向传播和模型训练。不保存梯度，减少 GPU 显存，加快预测。

        for key, motif in input_motifs_solo.items():

            logger.critical(f"\t\t{key}: {motif}")

            X_after_solo = substitute(X_nonpeaks, motif)

            y_before_list, y_after_list = [], []

            for model in models:

                y_before_fold = predict(model, X_nonpeaks)
                y_after_fold = predict(model, X_after_solo)

                y_before_list.append((y_before_fold[0], y_before_fold[1]))
                y_after_list.append((y_after_fold[0], y_after_fold[1]))

            del X_after_solo #当前 motif 已经预测完成，所以删除，避免多个 motif 的序列张量同时占用内存
            torch.cuda.empty_cache() #释放 PyTorch 缓存中不再使用的显存

            y_before = reshape_folds(y_before_list)
            y_after = reshape_folds(y_after_list)

            y_after_mean, y_after_sd = summarize_across_folds(
                mean_across_sequences(y_after[1])
            )#把当前 motif 插入100条背景序列后，先在每个 fold 内对100条序列取平均，再对5个 fold 取平均得到的预测 counts。主要表示不同模型 fold 之间的不确定性，而不是100条背景序列之间的离散程度

            y_effect = y_after[1] - y_before[1]#对于每个 fold、每条背景序列都单独做配对相减

            y_effect_mean, y_effect_sd = summarize_across_folds(
                mean_across_sequences(y_effect)
            )#每个 fold 内对100条序列求平均，再对5个 fold 求均值和标准差

            #保存完整预测：所有 fold，所有背景序列
            solo_preds[key] = {
                "y_before": y_before,
                "y_after": y_after
            }
            #为当前 motif 创建汇总表格
            df = pd.DataFrame({
                "key": [key],
                "motif": [motif],
                "counts_after_mean": y_after_mean.squeeze().tolist(),
                "counts_after_sd": y_after_sd.squeeze().tolist(),
                "effect_mean": y_effect_mean.squeeze().tolist(),
                "effect_sd": y_effect_sd.squeeze().tolist()
            })

            solo_dfs.append(df)

        solo_dfs = pd.concat(solo_dfs, ignore_index=True)#把多个单行 DataFrame 合并

        with open(path_solo_pkl, "wb") as f:
            pickle.dump(solo_preds, f)

        solo_dfs.to_csv(path_solo_df, sep="\t", index=False)

        logger.critical(f"Saved: {path_solo_pkl}")
        logger.critical(f"Saved: {path_solo_df}")


    # average the profiles -----------------------------------------------------

    solo_profile_preds_summarized = {}

    for i, key in enumerate(solo_preds.keys()):

        y_before_profile_means, y_before_profile_sd = summarize_across_folds(
            mean_across_sequences(solo_preds[key]["y_before"][0])
        )

        y_after_profile_means, y_after_profile_sd = summarize_across_folds(
            mean_across_sequences(solo_preds[key]["y_after"][0])
        )

        solo_profile_preds_summarized[key] = {
            "y_before_means": y_before_profile_means,
            "y_before_sd": y_before_profile_sd,
            "y_after_means": y_after_profile_means,
            "y_after_sd": y_after_profile_sd
        }

        logger.critical(
            f"{key}: "
            f"{solo_profile_preds_summarized[key]['y_before_means'].shape}"
        )


    # collect footprints for the final combined plots -------------------------
    # 这里不再立即画每个细胞类型的图，而是先保存当前细胞类型的结果。
    # 等所有细胞类型都处理完后，再按 motif 统一计算 ylim 并绘图。

    for motif_name in solo_profile_preds_summarized.keys():

        before = (
            solo_profile_preds_summarized[motif_name]["y_before_means"]
            .detach()#从 PyTorch 计算图中分离
            .cpu()#将张量从 GPU 移到 CPU
            .squeeze()#删除长度为 1 的维度
            .numpy()#转换成 NumPy 数组，便于后续处理和绘图
        )

        after = (
            solo_profile_preds_summarized[motif_name]["y_after_means"]
            .detach()
            .cpu()
            .squeeze()
            .numpy()
        )

        delta = after - before

        x = np.arange(len(delta))#delta 长度是1000时，x = 0, 1, 2, ..., 999
        center = len(delta) // 2 #500
        x_relative = x - center #转换成相对中心点的坐标，x-relative=-500, -499, ..., -1, 0, 1, ..., 499

        # 只保留 motif 中心 -200 到 +200 bp，包含两端，共 401 个位置
        keep = (
            (x_relative >= -PLOT_WINDOW)
            & (x_relative <= PLOT_WINDOW)
        )

        all_motif_delta[motif_name].append({
            "cell_class": cell_class,
            "cluster": cluster,
            "x": x_relative[keep].copy(),
            "delta": delta[keep].copy()
        })

    logger.critical(f"Done with {cell_class} {cluster}")

    #删除变量
    del base_models
    del models
    del X_nonpeaks
    del solo_preds
    del solo_profile_preds_summarized

    torch.cuda.empty_cache() #释放 PyTorch 中没有变量继续引用的 GPU 缓存
    gc.collect() #触发 Python 垃圾回收，释放 CPU 内存中不再使用的对象


# COMBINED PLOTTING ------------------------------------------------------------
# 到这里时，全部细胞类型已经处理完成。
# 现在每个 motif 单独计算一个共享 ylim，然后把所有细胞类型画在同一页。

combined_plot_dir = os.path.join(
    BASE_DIR,
    "all_celltypes_solo_motif_plots"
)
os.makedirs(combined_plot_dir, exist_ok=True)

combined_pdf_path = os.path.join(
    combined_plot_dir,
    "all_motifs_all_celltypes_delta_200bp.pdf"
)

#PdfPages 可以把多个 matplotlib figure 写入同一个 PDF
#后面每运行一次：pdf.savefig(fig)就向总 PDF 添加一页。离开 with 后，PDF 自动关闭和保存
with PdfPages(combined_pdf_path) as pdf:

    # 按 common_motif_dict 中定义的顺序输出 motif
    for motif_name in common_motif_dict.keys():

        motif_results = all_motif_delta.get(motif_name, [])#取出当前 motif 的所有细胞类型中的结果

        if len(motif_results) == 0:
            logger.warning(f"No result found for motif: {motif_name}")
            continue

        # 合并当前 motif 在所有细胞类型中的 -200 至 +200 bp 数值。全部连接成一个一维数组。
        # 使用整体最小值和最大值计算这个 motif 的统一 ylim。
        all_values = np.concatenate([
            result["delta"]
            for result in motif_results
        ])

        # 去除 NaN 和 Inf，防止 ylim 计算报错
        all_values = all_values[np.isfinite(all_values)]

        if all_values.size == 0:
            logger.warning(
                f"No finite delta values found for motif: {motif_name}"
            )
            continue

        # y 轴始终包含 0，便于比较正向和负向 footprint effect
        y_min = min(0.0, float(np.min(all_values)))
        y_max = max(0.0, float(np.max(all_values)))

        # 上下各留 5% 空白，避免曲线贴到图框上
        y_range = y_max - y_min

        if y_range == 0:
            padding = max(abs(y_max) * YLIM_PADDING_FRAC, 1e-6)
        else:
            padding = y_range * YLIM_PADDING_FRAC

        motif_ylim = (
            y_min - padding,
            y_max + padding
        )

        n_celltypes = len(motif_results)
        nrows = int(np.ceil(n_celltypes / COMBINED_NCOLS))

        fig, axes = plt.subplots(
            nrows=nrows,
            ncols=COMBINED_NCOLS,
            figsize=(3.4 * COMBINED_NCOLS, 2.45 * nrows),
            sharex=True, #所有子图共享 x 轴设置
            sharey=True, #所有子图共享 y 轴设置
            squeeze=False #即使只有一行或一列，axes 仍然保持二维数组形式
        )

        axes = axes.ravel() #把二维 axes 数组拉平成一维，方便逐个循环

        for ax, result in zip(axes, motif_results):

            ax.plot(
                result["x"], #横坐标−200～+200 bp
                result["delta"], #纵坐标after profile − before profile
                linewidth=1.0
            )

            # motif 插入中心画一条垂直虚线
            ax.axvline(
                0,
                color="grey",
                linestyle="--",
                linewidth=0.6
            )

            # delta = 0 参考线
            ax.axhline(
                0,
                color="grey",
                linestyle="--",
                linewidth=0.6
            )

            ax.set_xlim(-PLOT_WINDOW, PLOT_WINDOW)
            ax.set_ylim(motif_ylim)
            ax.set_xticks([-200, -100, 0, 100, 200])

            # 去掉 TS_C数字 前缀，并将下划线替换为空格，标题更适合展示
            celltype_label = re.sub(
                r"^TS_C[0-9]+",
                "",
                result["cell_class"]
            ).replace("_", " ")

            # 每个子图代表一个细胞类型/cluster
            ax.set_title(
                f'{result["cluster"]}\n{celltype_label}',
                fontsize=7.5
            )

            #ax.tick_params(
             #   axis="both",
              #  which="major",
               # direction="out",
                #length=3,
                #labelsize=7
            #)

            # 明确打开并设置刻度短线
            ax.tick_params(
                axis="both",
                which="major",
                reset=True,
                bottom=True,
                top=False,
                left=True,
                right=False,
                direction="out",
                length=4,
                width=0.8,
                colors="black",
                labelsize=7
            )

            # 明确指定刻度位置
            ax.xaxis.set_ticks_position("bottom")
            ax.yaxis.set_ticks_position("left")

        # 最后一行可能有多余的空白子图，将其隐藏。这段代码关闭多余子图的：坐标轴；边框；刻度。
        for ax in axes[n_celltypes:]:
            ax.axis("off")

        #整页标题
        fig.suptitle(
            (
                f"{motif_name} marginal footprint across cell types\n"
                f"shared ylim: {motif_ylim[0]:.4g} to {motif_ylim[1]:.4g}" #.4g 表示使用大约 4 位有效数字，数值较小时自动使用科学计数法
            ),
            fontsize=13,
            y=0.995
        )

        #设置全局坐标轴标题，这样不需要在每个小图上重复写坐标轴标题
        fig.supxlabel(
            "Position relative to inserted motif center (bp)",
            fontsize=10
        )

        fig.supylabel(
            "Predicted profile change (After - Before)",
            fontsize=10
        )

        #调整排版。tight_layout() 自动减少子图之间的重叠。rect 表示子图可使用的区域：左边留 3%底部留 3%右边到 100%顶部只用到 96%
        fig.tight_layout(
            rect=[0.03, 0.03, 1.0, 0.96]
        )

        # 保存到总的多页 PDF：每个 motif 一页
        pdf.savefig(
            fig,
            bbox_inches="tight"
        )

        # 同时为每个 motif 单独保存一份 PDF
        safe_motif_name = (
            motif_name
            .replace("/", "_")
            .replace("+", "_plus_")
        )

        motif_pdf_path = os.path.join(
            combined_plot_dir,
            f"{safe_motif_name}_all_celltypes_delta_200bp.pdf"
        )

        fig.savefig(
            motif_pdf_path,
            bbox_inches="tight"
        )

        plt.close(fig)

        logger.critical(
            f"Saved combined plot for {motif_name}; "
            f"shared ylim={motif_ylim}"
        )

logger.critical(
    f"Saved all combined motif plots to: {combined_pdf_path}"
)
logger.critical("All done.")