#pyBigwig environment
import os
import numpy as np
import pandas as pd
import pyBigWig

import subprocess
import tempfile
import pyfaidx

# --------------------------------------------------------------------------1 先从fragment file获得cut-site bigwig files --------------------------------------------------------------------------
# =========================
# 0. 函数
# =========================
def is_gz_file(filepath):
    # https://stackoverflow.com/questions/3703276/how-to-tell-if-a-file-is-gzip-compressed
    with open(filepath, 'rb') as test_f:
        return test_f.read(2) == b'\x1f\x8b'


def fragment_to_tagalign_stream(fragment_file_path):
    """
    Expected format for fragment file: tsv with columns chr, start, end and optionally more columns
    """
    frag_is_gz = is_gz_file(fragment_file_path)
    read_method = "zcat " if frag_is_gz else "cat "
    cmd = read_method + fragment_file_path + """ | awk -v OFS="\\t" '{print $1,$2,$3,1000,0,"+"; print $1,$2,$3,1000,0,"-"}'"""
    p = subprocess.Popen([cmd], stdout=subprocess.PIPE, shell=True)
    return p

def stream_filtered_tagaligns(src_tagaligns_stream, genome_file, out_stream):
    '''
    Given a tagalign subprocess stream and reference genome file, filters
    out any reads in chromosomes not included in the reference. Reads in the
    reference chromosomes are sent to the specified output stream.
    
    Returns:
        Boolean. Indicates whether any reads not in the reference fasta were 
        detected.
    '''
    has_unknown_chroms = False
    with pyfaidx.Fasta(genome_file) as g:
        for line in iter(src_tagaligns_stream.stdout.readline, b''):
            tagalign_chrom = line.decode('utf-8').strip().split('\t')[0]
            if tagalign_chrom in g.keys():
                out_stream.stdin.write(line)
            else:
                has_unknown_chroms = True
        src_tagaligns_stream.stdout.close()
    return has_unknown_chroms

def generate_bigwig_simple(input_fragment_file, output_prefix, genome_fasta_file, bsort, tmpdir, no_st, chrom_sizes_file, plus_shift_delta, minus_shift_delta):
    
    p1 = fragment_to_tagalign_stream(input_fragment_file)

    if tmpdir is None:
        if bsort:
            cmd = """awk -v OFS="\\t" '{{if ($6=="+"){{print $1,$2{0:+},$3,$4,$5,$6}} else if ($6=="-") {{print $1,$2,$3{1:+},$4,$5,$6}}}}' | sort -k1,1 | bedtools genomecov -bg -5 -i stdin -g {2} | bedtools sort -i stdin """.format(plus_shift_delta, minus_shift_delta, chrom_sizes_file)
        else:
            cmd = """awk -v OFS="\\t" '{{if ($6=="+"){{print $1,$2{0:+},$3,$4,$5,$6}} else if ($6=="-") {{print $1,$2,$3{1:+},$4,$5,$6}}}}' | sort -k1,1 | bedtools genomecov -bg -5 -i stdin -g {2} | LC_COLLATE="C" sort -k1,1 -k2,2n """.format(plus_shift_delta, minus_shift_delta, chrom_sizes_file)
    else:
        assert(os.path.isdir(tmpdir)) # tmp dir path does not exsist
        if bsort:
            cmd = """awk -v OFS="\\t" '{{if ($6=="+"){{print $1,$2{0:+},$3,$4,$5,$6}} else if ($6=="-") {{print $1,$2,$3{1:+},$4,$5,$6}}}}' | sort -T {3} -k1,1 | bedtools genomecov -bg -5 -i stdin -g {2} | bedtools sort -i stdin """.format(plus_shift_delta, minus_shift_delta, chrom_sizes_file, tmpdir)
        else:
            cmd = """awk -v OFS="\\t" '{{if ($6=="+"){{print $1,$2{0:+},$3,$4,$5,$6}} else if ($6=="-") {{print $1,$2,$3{1:+},$4,$5,$6}}}}' | sort -T {3} -k1,1 | bedtools genomecov -bg -5 -i stdin -g {2} | LC_COLLATE="C" sort -T {3} -k1,1 -k2,2n """.format(plus_shift_delta, minus_shift_delta, chrom_sizes_file, tmpdir)

    print(cmd)


    tmp_bedgraph = tempfile.NamedTemporaryFile()
    if no_st:
        print("Making BedGraph (Do not filter chromosomes not in reference fasta)")

        with open(tmp_bedgraph.name, 'w') as f:
            p2 = subprocess.Popen([cmd], stdin=p1.stdout, stdout=f, shell=True)
            p1.stdout.close()
            p2.communicate()
    else:
        print("Making BedGraph (Filter chromosomes not in reference fasta)")

        with open(tmp_bedgraph.name, 'w') as f:
            p2 = subprocess.Popen([cmd], stdin=subprocess.PIPE, stdout=f, shell=True)
            stream_filtered_tagaligns(p1, genome_fasta_file, p2)
            p2.communicate()

    print("Making Bigwig")
    subprocess.run(["bedGraphToBigWig", tmp_bedgraph.name, chrom_sizes_file, output_prefix + "_unstranded.bw"])

    tmp_bedgraph.close()

# =========================
# 1. cell type lists
# =========================

hs_immune_cell_list = [
    "C10.1", "C10.2", "C10.3", "C10.5", "C10.6", "C10.7",
    "C10.8", "C10.9", "C10.10", "C10.11", "C10.14",
    "C3.1", "C3.3"
]

TS_immune_cell_list = [
    "C24.1", "C24.3", "C24.4", "C24.8", "C24.10",
    "C24.11", "C24.12", "C24.14", "C24.15", "C24.16",
    "C24.17", "C24.18", "C24.21"
]


# =========================
# 2. input/output paths
# =========================

# human fragment / pseudo-cell BED dir
hs_input_dir = "~/Human_embryo_project_T250614/FileSave/sub_ID_bed_PseudoCel/bed_after_sort"

# tree shrew fragment / pseudo-cell BED dir
TS_input_dir = "~/Treeshrew/FileSave/bed_by_celltype_for_traning/bed_after_sort"

# output bigWig dirs
hs_output_dir = "~/chrombpnet/Hs_TS_immune_chrombpnet/observe_bigwig/human_observe_bigwig"
TS_output_dir = "~/chrombpnet/Hs_TS_immune_chrombpnet/observe_bigwig/treeshrew_observe_bigwig"

os.makedirs(hs_output_dir, exist_ok=True)
os.makedirs(TS_output_dir, exist_ok=True)

# =========================
# 3. genome files
# =========================

hs_genome_fasta_file = "~/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.fa"
hs_chrom_sizes_file = "~/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.chrom.sizes"

# tree shrew genome files
TS_genome_fasta_file = "~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.fa"
TS_chrom_sizes_file = "~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.chrom.sizes"


# =========================
# 4. shared parameters
# =========================

tmpdir = None
bsort = False
no_st = False

#bed文件已经考虑了Tn5 shift
plus_shift_delta = 0
minus_shift_delta = 0

# 如果你想按 ATAC-seq 常见 Tn5 shift，可考虑：
# plus_shift_delta = 4
# minus_shift_delta = -5

'''
# =========================
# 5. human loop
# =========================

for celltype in hs_immune_cell_list:
    input_fragment_file = os.path.join(hs_input_dir, f"{celltype}.sort.bed")
    output_prefix = os.path.join(hs_output_dir, celltype)

    if not os.path.exists(input_fragment_file):
        print(f"[SKIP] Missing input file: {input_fragment_file}")
        continue

    print("=" * 80)
    print(f"[Human] Generating cut-site bigWig for {celltype}")
    print(f"Input : {input_fragment_file}")
    print(f"Output: {output_prefix}_unstranded.bw")

    generate_bigwig_simple(
        input_fragment_file=input_fragment_file,
        output_prefix=output_prefix,
        genome_fasta_file=hs_genome_fasta_file,
        bsort=bsort,
        tmpdir=tmpdir,
        no_st=no_st,
        chrom_sizes_file=hs_chrom_sizes_file,
        plus_shift_delta=plus_shift_delta,
        minus_shift_delta=minus_shift_delta
    )


# =========================
# 6. tree shrew loop
# =========================

for celltype in TS_immune_cell_list:
    input_fragment_file = os.path.join(TS_input_dir, f"{celltype}.sort.bed")
    output_prefix = os.path.join(TS_output_dir, celltype)

    if not os.path.exists(input_fragment_file):
        print(f"[SKIP] Missing input file: {input_fragment_file}")
        continue

    print("=" * 80)
    print(f"[Tree shrew] Generating cut-site bigWig for {celltype}")
    print(f"Input : {input_fragment_file}")
    print(f"Output: {output_prefix}_unstranded.bw")

    generate_bigwig_simple(
        input_fragment_file=input_fragment_file,
        output_prefix=output_prefix,
        genome_fasta_file=TS_genome_fasta_file,
        bsort=bsort,
        tmpdir=tmpdir,
        no_st=no_st,
        chrom_sizes_file=TS_chrom_sizes_file,
        plus_shift_delta=plus_shift_delta,
        minus_shift_delta=minus_shift_delta
    )

'''
# --------------------------------------------------------------------------2 再从cut-site bigwig files计算observed counts ---------------------------------------------------------------
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


hs_bigwig_dir="~/chrombpnet/Hs_TS_immune_chrombpnet/observe_bigwig/human_observe_bigwig"
TS_bigwig_dir="~/chrombpnet/Hs_TS_immune_chrombpnet/observe_bigwig/treeshrew_observe_bigwig"

outdir="~/chrombpnet/Hs_TS_immune_chrombpnet"

# =========================
# 3. shared parameters
# =========================

width=502

# =========================
# 4. human loop
# =========================

#统计人的观测值
#hs_obs_cts_dict={}
hs_obs_total_counts_dict = {}
#hs_obs_logcounts_dict = {}

for hs_immune_cell in hs_immune_cell_list:
    bw_path = os.path.join(hs_bigwig_dir, f"{hs_immune_cell}_unstranded.bw")
    with pyBigWig.open(bw_path) as bw:
        cts = get_cts(hs_peak_df, bw, width=width)
        #hs_obs_cts_dict[hs_immune_cell] = cts
        hs_obs_total_counts_dict[hs_immune_cell] = cts.sum(axis=1)
        #hs_obs_logcounts_dict[hs_immune_cell] = np.log1p(cts.sum(axis=1))


hs_obs_total_df = hs_peak_df[["chr", "start", "end"]].copy()

for celltype, total_counts in hs_obs_total_counts_dict.items():
    hs_obs_total_df[f"{celltype}_obs_counts"] = total_counts
    hs_obs_total_df[f"{celltype}_obs_logcounts"] = np.log1p(total_counts)

hs_obs_total_df.to_csv(os.path.join(outdir,"hs_observed_counts_immune_celltypes.csv"), index=False)

# =========================
# 6. tree shrew loop
# =========================

#统计树鼩的观测值
#TS_obs_cts_dict={}
TS_obs_total_counts_dict = {}
#TS_obs_logcounts_dict = {}

for TS_immune_cell in TS_immune_cell_list:
    bw_path = os.path.join(TS_bigwig_dir, f"{TS_immune_cell}_unstranded.bw")
    with pyBigWig.open(bw_path) as bw:
        cts = get_cts(TS_peak_df, bw, width=width)
        #TS_obs_cts_dict[TS_immune_cell] = cts
        TS_obs_total_counts_dict[TS_immune_cell] = cts.sum(axis=1)
        #TS_obs_logcounts_dict[TS_immune_cell] = np.log1p(cts.sum(axis=1))


TS_obs_total_df = TS_peak_df[["chr", "start", "end"]].copy()

for celltype, total_counts in TS_obs_total_counts_dict.items():
    TS_obs_total_df[f"{celltype}_obs_counts"] = total_counts
    TS_obs_total_df[f"{celltype}_obs_logcounts"] = np.log1p(total_counts)

TS_obs_total_df.to_csv(os.path.join(outdir,"TS_observed_counts_immune_celltypes.csv"), index=False)