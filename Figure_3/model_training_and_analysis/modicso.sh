#!/bin/bash
source /miniconda3/etc/profile.d/conda.sh 
conda activate modisco-lite

# 获取当前日期和时间（格式：YYYY-MM-DD_HH-MM-SS）
current_time=$(date +"%Y-%m-%d_%H-%M-%S")
# 定义日志文件路径
log_file="~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/log_$current_time.log"


# 将所有输出（stdout 和 stderr）重定向到日志文件
#/dev/null 是一个特殊的虚拟设备文件，通常被称为 "位桶" 或 "黑洞",丢弃一切写入它的数据
#exec > >(tee -a "$log_file" > /dev/null) 2>&1 #所有输出（正常和错误）都被完美地记录到 $log_file 日志文件中，屏幕上不会显示任何内容
#exec > >(tee -a "$log_file") 2>&1
exec > "$log_file" 2>&1


# 获取当前脚本的 PID
script_pid=$$

# 输出日志信息
echo "Comment:chrombpnet,test,finding motifs"
echo "Processing started at $(date)"
echo "Log file: $log_file"
echo "Script PID: $script_pid"  # 记录PID，便于查找和杀死命令

#设置环境变量,发现chrombpnet可执行脚本被安装到了~/.local/bin里
#PATH：查找可执行文件；LD_LIBRARY_PATH：查找共享库文件（.so 文件）
#echo $PATH
#export PATH="~/.local/bin:$PATH"
export PATH="$CONDA_PREFIX/bin:$PATH"
echo $PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
echo $LD_LIBRARY_PATH

base_dir="~/chrombpnet/chrombpnet_use_snakefile/Treeshrew"
counts_score_h5_dir="${base_dir}/TS_counts_scores_h5"

genome_file="~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.fa"
chrom_sizes="~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.chrom.sizes"

meme_file="~/chrombpnet/chrombpnet_tutorial/data/meme/motifs.meme.txt"



modisco_function(){
    local majortype=$1 #声明局部变量。把函数的第一个位置参数赋给它。
    local subtype=$2

    mkdir -p "${base_dir}/${majortype}/output/predict"
    mkdir -p "${base_dir}/${majortype}/output/predict/${subtype}"
    tfmodisco_dir="${base_dir}/${majortype}/output/predict/${subtype}/tfmodisco"
    mkdir -p ${tfmodisco_dir}

    modisco motifs -i ${counts_score_h5_dir}/average_${subtype}_shaps.counts.h5 -n 1000000 -o ${tfmodisco_dir}/${subtype}_counts_modisco_results.h5
    echo "complete finding motifs !congradulations!!!"


    mkdir -p ${tfmodisco_dir}/report
    modisco report -i ${tfmodisco_dir}/${subtype}_counts_modisco_results.h5 -o ${tfmodisco_dir}/report/ -s ./ -m $meme_file
    echo "reports completed!congradulations!!!"


    mkdir -p ${tfmodisco_dir}/output_memedb
    modisco meme -i ${tfmodisco_dir}/${subtype}_counts_modisco_results.h5 -t PFM -o ${tfmodisco_dir}/output_memedb/output_counts_memedb.meme
    echo "reports completed!congradulations!!!"
}
#C4 C6 C10没跑完
main(){
    #modisco_function "TS_C1Adrenal_cells" "C1.1"
    #modisco_function "TS_C2PNS" "C2.5"
    #modisco_function "TS_C3Cardiomyocytes" "C3.1"
    #modisco_function "TS_C3Cardiomyocytes" "C3.2"
    modisco_function "TS_C4NAG" "C4.1"
    #modisco_function "TS_C5CNS_neurons" "C5.16"
    modisco_function "TS_C6Definitive_erythroid_cells" "C6.1"
    #modisco_function "TS_C7Endothelial_cells" "C7.1"
    #modisco_function "TS_C8Epithelial_cells" "C8.1"
    #modisco_function "TS_C9Gastrointestinal_tract" "C9.1"
    modisco_function "TS_C10Hepatocytes" "C10.1"
    #modisco_function "TS_C11Intermediate_neuronal_progenitors" "C11.2"
    #modisco_function "TS_C12Kidney_cells" "C12.4"
    #modisco_function "TS_C13Lung_cells" "C13.1"
    #modisco_function "TS_C14Mesoderm" "C14.7"
    #modisco_function "TS_C15Muscle_cells" "C15.1"
    #modisco_function "TS_C16Neural_crest_PNS_glia" "C16.5"
    #modisco_function "TS_C17Neural_crest_PNS_neurons" "C17.1"
    #modisco_function "TS_C18Neuroectoderm_and_glia" "C18.20"
    #modisco_function "TS_C19Olfactory_sensory_neurons" "C19.1"
    #modisco_function "TS_C20Oligodendrocytes" "C20.1"
    #modisco_function "TS_C21Otic_cells" "C21.1"
    #modisco_function "TS_C22Primitive_erythroid_cells" "C22.1"
    #modisco_function "TS_C23Retinal_cells" "C23.6"
    #modisco_function "TS_C24white_blood_cells" "C24.3"
}

main