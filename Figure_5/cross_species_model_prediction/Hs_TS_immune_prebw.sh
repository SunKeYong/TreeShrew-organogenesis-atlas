#!/bin/bash
source /miniconda3/etc/profile.d/conda.sh 
conda activate chrombpnetv8 

# 获取当前日期和时间（格式：YYYY-MM-DD_HH-MM-SS）
current_time=$(date +"%Y-%m-%d_%H-%M-%S")
# 定义日志文件路径
log_file="~/chrombpnet/Hs_TS_immune_chrombpnet/log_$current_time.log"


# 将所有输出（stdout 和 stderr）重定向到日志文件
#/dev/null 是一个特殊的虚拟设备文件，通常被称为 "位桶" 或 "黑洞",丢弃一切写入它的数据
#exec > >(tee -a "$log_file" > /dev/null) 2>&1 #所有输出（正常和错误）都被完美地记录到 $log_file 日志文件中，屏幕上不会显示任何内容
#exec > >(tee -a "$log_file") 2>&1
exec > "$log_file" 2>&1


# 获取当前脚本的 PID
script_pid=$$

# 输出日志信息
echo "Comment:chrombpnet,test, prediction bigwig"
echo "Processing started at $(date)"
echo "Log file: $log_file"
echo "Script PID: $script_pid"  # 记录PID，便于查找和杀死命令

#设置环境变量,发现chrombpnet可执行脚本被安装到了~/.local/bin里
#PATH：查找可执行文件；LD_LIBRARY_PATH：查找共享库文件（.so 文件）
echo $PATH
export PATH="~/.local/bin:$PATH"
echo $PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
echo $LD_LIBRARY_PATH
#先在镜像中安装 pycairo
#conda install -y -c conda-forge pycairo

#先简单测试一下gpu是否可用

python -c "
import tensorflow as tf
print('GPU 设备:', tf.config.list_physical_devices('GPU'))"

#python -c "
#import tensorflow as tf
#print('GPU 设备:', tf.config.list_physical_devices('GPU'))

# 创建一个简单的计算任务在GPU上运行
#with tf.device('/GPU:0'):
#    a = tf.constant([[1.0, 2.0], [3.0, 4.0]])
 #   b = tf.constant([[1.0, 1.0], [0.0, 1.0]])
  #  c = tf.matmul(a, b)
   # print('矩阵乘法结果:', c)
    #print('张量所在设备:', c.device)
#"

#第一次搭建环境时需要做的一次性操作
#cd ~/chrombpnet/
#git clone https://github.com/kundajelab/chrombpnet.git
#pip install -e chrombpnet

cd ~/chrombpnet/Hs_TS_immune_chrombpnet/human_model_predict
hs_C3model_dir=~/chrombpnet/chrombpnet_use_snakefile/Hs_C3Erythroid/output/models
hs_C10model_dir="~/chrombpnet/chrombpnet_use_snakefile/Hs_C10white_blood_cells/output/models"
TS_C24model_dir="~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/TS_C24white_blood_cells/output/models"

#人model预测人的免疫细胞peak
hs_genome_file="~/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.fa"
hs_chrom_sizes="~/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.chrom.sizes"
hs_peak_file="~/chrombpnet/Hs_TS_immune_chrombpnet/human_immune_peak_valid.modified10col.bed" #valid peak(去掉无法扩展成2114bp的区间)
hs_bigwig_dir="~/chrombpnet/Hs_TS_immune_chrombpnet/observe_bigwig/human_observe_bigwig" #应该要用cut-site bigwig

#-bw If bigwig is provided, individual predictions at each of the input regions is output to a .h5 along with metrics comparing with ground truth - pearsonr on count and jsd on profile head.
for fold in 0 1 2 3 4; do
    for celltype in C10.1 C10.2 C10.3 C10.5 C10.6 C10.7 C10.8 C10.9 C10.10 C10.11 C10.14; do
        echo "正在用 ${celltype}_fold_${fold}模型预测human_immune_peak..."
        mkdir -p "pre_bw"
        mkdir -p "pre_bw/fold_${fold}"
        mkdir -p "pre_bw/fold_${fold}/${celltype}"
        output_dir="pre_bw/fold_${fold}/${celltype}/human_immune_peak"
        mkdir -p $output_dir

        chrombpnet pred_bw \
            -bm "${hs_C10model_dir}/fold_${fold}/${celltype}/models/bias_model_scaled.h5" \
            -cm "${hs_C10model_dir}/fold_${fold}/${celltype}/models/chrombpnet.h5" \
            -cmb "${hs_C10model_dir}/fold_${fold}/${celltype}/models/chrombpnet_nobias.h5" \
            -r $hs_peak_file \
            -g $hs_genome_file \
            -c $hs_chrom_sizes \
            -op "${output_dir}/${celltype}_fold_${fold}_human_immune_peak_prebw" \
            -bs 32 \
            -bw "${hs_bigwig_dir}/${celltype}_unstranded.bw"
    done
    echo "fold_${fold} 训练完成"
    echo "----------------------------------------"
done
echo "Prediction completed!congradulations!!!"

for fold in 0 1 2 3 4; do
    for celltype in C3.1 C3.3; do
        echo "正在用 ${celltype}_fold_${fold}模型预测human_immune_peak..."
        mkdir -p "pre_bw"
        mkdir -p "pre_bw/fold_${fold}"
        mkdir -p "pre_bw/fold_${fold}/${celltype}"
        output_dir="pre_bw/fold_${fold}/${celltype}/human_immune_peak"
        mkdir -p $output_dir

        chrombpnet pred_bw \
            -bm "${hs_C3model_dir}/fold_${fold}/${celltype}/models/bias_model_scaled.h5" \
            -cm "${hs_C3model_dir}/fold_${fold}/${celltype}/models/chrombpnet.h5" \
            -cmb "${hs_C3model_dir}/fold_${fold}/${celltype}/models/chrombpnet_nobias.h5" \
            -r $hs_peak_file \
            -g $hs_genome_file \
            -c $hs_chrom_sizes \
            -op "${output_dir}/${celltype}_fold_${fold}_human_immune_peak_prebw" \
            -bs 32 \
            -bw "${hs_bigwig_dir}/${celltype}_unstranded.bw"
    done
    echo "fold_${fold} 训练完成"
    echo "----------------------------------------"
done
echo "Prediction completed!congradulations!!!"

#--------------------------------------------------------------------------------------------
#人model预测树鼩的免疫细胞peak
TS_genome_file="~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.fa"
TS_chrom_sizes="~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.chrom.sizes"
TS_peak_file="~/chrombpnet/Hs_TS_immune_chrombpnet/Treeshrew_immune_peak_valid.modified10col.bed"
TS_bigwig_dir="~/chrombpnet/Hs_TS_immune_chrombpnet/observe_bigwig/treeshrew_observe_bigwig"

#人model celltype 到 树鼩bigwig celltype 的对应关系
declare -A HS_TO_TS_CELLTYPE=(
    ["C10.14"]="C24.1"
    ["C10.2"]="C24.3"
    ["C10.3"]="C24.4"
    ["C3.1"]="C24.8"
    ["C10.5"]="C24.10"
    ["C10.7"]="C24.11"
    ["C10.8"]="C24.12"
    ["C3.3"]="C24.14"
    ["C10.9"]="C24.15"
    ["C10.10"]="C24.16"
    ["C10.6"]="C24.17"
    ["C10.11"]="C24.18"
    ["C10.1"]="C24.21"
)
for fold in 0 1 2 3 4; do
    for celltype in C10.1 C10.2 C10.3 C10.5 C10.6 C10.7 C10.8 C10.9 C10.10 C10.11 C10.14; do
        echo "正在用 ${celltype}_fold_${fold}模型预测treeshrew_immune_peak..."
        TS_celltype="${HS_TO_TS_CELLTYPE[$celltype]}"
        mkdir -p "pre_bw"
        mkdir -p "pre_bw/fold_${fold}"
        mkdir -p "pre_bw/fold_${fold}/${celltype}"
        output_dir="pre_bw/fold_${fold}/${celltype}/treeshrew_immune_peak"
        mkdir -p $output_dir

        chrombpnet pred_bw \
            -bm "${hs_C10model_dir}/fold_${fold}/${celltype}/models/bias_model_scaled.h5" \
            -cm "${hs_C10model_dir}/fold_${fold}/${celltype}/models/chrombpnet.h5" \
            -cmb "${hs_C10model_dir}/fold_${fold}/${celltype}/models/chrombpnet_nobias.h5" \
            -r $TS_peak_file \
            -g $TS_genome_file \
            -c $TS_chrom_sizes \
            -op "${output_dir}/${celltype}_fold_${fold}_treeshrew_immune_peak_prebw" \
            -bs 32 \
            -bw "${TS_bigwig_dir}/${TS_celltype}_unstranded.bw"
    done
    echo "fold_${fold} 训练完成"
    echo "----------------------------------------"
done
echo "Prediction completed!congradulations!!!"

for fold in 0 1 2 3 4; do
    for celltype in C3.1 C3.3; do
        echo "正在用 ${celltype}_fold_${fold}模型预测treeshrew_immune_peak..."
        TS_celltype="${HS_TO_TS_CELLTYPE[$celltype]}"
        mkdir -p "pre_bw"
        mkdir -p "pre_bw/fold_${fold}"
        mkdir -p "pre_bw/fold_${fold}/${celltype}"
        output_dir="pre_bw/fold_${fold}/${celltype}/treeshrew_immune_peak"
        mkdir -p $output_dir

        chrombpnet pred_bw \
            -bm "${hs_C3model_dir}/fold_${fold}/${celltype}/models/bias_model_scaled.h5" \
            -cm "${hs_C3model_dir}/fold_${fold}/${celltype}/models/chrombpnet.h5" \
            -cmb "${hs_C3model_dir}/fold_${fold}/${celltype}/models/chrombpnet_nobias.h5" \
            -r $TS_peak_file \
            -g $TS_genome_file \
            -c $TS_chrom_sizes \
            -op "${output_dir}/${celltype}_fold_${fold}_treeshrew_immune_peak_prebw" \
            -bs 32 \
            -bw "${TS_bigwig_dir}/${TS_celltype}_unstranded.bw"
    done
    echo "fold_${fold} 训练完成"
    echo "----------------------------------------"
done
echo "Prediction completed!congradulations!!!"

# ------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cd ~/chrombpnet/Hs_TS_immune_chrombpnet/treeshrew_model_predict
# 树鼩 model celltype 到 human bigwig celltype 的对应关系
declare -A TS_TO_HS_CELLTYPE=(
    ["C24.1"]="C10.14"
    ["C24.3"]="C10.2"
    ["C24.4"]="C10.3"
    ["C24.8"]="C3.1"
    ["C24.10"]="C10.5"
    ["C24.11"]="C10.7"
    ["C24.12"]="C10.8"
    ["C24.14"]="C3.3"
    ["C24.15"]="C10.9"
    ["C24.16"]="C10.10"
    ["C24.17"]="C10.6"
    ["C24.18"]="C10.11"
    ["C24.21"]="C10.1"
)
#用树鼩model预测人的免疫细胞peak
for fold in 0 1 2 3 4; do
    for celltype in C24.1 C24.3 C24.4 C24.8 C24.10 C24.11 C24.12 C24.14 C24.15 C24.16 C24.17 C24.18 C24.21; do
        echo "正在用 ${celltype}_fold_${fold}模型预测human_immune_peak..."

        hs_celltype="${TS_TO_HS_CELLTYPE[$celltype]}"

        mkdir -p "pre_bw"
        mkdir -p "pre_bw/fold_${fold}"
        mkdir -p "pre_bw/fold_${fold}/${celltype}"
        output_dir="pre_bw/fold_${fold}/${celltype}/human_immune_peak"
        mkdir -p $output_dir

        chrombpnet pred_bw \
            -bm "${TS_C24model_dir}/fold_${fold}/${celltype}/models/bias_model_scaled.h5" \
            -cm "${TS_C24model_dir}/fold_${fold}/${celltype}/models/chrombpnet.h5" \
            -cmb "${TS_C24model_dir}/fold_${fold}/${celltype}/models/chrombpnet_nobias.h5" \
            -r $hs_peak_file \
            -g $hs_genome_file \
            -c $hs_chrom_sizes \
            -op "${output_dir}/${celltype}_fold_${fold}_human_immune_peak_prebw" \
            -bs 32 \
            -bw "${hs_bigwig_dir}/${hs_celltype}_unstranded.bw"
    done
    echo "fold_${fold} 训练完成"
    echo "----------------------------------------"
done
echo "Prediction completed!congradulations!!!"


#----------------------------------------------------------------------------------------
#用树鼩model预测树鼩的免疫细胞peak
for fold in 0 1 2 3 4; do
    for celltype in C24.1 C24.3 C24.4 C24.8 C24.10 C24.11 C24.12 C24.14 C24.15 C24.16 C24.17 C24.18 C24.21; do
        echo "正在用 ${celltype}_fold_${fold}模型预测treeshrew_immune_peak..."
        mkdir -p "pre_bw"
        mkdir -p "pre_bw/fold_${fold}"
        mkdir -p "pre_bw/fold_${fold}/${celltype}"
        output_dir="pre_bw/fold_${fold}/${celltype}/treeshrew_immune_peak"
        mkdir -p $output_dir

        chrombpnet pred_bw \
            -bm "${TS_C24model_dir}/fold_${fold}/${celltype}/models/bias_model_scaled.h5" \
            -cm "${TS_C24model_dir}/fold_${fold}/${celltype}/models/chrombpnet.h5" \
            -cmb "${TS_C24model_dir}/fold_${fold}/${celltype}/models/chrombpnet_nobias.h5" \
            -r $TS_peak_file \
            -g $TS_genome_file \
            -c $TS_chrom_sizes \
            -op "${output_dir}/${celltype}_fold_${fold}_treeshrew_immune_peak_prebw" \
            -bs 32 \
            -bw "${TS_bigwig_dir}/${celltype}_unstranded.bw"
    done
    echo "fold_${fold} 训练完成"
    echo "----------------------------------------"
done
echo "Prediction completed!congradulations!!!"