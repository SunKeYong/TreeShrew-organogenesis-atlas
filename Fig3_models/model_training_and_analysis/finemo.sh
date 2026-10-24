#!/bin/bash
source /miniconda3/etc/profile.d/conda.sh 
conda activate Finemo 

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
echo "Comment:chrombpnet,test,preprocessing"
echo "Processing started at $(date)"
echo "Log file: $log_file"
echo "Script PID: $script_pid"  # 记录PID，便于查找和杀死命令

#设置环境变量,发现chrombpnet可执行脚本被安装到了~/.local/bin里
#PATH：查找可执行文件；LD_LIBRARY_PATH：查找共享库文件（.so 文件）
echo $PATH
export PATH="~/.local/bin:$CONDA_PREFIX/bin:$PATH"
echo $PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
echo $LD_LIBRARY_PATH
#先在镜像中安装 pycairo
#conda install -y -c conda-forge pycairo

#先简单测试一下gpu是否可用

python -c "
import torch
print('PyTorch 版本:', torch.__version__)
print('CUDA 是否可用:', torch.cuda.is_available())
print('GPU 数量:', torch.cuda.device_count())
if torch.cuda.is_available():
    for i in range(torch.cuda.device_count()):
        print(f'GPU {i}: {torch.cuda.get_device_name(i)}')
else:
    print('未检测到 GPU 设备')
"
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

base_dir="~/chrombpnet/chrombpnet_use_snakefile/Treeshrew"
counts_score_h5_dir="${base_dir}/TS_counts_scores_h5"

genome_file="~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.fa"
chrom_sizes="~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.chrom.sizes"

peak_dir="~/Treeshrew/FileSave/TreeShrew_peak_1000bp_OBO/modified_peak"


finemo_function(){
    local majortype=$1 #声明局部变量。把函数的第一个位置参数赋给它。
    local subtype=$2

    mkdir -p "${base_dir}/${majortype}/output/predict"
    mkdir -p "${base_dir}/${majortype}/output/predict/${subtype}"
    finemo_dir="${base_dir}/${majortype}/output/predict/${subtype}/finemo"
    mkdir -p ${finemo_dir}

    tfmodisco_dir="${base_dir}/${majortype}/output/predict/${subtype}/tfmodisco"

    #Preprocessing,transform input contributions and sequences into a compressed .npz file for quick loading
    #treeshrew的peak文件没有经过blacklist过滤
    #使用${peak_dir}/${subtype}_W1000bp_Peak.modified.bed时出现数量不一致，故使用interpreted_regions.bed
    finemo extract-regions-chrombpnet-h5 \
        -c ${counts_score_h5_dir}/average_${subtype}_shaps.counts.h5 \
        -o ${finemo_dir}/${subtype}_counts_finemo \
        -p ${base_dir}/${majortype}/output/models/fold_0/${subtype}/interpret_counts/${subtype}.interpreted_regions.bed \
        -w 1000
    echo "finish preprocessing"

    #Hit Calling(只有这一步用gpu)
    mkdir -p ${finemo_dir}/hit_call
    finemo call-hits -r ${finemo_dir}/${subtype}_counts_finemo.npz \
        -m ${tfmodisco_dir}/${subtype}_counts_modisco_results.h5 \
        -o ${finemo_dir}/hit_call \
        -t 0.3 \
        -l 0.8
    echo "finish hit calling"

    #Output reporting
    mkdir -p ${finemo_dir}/report
    finemo report -r ${finemo_dir}/${subtype}_counts_finemo.npz \
        -H ${finemo_dir}/hit_call \
        -o ${finemo_dir}/report \
        -m ${tfmodisco_dir}/${subtype}_counts_modisco_results.h5 \
        -W 400
    echo "finish reporting"

}


main(){
    finemo_function "TS_C1Adrenal_cells" "C1.1"
    finemo_function "TS_C2PNS" "C2.5"
    #finemo_function "TS_C3Cardiomyocytes" "C3.1"
    #finemo_function "TS_C3Cardiomyocytes" "C3.2"
    finemo_function "TS_C4NAG" "C4.1"
    finemo_function "TS_C5CNS_neurons" "C5.16"
    finemo_function "TS_C6Definitive_erythroid_cells" "C6.1"
    finemo_function "TS_C7Endothelial_cells" "C7.1"
    finemo_function "TS_C8Epithelial_cells" "C8.1"
    finemo_function "TS_C9Gastrointestinal_tract" "C9.1"
    finemo_function "TS_C10Hepatocytes" "C10.1"
    finemo_function "TS_C11Intermediate_neuronal_progenitors" "C11.2"
    finemo_function "TS_C12Kidney_cells" "C12.4"
    finemo_function "TS_C13Lung_cells" "C13.1"
    finemo_function "TS_C14Mesoderm" "C14.7"
    finemo_function "TS_C15Muscle_cells" "C15.1"
    finemo_function "TS_C16Neural_crest_PNS_glia" "C16.5"
    finemo_function "TS_C17Neural_crest_PNS_neurons" "C17.1"
    finemo_function "TS_C18Neuroectoderm_and_glia" "C18.20"
    finemo_function "TS_C19Olfactory_sensory_neurons" "C19.1"
    finemo_function "TS_C20Oligodendrocytes" "C20.1"
    finemo_function "TS_C21Otic_cells" "C21.1"
    finemo_function "TS_C22Primitive_erythroid_cells" "C22.1"
    finemo_function "TS_C23Retinal_cells" "C23.6"
    finemo_function "TS_C24white_blood_cells" "C24.3"
}

main