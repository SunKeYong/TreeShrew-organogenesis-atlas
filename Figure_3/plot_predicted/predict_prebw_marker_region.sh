#!/bin/bash
source /miniconda3/etc/profile.d/conda.sh 
conda activate chrombpnetv8

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
base_dir="~/chrombpnet/chrombpnet_use_snakefile/Treeshrew"
genome_file="~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.fa"
chrom_sizes="~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.chrom.sizes"
predict_region_bed="~/chrombpnet/chrombpnet_use_snakefile/Treeshrew/predict_marker_region.window.modified.bed"
predict_region_name="marker_region"

predict_marker_region(){
    local majortype=$1 #声明局部变量。把函数的第一个位置参数赋给它。
    local celltype=$2

    cd ${base_dir}/${majortype}
    models=./output/models
    predict=./output/predict
    data=./data

    celltype_dir="${predict}/${celltype}"
    mkdir -p $celltype_dir

    mkdir -p ${celltype_dir}/pre_bw

    for fold in 0 1 2 3 4; do
        echo "正在预测 ${celltype}_fold_${fold}..."
        output_dir_prebw="${celltype_dir}/pre_bw/fold_${fold}"
        mkdir -p $output_dir_prebw

        chrombpnet pred_bw \
            -bm "${models}/fold_${fold}/${celltype}/models/bias_model_scaled.h5" \
            -cm "${models}/fold_${fold}/${celltype}/models/chrombpnet.h5" \
            -cmb "${models}/fold_${fold}/${celltype}/models/chrombpnet_nobias.h5" \
            -r $predict_region_bed \
            -g $genome_file \
            -c $chrom_sizes \
            -op "$output_dir_prebw/${predict_region_name}_pred_bw" \
            -bs 32

        #echo "${celltype}_fold_${fold}_pred_bw预测完成"
    done
    echo "${celltype} 预测完成"
    echo "----------------------------------------"
    
}

main(){
    predict_marker_region "TS_C1Adrenal_cells" "C1.1"
    predict_marker_region "TS_C2PNS" "C2.5"
    predict_marker_region "TS_C3Cardiomyocytes" "C3.1"
    predict_marker_region "TS_C3Cardiomyocytes" "C3.2"
    predict_marker_region "TS_C4NAG" "C4.1"
    predict_marker_region "TS_C5CNS_neurons" "C5.16"
    predict_marker_region "TS_C6Definitive_erythroid_cells" "C6.1"
    predict_marker_region "TS_C7Endothelial_cells" "C7.1"
    predict_marker_region "TS_C8Epithelial_cells" "C8.1"
    predict_marker_region "TS_C9Gastrointestinal_tract" "C9.1"
    predict_marker_region "TS_C10Hepatocytes" "C10.1"
    predict_marker_region "TS_C11Intermediate_neuronal_progenitors" "C11.2"
    predict_marker_region "TS_C12Kidney_cells" "C12.4"
    predict_marker_region "TS_C13Lung_cells" "C13.1"
    predict_marker_region "TS_C14Mesoderm" "C14.7"
    predict_marker_region "TS_C15Muscle_cells" "C15.1"
    predict_marker_region "TS_C16Neural_crest_PNS_glia" "C16.5"
    predict_marker_region "TS_C17Neural_crest_PNS_neurons" "C17.1"
    predict_marker_region "TS_C18Neuroectoderm_and_glia" "C18.20"
    predict_marker_region "TS_C19Olfactory_sensory_neurons" "C19.1"
    predict_marker_region "TS_C20Oligodendrocytes" "C20.1"
    predict_marker_region "TS_C21Otic_cells" "C21.1"
    predict_marker_region "TS_C22Primitive_erythroid_cells" "C22.1"
    predict_marker_region "TS_C23Retinal_cells" "C23.6"
    predict_marker_region "TS_C24white_blood_cells" "C24.3"
}

main