#!/bin/bash
source /miniconda3/etc/profile.d/conda.sh 
conda activate chrombpnetv8 

# 获取当前日期和时间（格式：YYYY-MM-DD_HH-MM-SS）
current_time=$(date +"%Y-%m-%d_%H-%M-%S")
# 定义日志文件路径
log_file="~/chrombpnet/chrombpnet_TS_brain_pericytes/log_$current_time.log"


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

echo "start at $current_time"
#----------------树鼩目前没有找到blacklist--代码需要调整一下-genome,chromsizes,split---------------------------------
cd ~/chrombpnet/chrombpnet_TS_brain_pericytes
mkdir -p data
mkdir -p bias_model

data="./data"
bias_model="./bias_model"

bed_file="~/Treeshrew/FileSave/bed_by_celltype_for_traning/bed_after_sort/C14.1.sort.bed"
macs2_peak_bed="~/Treeshrew/FileSave/TreeShrew_peak_1000bp_OBO/modified_peak/C14.1_W1000bp_Peak.modified.bed"

genome_file="~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.fa"
chrom_sizes="~/software/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.chrom.sizes"
#blacklist="~/chrombpnet/chrombpnet_tutorial/data/downloads/blacklist.bed.gz"
#black_tmp_bed="~/chrombpnet/chrombpnet_tutorial/data/downloads/temp.bed"

splits_dir="~/chrombpnet/old/chrombpnet_TS_Mono/fold_split_version2/data/splits"

#preprocessing the data

#Ensure that the peak regions do not intersect with the blacklist regions (extended by 1057 bp on both sides) by running the command below.
#bedtools slop -i ~/chrombpnet/chrombpnet_tutorial/data/downloads/blacklist.bed.gz -g ~/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.chrom.sizes -b 1057 > ~/chrombpnet/chrombpnet_Mono/data/temp.bed
#bedtools intersect -v -a $macs2_peak_bed -b $black_tmp_bed  > $data/peaks_no_blacklist.bed

#Define train, validation and test chromosome splits
#mkdir ~/chrombpnet_tutorial/data/splits
#chrombpnet prep splits -c ~/chrombpnet_tutorial/data/downloads/hg38.chrom.subset.sizes -tcr chr1 chr3 chr6 -vcr chr8 chr20 -op ~/chrombpnet_tutorial/data/splits/fold_0

#Generate non-peaks (background regions)
#for fold in 0; do
 #   echo "Processing fold $fold..."
  #  
    # 设置输入输出路径
   # output_dir="${data}/fold_${fold}"
    #mkdir -p $output_dir
    
    # 运行chrombpnet命令
    #chrombpnet prep nonpeaks \
     #   -g $genome_file \
      #  -p $macs2_peak_bed \
       # -c $chrom_sizes \
        #-fl ${splits_dir}/fold_${fold}.json \
        #-o ${output_dir}/fold_${fold}
    #echo "Completed fold $fold"
    #echo "------------------------"
#done
#training bias model
for fold in 0; do
    echo "正在训练 fold_${fold}..."
    # 设置输入输出路径
    output_dir="${bias_model}/fold_${fold}"
    mkdir -p $output_dir
        for threshold in 0.8 0.9; do
            mkdir -p ${output_dir}/threshold_$(echo "$threshold" | sed 's/\.//')
            chrombpnet bias pipeline \
                -ifrag $bed_file \
                -d "ATAC" \
                -g $genome_file \
                -c $chrom_sizes \
                -p $macs2_peak_bed \
                -n ${data}/fold_${fold}/fold_${fold}_negatives.bed \
                -fl ${splits_dir}/fold_${fold}.json \
                -b ${threshold} \
                -o ${output_dir}/threshold_$(echo "$threshold" | sed 's/\.//') \
                -fp TS_brain_pericytes_1000bp_fold${fold}_$(echo "$threshold" | sed 's/\.//')
             echo "fold_${fold}_threshold_$(echo "$threshold" | sed 's/\.//') 训练完成"
        done
    echo "fold_${fold} 训练完成"
    echo "----------------------------------------"
done
echo "Completing training bias model! Congradulations!"
end_time=$(date +"%Y-%m-%d_%H-%M-%S")
echo "finish at $end_time"