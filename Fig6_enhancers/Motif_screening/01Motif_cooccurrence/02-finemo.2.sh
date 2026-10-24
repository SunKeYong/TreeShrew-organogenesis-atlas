#!/bin/bash
source /miniconda3/etc/profile.d/conda.sh 
conda activate Finemo 

# 获取当前日期和时间（格式：YYYY-MM-DD_HH-MM-SS）
current_time=$(date +"%Y-%m-%d_%H-%M-%S")
# 定义日志文件路径
log_file="~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/log_$current_time.log"


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

cd ~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp
mkdir -p data
mkdir -p chrombpnet_model_predicting
mkdir -p chrombpnet_model_training

mkdir -p ./chrombpnet_model_predicting/contribs_bw
mkdir -p ./chrombpnet_model_predicting/pre_bw
mkdir -p ./chrombpnet_model_predicting/tfmodisco
mkdir -p ./chrombpnet_model_predicting/finemo

data="./data"
chrombpnet_model_predicting="./chrombpnet_model_predicting"
chrombpnet_model_training="./chrombpnet_model_training"

bed_file="~/chrombpnet/chromnpnet_cardiomyocytes/macs2/hs_cardiomyocytes.sort.bed"
macs2_peak_bed="~/chrombpnet/chromnpnet_cardiomyocytes/peak_from_ArchR_1000bp/Hs_cardiomyocytes_1000bp_peak.bed"

genome_file="~/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.fa"
chrom_sizes="~/chrombpnet/chrombpnet_tutorial/data/downloads/hg38.chrom.sizes"
blacklist="~/chrombpnet/chrombpnet_tutorial/data/downloads/blacklist.bed.gz"
black_tmp_bed="~/chrombpnet/chrombpnet_tutorial/data/downloads/temp.bed"

splits_dir="~/chrombpnet/chrombpnet_tutorial/data/splits"

bias_model="~/chrombpnet/chrombpnet_fibroblast/peak_from_ArchR_1000bp/bias_model/fold_0/threshold_05/models/hs_fibroblast_1000bp_fold0_05_bias.h5"

predict_region="heart_tissue_VISTA"

#Preprocessing,transform input contributions and sequences into a compressed .npz file for quick loading
#finemo extract-regions-chrombpnet-h5
#-c/--h5s: One or more H5 files of contribution scores, with paths delimited by whitespace. Scores are averaged across files
#-w/--region-width: The width of the input region centered around each peak summit. Default is 1000
#-p/--peaks: A peak regions file in ENCODE NarrowPeak format. (option) 用contribs_bw过滤后的peak，不然数目不一致
mkdir -p ${chrombpnet_model_predicting}/finemo/${predict_region}
output_dir="${chrombpnet_model_predicting}/finemo/${predict_region}/counts"
mkdir -p $output_dir
finemo extract-regions-chrombpnet-h5 \
    -c ${chrombpnet_model_predicting}/contribs_bw/average_C6.2_heart_tissue_VISTA_shaps.counts.h5 \
    -o ${output_dir}/${predict_region}_counts_scores \
    -p "${data}/region_for_predict/vista/heart/${predict_region}.modified.bed" \
    -w 1000
echo "finish preprocessing"

#Hit Calling(只有这一步用gpu)
#finemo call-hits
#-r/--regions: A .npz file of input sequences, contributions, and coordinates. Created with a finemo extract-regions-* command.
#-m/--modisco-h5: A tfmodisco-lite output H5 file of motif patterns.
#-t/--cwm-trim-threshold: The threshold to determine motif start and end positions within the full CWMs. Default is 0.3.
#-l/--global-lambda: The L1 regularization weight determining the sparsity of hits. Default is 0.7.
mkdir -p $output_dir/hit_call
finemo call-hits -r ${output_dir}/${predict_region}_counts_scores.npz \
    -m ${chrombpnet_model_predicting}/tfmodisco/counts/peaks_modisco_results.h5 \
    -o $output_dir/hit_call \
    -t 0.3 \
    -l 0.8
echo "finish hit calling"

#Output reporting
#finemo report
#-W/--modisco-region-width: The width of the region around each peak summit used by tfmodisco-lite. Default is 400.
mkdir -p $output_dir/report
finemo report -r ${output_dir}/${predict_region}_counts_scores.npz \
    -H $output_dir/hit_call \
    -o $output_dir/report \
    -m ${chrombpnet_model_predicting}/tfmodisco/counts/peaks_modisco_results.h5 \
    -W 400
echo "finish reporting"