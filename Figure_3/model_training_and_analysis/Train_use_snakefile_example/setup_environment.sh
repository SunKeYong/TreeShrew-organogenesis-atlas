#!/bin/bash
#source /miniconda3/etc/profile.d/conda.sh
#module load cuda/11.2.2 
#module load cudnn/8.1.0.77

eval "$(conda shell.bash hook)"
conda activate chrombpnetv8

#设置环境变量,发现chrombpnet可执行脚本被安装到了~/.local/bin里
#PATH：查找可执行文件；LD_LIBRARY_PATH：查找共享库文件（.so 文件）
#echo $PATH
export PATH="~/.local/bin:$CONDA_PREFIX/bin:$PATH"
#echo $PATH
export LD_LIBRARY_PATH=$CONDA_PREFIX/lib
#echo $LD_LIBRARY_PATH
#export TF_GPU_ALLOCATOR=cuda_malloc_async


#module load system cairo
#module load system pango

#echo "Live"

#特殊参数，所有命令行参数作为一个数组。将传递给脚本的所有参数作为命令执行
"$@"