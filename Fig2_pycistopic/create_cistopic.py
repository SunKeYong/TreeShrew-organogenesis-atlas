import os
os.chdir("~/Treeshrew/pycistopic_260418")

import pickle
import gzip 
from scipy import sparse
from scipy.io import mmwrite
import pycisTopic 
pycisTopic.__version__

import scanpy as sc

#创建一个名为"outs"的目录
import os
out_dir = "outs"
os.makedirs(out_dir, exist_ok = True)

#创建一个名为"data"的目录
os.makedirs("data", exist_ok = True)

#cistopic_obj = pickle.load(open("~/Treeshrew/pycistopic/cistopic_obj.pkl", "rb"))
#print(cistopic_obj)


adata = sc.read_h5ad("~/Treeshrew/pycistopic_260418/TS_pseudoCell_KeepPeakMatrix.T260413.h5ad")

# 访问表达矩阵（稀疏矩阵格式）
sparse_matrix = adata.X

sparse_matrix_T= sparse_matrix.T
#print(sparse_matrix_T)

csr_sparse_matrix_T=sparse_matrix_T.tocsr()
#print(csr_sparse_matrix_T)

# 访问观测数据（样本/细胞注释）
#cell_metadata = adata.obs
cell_barcodes = adata.obs.index
cell_barcodes_list=cell_barcodes.tolist()

# 访问变量数据（基因/特征注释）
#peak_metadata = adata.var
region_names = adata.var_names
region_names_list=region_names.tolist()

from pycisTopic.cistopic_class import create_cistopic_object_from_fragments
# Create cisTopic object from sparse fragment matrix.
cistopic_obj = pycisTopic.cistopic_class.create_cistopic_object(
    fragment_matrix=csr_sparse_matrix_T,
    cell_names=cell_barcodes_list,
    region_names=region_names_list)

# 保存对象
import pickle
pickle.dump(
    cistopic_obj,
    open("cistopic_obj.pkl", "wb")
)
print("Finishing creating cistopic object!!!")

# get binary matrix(final_binary_accessibility_matrix.mtx应该是原来从h5ad读取的binary_accessibility_matrix.mtx转置后的矩阵)
binary_matrix = cistopic_obj.binary_matrix
# write to file
mmwrite("data/final_binary_accessibility_matrix.mtx", binary_matrix)
binary_matrix

# get region names and cell barcodes
regions = cistopic_obj.region_names
cells = cistopic_obj.cell_names

# write to files
with open('data/final_region_names.txt', 'w') as f:
    for name in regions:
        f.write(f"{name}\n")

with open('data/final_cell_barcodes.txt', 'w') as f:
    for name in cells:
        f.write(f"{name}\n")


binary_matrix_path='data/final_binary_accessibility_matrix.mtx'
outfile='outs/corpus'
mallet_path='~/study_test/pyCisTopic/Endo_DS_test2/Mallet-202108/bin/mallet'
## create corpus corpus（语料库） 是指将原始数据转换为Mallet可以理解的格式后的文件
os.system(f"pycistopic topic_modeling mallet create_corpus -i {binary_matrix_path} -o {outfile} -m 150 -b {mallet_path}")
print("Creating corpus!!!")

mallet_corpus_path='outs/corpus'
os.makedirs("outs/models", exist_ok = True)
output_prefix='outs/models/tm'
## run topic modeling: -t topic, -p threads, -n iterations, -m memory
os.system(f"pycistopic topic_modeling mallet run -i {mallet_corpus_path} -o {output_prefix} -t 20 50 100 150 190 195 200 205 210 230 250 -p 35 -n 200 -a 50 -A True -e 0.1 -E False -s 555 -m 220 -b {mallet_path} -v")


binary_matrix_path='data/final_binary_accessibility_matrix.mtx'
cells_path='data/final_cell_barcodes.txt'
regions_path='data/final_region_names.txt'
output_prefix='outs/models/tm'

## assemble topic models
os.system(f"pycistopic topic_modeling mallet stats -i {binary_matrix_path} -c {cells_path} -r{regions_path} -o {output_prefix} -t 20 50 100 150 190 195 200 205 210 230 250 -v")

print("finishing building models successfully!!!!!")

#models = pickle.load(open("~/Treeshrew/pycistopic_new/outs/models.pkl", "rb")) #已经有models了，直接接着添加进去
models = []
for file in os.listdir("outs/models"):
   # check wether file is a result from topic modelling, e.g. based on the name
   if file.endswith(".pkl"):
      model = pickle.load(open(os.path.join("outs/models", file), "rb"))
      models.append(model)

pickle.dump(
    models,
    open(os.path.join(out_dir, "models.pkl"), "wb")
)
print("Build models successfully!!!!!")