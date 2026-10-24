################################
## Process raw fastq files to get final fragment bed files
## fragment bed files were processed by sinto
# see details in github：https://github.com/Lan-lab/SPATAC-seq/
#  Treeshrew GTF is downloaded from from https://www.treeshrewdb.org/
#  Treeshrew Genome is downloaded from from https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_033439345.1/
#  The chromsome of Treeshrew Genome fasta were changed into standard format, such as chr1, chr2
# 
##########
# trim raw fastq files
for sample in SPATACseq-SubLibrary_R1.fq.gz
do
i=`echo $sample|cut -d '_' -f 1`
#
zcat ${i}_R1.fq.gz | fastx_trimmer -Q33 -f 66 -l 116 -i -o | gzip > ${i}_F_1.fastq.gz
zcat ${i}_R2.fq.gz | fastx_trimmer -Q33 -f 51 -l 100 -i -o | gzip > ${i}_F_2.fastq.gz
zcat ${i}_R1.fq.gz | fastx_trimmer -Q33 -f 39 -l 46 -i -o > ${i}_R1F.fastq
zcat ${i}_R2.fq.gz | fastx_trimmer -Q33 -f 24 -l 31 -i -o > ${i}_R1R.fastq
zcat ${i}_R1.fq.gz | fastx_trimmer -Q33 -f 1 -l 8 -i -o > ${i}_R2BC.fastq
zcat ${i}_R2.fq.gz | fastx_trimmer -Q33 -f 1 -l 8 -i -o > ${i}_R3BC.fastq
#
done
#
##
##
##########
# barcode correction
for sample in SPATACseq-SubLibrary_R1.fq.gz
do
i=`echo $sample|cut -d '_' -f 1`
python3 /home/sunkeyong/data_storage/SPATAC/Code_for_demultiplex/ss2_py.py -1  /home/sunkeyong/data_storage/SPATACseqv3_BC/nSPATAC_Round1F_BC_Final.csv -2 ${i}_R1F.fastq -o ${i}_R1F_correction.fastq
python3 /home/sunkeyong/data_storage/SPATAC/Code_for_demultiplex/ss2_py.py -1  /home/sunkeyong/data_storage/SPATACseqv3_BC/nSPATAC_Round1R_BC_Final.csv -2 ${i}_R1R.fastq -o ${i}_R1R_correction.fastq
python3 /home/sunkeyong/data_storage/SPATAC/Code_for_demultiplex/ss2_py.py -1  /home/sunkeyong/data_storage/SPATACseqv3_BC/NewR2.csv -2 ${i}_R2BC.fastq -o ${i}_R2BC_correction.fastq
python3 /home/sunkeyong/data_storage/SPATAC/Code_for_demultiplex/ss2_py.py -1  /home/sunkeyong/data_storage/SPATACseqv3_BC/nRound3num.correction.Final2.csv -2 ${i}_R3BC.fastq -o ${i}_R3BC_correction.fastq
done
#
#
for sample in SPATACseq-SubLibrary_R1.fq.gz
do
i=`echo $sample|cut -d '_' -f 1`
python /home/sunkeyong/data_storage/ZV2/ref/software/BC_process.py -i ${i}_R1F_correction.fastq -i2 ${i}_R1R_correction.fastq -i3 ${i}_R2BC_correction.fastq -i4 ${i}_R3BC_correction.fastq -o ${i}_BC.fastq
done
#
##########
# attach barcode into header
for sample in SPATACseq-SubLibrary_R1.fq.gz
do
i=`echo $sample|cut -d '_' -f 1`
sinto barcode -b 32 --barcode_fastq ${i}_BC.fastq --read1 ${i}_F_1.fastq.gz
sinto barcode -b 32 --barcode_fastq ${i}_BC.fastq --read1 ${i}_F_2.fastq.gz

done
#
#
#
#
#
##########
# mapping
for sample in SPATACseq-SubLibrary_R1.fq.gz
do
i=`echo $sample|cut -d '_' -f 1`
mv ${i}_F_1.barcoded.fastq.gz ${i}_1.fq.gz
mv ${i}_F_2.barcoded.fastq.gz ${i}_2.fq.gz
#
bwa mem -t 22 /home/sunkeyong/softeware/ref_tree_shrew/TreeShrew_scATAC_ref/fasta/genome.fa ${i}_1.fq.gz ${i}_2.fq.gz | samtools view -b - > ${i}.raw.bam
samtools sort -@ 22 ${i}.raw.bam -o ${i}.sort.bam
samtools index ${i}.sort.bam
samtools flagstat  ${i}.sort.bam > ${i}.sort.bam.stat
rm ${i}.raw.bam
done
#
#
##########
# get fragments bed files by sinto 
for sample in SPATACseq-SubLibrary.sort.bam
do
i=`echo $sample|cut -d '.' -f 1`
sinto fragments -b $sample -p 20 -m 10 --use_chrom "[0-9]*" -f ${i}.fragments.bed  --barcode_regex "[^:]*"
grep -F -f /home/sunkeyong/data_storage/SPATACseqv3_BC/BC.combination.new.8.12.96.48.csv ${i}.fragments.bed > ${i}.fragments.filter.bed
awk -v OFS='\t' '{$4=$4"_"gensub(/\..+$/,"",1,FILENAME); print}' ${i}.fragments.filter.bed > ${i}.fragments.filter.SN.bed
sort -k1,1 -k2,2n ${i}.fragments.filter.SN.bed > ${i}.bed
bgzip ${i}.bed
tabix -p bed ${i}.bed.gz
rm ${i}.fragments.bed
rm ${i}.fragments.filter.bed
done
#
