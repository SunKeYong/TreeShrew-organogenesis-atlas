
########
 # mkdir TS_Human_conserved_min010
#
#
#


cp ../human.sorted.bed .
cp ../TS_to_human.back.Min0.10.bed .
cp ../TS_to_human.Min0.10.bed .
cp ../TreeShrew.sorted.bed .



mkdir -p TS_Human_conserved_min010

# 先排序
sort -k1,1 -k2,2n TreeShrew.sorted.bed \
  > TS_Human_conserved_min010/TreeShrew.sorted.bed

sort -k1,1 -k2,2n TS_to_human.Min0.10.bed \
  > TS_Human_conserved_min010/TS_to_human.Min0.10.sorted.bed

sort -k1,1 -k2,2n TS_to_human.back.Min0.10.bed \
  > TS_Human_conserved_min010/TS_to_human.back.Min0.10.sorted.bed

sort -k1,1 -k2,2n human.sorted.bed \
  > TS_Human_conserved_min010/human.sorted.bed
#
 
# 
# 1. 然后筛选 reciprocal 一致的树鼩 peak：
  bedtools intersect \
  -a TS_Human_conserved_min010/TreeShrew.sorted.bed \
  -b TS_Human_conserved_min010/TS_to_human.back.Min0.10.sorted.bed \
  -wa -wb -f 0.5 -r \
| awk 'BEGIN{OFS="\t"} $4==$8 {print $1,$2,$3,$4}' \
| sort -u -k1,1 -k2,2n \
> TS_Human_conserved_min010/TreeShrew.conserved_sequence.reciprocal.Min0.10.bed
#
# 2. 得到这些树鼩 peak 对应的人同源坐标
cut -f4 TS_Human_conserved_min010/TreeShrew.conserved_sequence.reciprocal.Min0.10.bed \
  > TS_Human_conserved_min010/TreeShrew.conserved_sequence.ids

awk 'NR==FNR{keep[$1]=1; next} keep[$4]' \
  TS_Human_conserved_min010/TreeShrew.conserved_sequence.ids \
  TS_Human_conserved_min010/TS_to_human.Min0.10.sorted.bed \
| sort -k1,1 -k2,2n \
> TS_Human_conserved_min010/Human.ortholog_regions.from_TreeShrew.conserved_sequence.Min0.10.bed
#
# 3. 和 human.sorted.bed 求交，得到真正 conserved activity peaks
bedtools intersect \
  -a TS_Human_conserved_min010/Human.ortholog_regions.from_TreeShrew.conserved_sequence.Min0.10.bed \
  -b TS_Human_conserved_min010/human.sorted.bed \
  -wo -f 0.002 -r \
> TS_Human_conserved_min010/TS_to_Human.conserved_activity.allPairs.Min0.10.recip10.tsv

# 如果一个树鼩 peak 对应多个人 peak，可以保留 overlap 最大的那个：
sort -k4,4 -k9,9nr \
  TS_Human_conserved_min010/TS_to_Human.conserved_activity.allPairs.Min0.10.recip10.tsv \
| awk 'BEGIN{OFS="\t"} !seen[$4]++ {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' \
> TS_Human_conserved_min010/TS_to_Human.conserved_activity.bestPairs.Min0.10.recip10.tsv


# 4. 输出树鼩 conserved peak list
cut -f4 TS_Human_conserved_min010/TS_to_Human.conserved_activity.bestPairs.Min0.10.recip10.tsv \
  > TS_Human_conserved_min010/TreeShrew.conserved_activity.ids

awk 'NR==FNR{keep[$1]=1; next} keep[$4]' \
  TS_Human_conserved_min010/TreeShrew.conserved_activity.ids \
  TS_Human_conserved_min010/TreeShrew.sorted.bed \
| sort -k1,1 -k2,2n \
> TS_Human_conserved_min010/TreeShrew.conserved_activity.peak_list.Min0.10.recip10.bed

# 5. 输出对应的人同源 peak list
cut -f5-8 TS_Human_conserved_min010/TS_to_Human.conserved_activity.bestPairs.Min0.10.recip10.tsv \
| sort -u -k1,1 -k2,2n \
> TS_Human_conserved_min010/Human.conserved_activity.peak_list.matched_to_TreeShrew.Min0.10.recip10.bed

#####
#  树鼩中有可靠人类同源序列的 peaks。
TreeShrew.conserved_sequence.reciprocal.Min0.10.bed

# 这些树鼩 peaks liftOver 到人基因组后的同源坐标。
Human.ortholog_regions.from_TreeShrew.conserved_sequence.Min0.10.bed

# 树鼩 conserved accessible peaks，也就是对应位置在人里也是 peak 的树鼩 peaks。
TreeShrew.conserved_activity.peak_list.Min0.10.recip10.bed

# 对应的人 peaks。
Human.conserved_activity.peak_list.matched_to_TreeShrew.Min0.10.recip10.bed


# 其中最重要的对应关系表是, 它记录了每个树鼩 conserved peak 对应到哪个 human peak。这个文件后面做矩阵提取、画 conserved / divergent / tree-shrew-specific 分类都很有用。
TS_to_Human.conserved_activity.bestPairs.Min0.10.recip10.tsv


# 均是 534931
wc -l TreeShrew.conserved_activity.peak_list.Min0.10.recip10.bed
wc -l Human.conserved_activity.peak_list.matched_to_TreeShrew.Min0.10.recip10.bed
wc -l TS_to_Human.conserved_activity.bestPairs.Min0.10.recip10.tsv

head TS_to_Human.conserved_activity.bestPairs.Min0.10.recip10.tsv
