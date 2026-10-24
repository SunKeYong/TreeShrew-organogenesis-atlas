####################################################################
###### define TE-derived peaks
# the TE annotation files of human and treeshrew are downloaded from UCSC
# the TE annotation files of mouse are downloaded from dfam
# 
# first process human TE annotation file
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/database/rmsk.txt.gz
mv rmsk.txt.gz hg38.rmsk.txt.gz
zcat hg38.rmsk.txt.gz | awk 'BEGIN{OFS="\t"} {print $6,$7,$8,$11,0,$10,$12,$13}' > hg38.rmsk.bed8
awk '$7=="LINE" || $7=="SINE" || $7=="LTR" || $7=="DNA"' hg38.rmsk.bed8 > hg38.TE.only.bed8
#
# see Fig4_TE_TreeShrew_UCSC.sh and Fig4_TE_mm10_Dfam_3.8.sh for process TreeShrew and mouse TE annotation file
#
#
bedtools intersect -a human_keep_Final_bed_tmp.bed -b /home/sunkeyong/softeware/TE_RepeatMasker/hg38_/hg38.TE.only.bed8  -wo \
| awk 'BEGIN{OFS="\t"} { \
    a_len = $3 - $2; \
    b_len = $7 - $6; \
    overlap = $13; \
    print $0, a_len, b_len \
}' > human_keep_Final_bed_tmp.bed.UCSC.TEonly_overlap_final.txt
#
# 
bedtools intersect -a TreeShrew_keep_Final_bed_tmp.bed -b /home/sunkeyong/softeware/TE_RepeatMasker/TreeShrew_UCSC/GCA_033439345.1.TE.annot.only.sort.renameChr.resorted.bed11  -wo \
| awk 'BEGIN{OFS="\t"} { \
    a_len = $3 - $2; \
    b_len = $7 - $6; \
    overlap = $13; \
    print $0, a_len, b_len \
}' > TreeShrew_keep_Final_bed_tmp.UCSC.TEonly_overlap_final.txt


bedtools intersect -a Mouse_keep_Final_bed_tmp.bed -b /home/sunkeyong/softeware/TE_RepeatMasker/mm10_Dfam_3.8/mm10.dfam.with_class.sorted.bed12  -wo \
| awk 'BEGIN{OFS="\t"} { \
    a_len = $3 - $2; \
    b_len = $7 - $6; \
    overlap = $13; \
    print $0, a_len, b_len \
}' > Mouse_keep_Final_bed_tmp.dfam.TEonly_overlap_final.txt


