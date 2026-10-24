wget https://hgdownload.soe.ucsc.edu/hubs/GCA/033/439/345/GCA_033439345.1/GCA_033439345.1.repeatMasker.out.gz
wget https://hgdownload.soe.ucsc.edu/hubs/GCA/033/439/345/GCA_033439345.1/GCA_033439345.1.repeatMasker.version.txt
#

zcat GCA_033439345.1.repeatMasker.out.gz \
| awk 'BEGIN{OFS="\t"}
$1 ~ /^[0-9]+$/ {
    chrom  = $5
    start  = $6 - 1
    end    = $7
    strand = ($9 == "C" ? "-" : "+")
    repName = $10
    classfam = $11

    repClass = classfam
    repFamily = classfam
    if (index(classfam, "/") > 0) {
        repClass  = substr(classfam, 1, index(classfam, "/") - 1)
        repFamily = substr(classfam, index(classfam, "/") + 1)
    }

    print chrom, start, end, repName, 0, strand, repClass, repFamily
}' > GCA_033439345.1.repeatMasker.bed8

awk 'BEGIN{FS=OFS="\t"} $7=="LINE" || $7=="SINE" || $7=="LTR" || $7=="DNA"' \
GCA_033439345.1.repeatMasker.bed8 > GCA_033439345.1.TE.only.bed8

sort -k1,1 -k2,2n GCA_033439345.1.TE.only.bed8 > GCA_033439345.1.TE.only.sorted.bed8

head GCA_033439345.1.TE.only.sorted.bed8


zcat GCA_033439345.1.repeatMasker.out.gz \
| awk 'BEGIN{OFS="\t"}
$1 ~ /^[0-9]+$/ {
    chrom  = $5
    start  = $6 - 1
    end    = $7
    strand = ($9 == "C" ? "-" : "+")
    repName = $10
    classfam = $11

    repClass = classfam
    repFamily = classfam
    if (index(classfam, "/") > 0) {
        repClass  = substr(classfam, 1, index(classfam, "/") - 1)
        repFamily = substr(classfam, index(classfam, "/") + 1)
    }

    print chrom, start, end, repName, $1, strand, repClass, repFamily, $2, $3, $4
}' > GCA_033439345.1.TE.annot.bed11

awk 'BEGIN{FS=OFS="\t"} $7=="LINE" || $7=="SINE" || $7=="LTR" || $7=="DNA"' \
GCA_033439345.1.TE.annot.bed11> GCA_033439345.1.TE.annot.only.bed11

sort -k1,1 -k2,2n GCA_033439345.1.TE.annot.only.bed11 > GCA_033439345.1.TE.annot.only.sort.bed11

head GCA_033439345.1.TE.annot.only.sort.bed11

awk 'BEGIN{FS=OFS="\t"}
NR==FNR {map[$1]=$2; next}
{
  if ($1 in map) $1=map[$1]
  print
}' /home/sunkeyong/softeware/ref_tree_shrew/Download_from_NCBI/chr.txt \
GCA_033439345.1.TE.annot.only.sort.bed11 \
> GCA_033439345.1.TE.annot.only.sort.renameChr.bed11

sort -k1,1 -k2,2n GCA_033439345.1.TE.annot.only.sort.renameChr.bed11 \
> GCA_033439345.1.TE.annot.only.sort.renameChr.resorted.bed11

awk 'BEGIN{FS=OFS="\t"} $1 ~ /^chr/' GCA_033439345.1.TE.annot.only.sort.renameChr.bed11 > GCA_033439345.1.TE.annot.only.sort.renameChr.chrOnly.bed11

head GCA_033439345.1.TE.annot.only.sort.renameChr.resorted.bed11
# 