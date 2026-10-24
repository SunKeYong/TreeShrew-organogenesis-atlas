wget https://www.dfam.org/releases/Dfam_3.8/annotations/mm10/mm10.nrph.hits.gz
wget https://www.dfam.org/releases/Dfam_3.8/families/Dfam_curatedonly.embl.gz

zcat mm10.nrph.hits.gz | head -n 5

zcat Dfam_curatedonly.embl.gz \
| awk '
BEGIN{OFS="\t"}
/^NM[[:space:]]/{
  name=$2
}
/^AC[[:space:]]/{
  acc=$2
  sub(/;$/, "", acc)
}
/^KW[[:space:]]/{
  kw = substr($0, 6)
  sub(/^[[:space:]]+/, "", kw)
  sub(/\.[[:space:]]*$/, "", kw)
}
/^\/\/$/{
  if (acc != "") print acc, name, kw
  acc=""; name=""; kw=""
}' > mm10_dfam_acc_name_kw.tsv

head mm10_dfam_acc_name_kw.tsv


awk '
BEGIN{FS=OFS="\t"}
{
  acc=$1
  name=$2
  kw=$3

  repClass="Other"
  repFamily="Other"

  gsub(/\.$/, "", kw)

  if (kw ~ /^[A-Za-z0-9_-]+\/[A-Za-z0-9_.-]+$/) {
    split(kw, a, "/")
    repClass=a[1]
    repFamily=a[2]
  } else if (kw ~ /^Long terminal repeat of retrovirus-like element/) {
    repClass="LTR"
    n=split(kw, b, /;[ ]*/)
    if (n >= 2 && b[2] != "") repFamily=b[2]
    else repFamily=name
  } else {
    repClass=kw
    repFamily=kw
  }

  print acc, name, repClass, repFamily, kw
}' mm10_dfam_acc_name_kw.tsv > mm10_dfam_acc_to_class.clean.tsv

head mm10_dfam_acc_to_class.clean.tsv

awk 'BEGIN{FS=OFS="\t"} $3=="LINE" || $3=="SINE" || $3=="LTR" || $3=="DNA"' \
mm10_dfam_acc_to_class.clean.tsv > mm10_dfam_acc_to_class.TEonly.tsv
head mm10_dfam_acc_to_class.TEonly.tsv

zcat mm10.nrph.hits.gz | awk 'BEGIN{FS=OFS="\t"} !/^#/ {print}' > mm10.nrph.hits.tsv
head mm10.nrph.hits.tsv

join -t $'\t' -1 2 -2 1 \
  <(sort -k2,2 mm10.nrph.hits.tsv) \
  <(sort -k1,1 mm10_dfam_acc_to_class.TEonly.tsv) \
  > mm10.nrph.with_class.tsv

head mm10.nrph.with_class.tsv

awk '
BEGIN{FS=OFS="\t"}
{
  acc=$1
  chr=$2
  repName=$3
  bits=$4
  strand=$9
  envst=$12
  enven=$13
  kimura=$15
  repClass=$17
  repFamily=$18
  kw=$19

  s=(envst<enven?envst:enven)
  e=(envst>enven?envst:enven)

  print chr, s-1, e, repName, bits, strand, acc, kimura, repClass, repFamily, kw, e-s+1
}' mm10.nrph.with_class.tsv > mm10.dfam.with_class.bed12

head mm10.dfam.with_class.bed12

sort -k1,1 -k2,2n mm10.dfam.with_class.bed12 > mm10.dfam.with_class.sorted.bed12
head mm10.dfam.with_class.sorted.bed12


awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,0,$6,$9,$10}' \
mm10.dfam.with_class.sorted.bed12 > mm10.dfam.TE.only.bed8

awk 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,$4,$5,$6,$9,$10}' \
mm10.dfam.with_class.sorted.bed12 > mm10.dfam.TE.only.bits.bed8

head mm10.dfam.TE.only.bed8
head mm10.dfam.TE.only.bits.bed8