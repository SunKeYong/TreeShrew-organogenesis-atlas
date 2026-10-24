####################################################################
# define the CA conserved, CA divergent, and species-specific peaks
# liftOver the peak from species A to species B, reback to species A
# Taking human peak to mouse peak as an example
# vice versa. You can process human peak to threeshrew peak, mouse peak to human peak,  peak to threeshrew peak, threeshrew peak to  mouse peak and threeshrew peak to  human peak
# 
# minMatch=0.1
# overlap > 0.5
# 
####################################################################
# First fix chain files
awk '
BEGIN{FS="[ \t]+"; OFS=" "}
FNR==NR {gsub(/\r/,"",$1); gsub(/\r/,"",$2); map[$1]=$2; next}

$1=="chain"{
  if ($3 in map) $3 = map[$3]   # tName
  if ($8 in map) $8 = map[$8]   # qName
}
{print}
' /home/sunkeyong/softeware/ref_tree_shrew/Download_from_NCBI/chr.txt mm10ToGCA_033439345.1.over.chain > mm10ToTreeShrew.over.fixed.chain

awk '
BEGIN{FS="[ \t]+"; OFS=" "}
FNR==NR {gsub(/\r/,"",$1); gsub(/\r/,"",$2); map[$1]=$2; next}

$1=="chain"{
  if ($3 in map) $3 = map[$3]   # tName
  if ($8 in map) $8 = map[$8]   # qName
}
{print}
' /home/sunkeyong/softeware/ref_tree_shrew/Download_from_NCBI/chr.txt GCA_033439345.1ToMm10.over.chain > TreeShrewToMm10.over.chain
#
#
#
####################################################################
# liftover and define the CA conserved, CA divergent, and species-specific peaks 
sort -k1,1 -k2,2n human_keep_Final_bed_tmp.bed > human.sorted.bed
sort -k1,1 -k2,2n Mouse_keep_Final_bed_tmp.bed > mouse.sorted.bed
#
#
LiftOver=/home/sunkeyong/softeware/LiftOver
Chain_files=/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/chain_files
#
#
$LiftOver/liftOver  ./human.sorted.bed   $Chain_files/hg38ToMm10.over.chain  human.h2m.bed      human.h2m.unmapped.bed       -minMatch=0.1 &
$LiftOver/liftOver  ./human.h2m.bed      $Chain_files/mm10ToHg38.over.chain  human.h2m.back.bed human.h2m.back.unmapped.bed  -minMatch=0.1 &
#
awk 'BEGIN{OFS="\t"}
NR==FNR{
    chr[$4]=$1; s[$4]=$2; e[$4]=$3; len[$4]=$3-$2;
    next
}
($4 in len) && ($1 == chr[$4]){
    ov = ((e[$4] < $3 ? e[$4] : $3) - (s[$4] > $2 ? s[$4] : $2));
    if (ov < 0) ov = 0;
    if (ov > maxov[$4]) {
        maxov[$4] = ov;
        bchr[$4] = $1; bs[$4] = $2; be[$4] = $3;
    }
}
END{
    for(id in len){
        frac = (id in maxov ? maxov[id]/len[id] : 0);
        print id, len[id], (id in maxov ? maxov[id] : 0), frac, \
              (id in bchr ? bchr[id] : "."), \
              (id in bs ? bs[id] : -1), \
              (id in be ? be[id] : -1);
    }
}' human.sorted.bed human.h2m.back.bed \
> human.reciprocal_similarity.tsv
#
awk '$4 > 0.5 {print $1}' human.reciprocal_similarity.tsv | sort -u > human.seq_conserved.ids
#
awk 'NR==FNR{ok[$1]=1; next} ($4 in ok)' \
  human.seq_conserved.ids human.h2m.bed \
  | sort -k1,1 -k2,2n \
  > human.h2m.seq_conserved.bed
#
bedtools intersect -u \
  -a human.h2m.seq_conserved.bed \
  -b mouse.sorted.bed \
  > human.h2m.seq_conserved.open_in_mouse.bed
#
cut -f4 human.h2m.seq_conserved.open_in_mouse.bed | sort -u > human.CA_conserved.ids
#
awk 'BEGIN{OFS="\t"}
FILENAME==ARGV[1]{ca[$1]=1; next}
FILENAME==ARGV[2]{seq[$1]=1; next}
{
    cls = (($4 in ca) ? "CA_conserved" : (($4 in seq) ? "CA_divergent" : "Human_specific"));
    print $0, cls
}' human.CA_conserved.ids human.seq_conserved.ids human.sorted.bed \
> human.classified.bed
#
awk 'BEGIN{
    cons=0; div=0; spec=0
}
$5=="CA_conserved"{cons++}
$5=="CA_divergent"{div++}
$5=="Human_specific"{spec++}
END{
    total=cons+div+spec
    print "CA_conserved\t" cons "\t" cons/total
    print "CA_divergent\t" div "\t" div/total
    print "Human_specific\t" spec "\t" spec/total
}' human.classified.bed
#

