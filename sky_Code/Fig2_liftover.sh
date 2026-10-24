####################################################################
#  Fig2B
# liftover Treeshrew genome to human and mouse genome
# 
# https://hgdownload.soe.ucsc.edu/hubs/GCA/033/439/345/GCA_033439345.1/liftOver/GCA_033439345.1ToHg38.over.chain.gz
# https://hgdownload.soe.ucsc.edu/hubs/GCA/033/439/345/GCA_033439345.1/liftOver/GCA_033439345.1ToMm10.over.chain.gz
# https://hgdownload.gi.ucsc.edu/goldenPath/hg38/liftOver/hg38ToGCA_033439345.1.over.chain.gz
# https://hgdownload.gi.ucsc.edu/goldenPath/mm10/liftOver/mm10ToGCA_033439345.1.over.chain.gz
# ####################################################################
# process random region
sort -k1,1 -k2,2n TS_site.random.bed > TS_site.random.sorted.bed
#
LiftOver=/home/sunkeyong/softeware/LiftOver
Chain_files=/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/chain_files
#
$LiftOver/liftOver ./TS_site.random.sorted.bed       $Chain_files/TS.ToHg38.over.fixed.chain   Random_TS_to_human.Min0.10.bed  TS_to_human.Min0.10.unmapped.bed  -minMatch=0.1 &
$LiftOver/liftOver ./TS_site.random.sorted.bed       $Chain_files/TS.ToHg38.over.fixed.chain   Random_TS_to_human.Min0.25.bed  TS_to_human.Min0.25.unmapped.bed  -minMatch=0.25 &
$LiftOver/liftOver ./TS_site.random.sorted.bed       $Chain_files/TS.ToHg38.over.fixed.chain   Random_TS_to_human.Min0.50.bed  TS_to_human.Min0.50.unmapped.bed  -minMatch=0.50 &
$LiftOver/liftOver ./TS_site.random.sorted.bed       $Chain_files/TS.ToHg38.over.fixed.chain   Random_TS_to_human.Min0.75.bed  TS_to_human.Min0.75.unmapped.bed  -minMatch=0.75 &

$LiftOver/liftOver ./TS_site.random.sorted.bed       $Chain_files/TreeShrewToMm10.over.fixed.chain   Random_TS_to_Mouse.Min0.10.bed  TS_to_Mouse.Min0.10.unmapped.bed  -minMatch=0.1 &
$LiftOver/liftOver ./TS_site.random.sorted.bed       $Chain_files/TreeShrewToMm10.over.fixed.chain   Random_TS_to_Mouse.Min0.25.bed  TS_to_Mouse.Min0.25.unmapped.bed  -minMatch=0.25 &
$LiftOver/liftOver ./TS_site.random.sorted.bed       $Chain_files/TreeShrewToMm10.over.fixed.chain   Random_TS_to_Mouse.Min0.50.bed  TS_to_Mouse.Min0.50.unmapped.bed  -minMatch=0.50 &
$LiftOver/liftOver ./TS_site.random.sorted.bed       $Chain_files/TreeShrewToMm10.over.fixed.chain   Random_TS_to_Mouse.Min0.75.bed  TS_to_Mouse.Min0.75.unmapped.bed  -minMatch=0.75 &
#  
#
# ####################################################################
# process TreeShrew peak
sort -k1,1 -k2,2n TreeShrew_keep_Final_bed_tmp.bed > TreeShrew.sorted.bed
#
LiftOver=/home/sunkeyong/softeware/LiftOver
Chain_files=/home/sunkeyong/Human_em/Cross_species/3_species/liftover_/chain_files
#
$LiftOver/liftOver ./TreeShrew.sorted.bed       $Chain_files/TS.ToHg38.over.fixed.chain   TS_to_human.Min0.10.bed  TS_to_human.Min0.10.unmapped.bed  -minMatch=0.1 &
$LiftOver/liftOver ./TreeShrew.sorted.bed       $Chain_files/TS.ToHg38.over.fixed.chain   TS_to_human.Min0.25.bed  TS_to_human.Min0.25.unmapped.bed  -minMatch=0.25 &
$LiftOver/liftOver ./TreeShrew.sorted.bed       $Chain_files/TS.ToHg38.over.fixed.chain   TS_to_human.Min0.50.bed  TS_to_human.Min0.50.unmapped.bed  -minMatch=0.5 &
$LiftOver/liftOver ./TreeShrew.sorted.bed       $Chain_files/TS.ToHg38.over.fixed.chain   TS_to_human.Min0.75.bed  TS_to_human.Min0.75.unmapped.bed  -minMatch=0.75 &
#
$LiftOver/liftOver ./TreeShrew.sorted.bed       $Chain_files/TreeShrewToMm10.over.fixed.chain   TS_to_Mouse.Min0.10.bed  TS_to_Mouse.Min0.10.unmapped.bed  -minMatch=0.1 &
$LiftOver/liftOver ./TreeShrew.sorted.bed       $Chain_files/TreeShrewToMm10.over.fixed.chain   TS_to_Mouse.Min0.25.bed  TS_to_Mouse.Min0.25.unmapped.bed  -minMatch=0.25 &
$LiftOver/liftOver ./TreeShrew.sorted.bed       $Chain_files/TreeShrewToMm10.over.fixed.chain   TS_to_Mouse.Min0.50.bed  TS_to_Mouse.Min0.50.unmapped.bed  -minMatch=0.5 &
$LiftOver/liftOver ./TreeShrew.sorted.bed       $Chain_files/TreeShrewToMm10.over.fixed.chain   TS_to_Mouse.Min0.75.bed  TS_to_Mouse.Min0.75.unmapped.bed  -minMatch=0.75 &
#
$LiftOver/liftOver ./TS_to_Mouse.Min0.10.bed    $Chain_files/mm10ToTreeShrew.over.fixed.chain   TS_to_Mouse.back.Min0.10.bed  TS_to_Mouse.back.Min0.10.unmapped.bed  -minMatch=0.1 &
$LiftOver/liftOver ./TS_to_Mouse.Min0.25.bed    $Chain_files/mm10ToTreeShrew.over.fixed.chain   TS_to_Mouse.back.Min0.25.bed  TS_to_Mouse.back.Min0.25.unmapped.bed  -minMatch=0.25 &
$LiftOver/liftOver ./TS_to_Mouse.Min0.50.bed    $Chain_files/mm10ToTreeShrew.over.fixed.chain   TS_to_Mouse.back.Min0.50.bed  TS_to_Mouse.back.Min0.50.unmapped.bed  -minMatch=0.5 &
$LiftOver/liftOver ./TS_to_Mouse.Min0.75.bed    $Chain_files/mm10ToTreeShrew.over.fixed.chain   TS_to_Mouse.back.Min0.75.bed  TS_to_Mouse.back.Min0.75.unmapped.bed  -minMatch=0.75 &
#
$LiftOver/liftOver ./TS_to_human.Min0.10.bed    $Chain_files/hg38.To.TS.over.fixed.chain   TS_to_human.back.Min0.10.bed  TS_to_human.back.Min0.10.unmapped.bed  -minMatch=0.1 &
$LiftOver/liftOver ./TS_to_human.Min0.25.bed    $Chain_files/hg38.To.TS.over.fixed.chain   TS_to_human.back.Min0.25.bed  TS_to_human.back.Min0.25.unmapped.bed  -minMatch=0.25 &
$LiftOver/liftOver ./TS_to_human.Min0.50.bed    $Chain_files/hg38.To.TS.over.fixed.chain   TS_to_human.back.Min0.50.bed  TS_to_human.back.Min0.50.unmapped.bed  -minMatch=0.5 &
$LiftOver/liftOver ./TS_to_human.Min0.75.bed    $Chain_files/hg38.To.TS.over.fixed.chain   TS_to_human.back.Min0.75.bed  TS_to_human.back.Min0.75.unmapped.bed  -minMatch=0.75 &

