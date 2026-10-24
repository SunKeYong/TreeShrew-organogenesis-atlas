############################################# 
# calculate the topic cCREs enrichment in representative cell types
# by deeptools
# 
cd deeptools_enrichement
# 
cp ../BigWig_files/C1.1.bw .
cp ../BigWig_files/C3.2.bw . 
cp ../BigWig_files/C4.1.bw . 
cp ../BigWig_files/C6.1.bw . 
cp ../BigWig_files/C7.1.bw . 
cp ../BigWig_files/C9.1.bw . 
cp ../BigWig_files/C10.1.bw . 
cp ../BigWig_files/C11.2.bw . 
cp ../BigWig_files/C13.1.bw . 
cp ../BigWig_files/C14.7.bw . 
cp ../BigWig_files/C15.1.bw . 
cp ../BigWig_files/C16.5.bw . 
cp ../BigWig_files/C17.1.bw .
cp ../BigWig_files/C19.1.bw . 
cp ../BigWig_files/C23.6.bw . 
cp ../BigWig_files/C24.3.bw . 
#
#
# 
#
# calculate
computeMatrix reference-point  --referencePoint center  -p 16 -b 1000 -a 1000 -R  ../Topics_otsu/topic_62.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C1.1.deeptools.gz --outFileSortedRegions C1.1.deeptools.bed
computeMatrix reference-point  --referencePoint center  -p 4 -b 1000 -a 1000 -R  ../Topics_otsu/topic_37.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C3.2.deeptools.gz --outFileSortedRegions C3.2.deeptools.bed
computeMatrix reference-point  --referencePoint center  -p 4 -b 1000 -a 1000 -R  ../Topics_otsu/topic_85.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C4.1.deeptools.gz --outFileSortedRegions C4.1.deeptools.bed 
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_106.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C6.1.deeptools.gz --outFileSortedRegions C6.1.deeptools.bed 
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_149.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C7.1.deeptools.gz --outFileSortedRegions C7.1.deeptools.bed
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_107.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C9.1.deeptools.gz --outFileSortedRegions C9.1.deeptools.bed
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_134.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C10.1.deeptools.gz --outFileSortedRegions C10.1.deeptools.bed 
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_31.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C11.2.deeptools.gz --outFileSortedRegions C11.2.deeptools.bed
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_167.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C13.1.deeptools.gz --outFileSortedRegions C13.1.deeptools.bed 
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_194.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C14.7.deeptools.gz --outFileSortedRegions C14.7.deeptools.bed
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_64.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C15.1.deeptools.gz --outFileSortedRegions C15.1.deeptools.bed
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_126.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C16.5.deeptools.gz --outFileSortedRegions C16.5.deeptools.bed
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_182.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C17.1.deeptools.gz --outFileSortedRegions C17.1.deeptools.bed
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_196.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C19.1.deeptools.gz --outFileSortedRegions C19.1.deeptools.bed 
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_157.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C23.6.deeptools.gz --outFileSortedRegions C23.6.deeptools.bed
computeMatrix reference-point  --referencePoint center  -p 8 -b 1000 -a 1000 -R  ../Topics_otsu/topic_138.bed -S C1.1.bw C3.2.bw C4.1.bw C6.1.bw C7.1.bw C9.1.bw C10.1.bw C11.2.bw C13.1.bw C14.7.bw C15.1.bw C16.5.bw C17.1.bw C19.1.bw C23.6.bw C24.3.bw   --skipZeros --binSize 40  -o C24.3.deeptools.gz --outFileSortedRegions C24.3.deeptools.bed
#
#
#
#
#
# plot 
plotHeatmap -m C1.1.deeptools.gz  -out C1.1.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only'  
plotHeatmap -m C3.2.deeptools.gz  -out C3.2.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only' 
plotHeatmap -m C4.1.deeptools.gz  -out C4.1.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only' 
plotHeatmap -m C6.1.deeptools.gz  -out C6.1.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only'   
plotHeatmap -m C7.1.deeptools.gz  -out C7.1.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only'  
plotHeatmap -m C9.1.deeptools.gz  -out C9.1.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only'  
plotHeatmap -m C10.1.deeptools.gz  -out C10.1.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only'  
plotHeatmap -m C11.2.deeptools.gz  -out C11.2.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only' 
plotHeatmap -m C13.1.deeptools.gz  -out C13.1.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only' 
plotHeatmap -m C14.7.deeptools.gz  -out C14.7.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only' 
plotHeatmap -m C15.1.deeptools.gz  -out C15.1.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only'  
plotHeatmap -m C16.5.deeptools.gz  -out C16.5.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only'  
plotHeatmap -m C17.1.deeptools.gz  -out C17.1.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only' 
plotHeatmap -m C19.1.deeptools.gz  -out C19.1.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only' 
plotHeatmap -m C23.6.deeptools.gz  -out C23.6.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only' 
plotHeatmap -m C24.3.deeptools.gz  -out C24.3.Heatmap.pdf   --colorList "#f0f0f0",'#F7FBFF','#C6DBEF','#9ECAE1','#6BAED6','#8292C6','#2171B5','#084594'  --dpi 300   --zMin 10  --missingDataColor white  --whatToShow 'heatmap only'  
#
#
#
# 
#############################################
############################################# 
# liftover cCREs of each topic in treeshrew genome  to  human genome for performing GREAT analysis
# 
cd GREAT_analysis
### 
cp  ../Topics_otsu/topic_62.bed S1.topic_62.bed
cp  ../Topics_otsu/topic_37.bed S2.topic_37.bed
cp  ../Topics_otsu/topic_85.bed S3.topic_85.bed
cp  ../Topics_otsu/topic_106.bed S4.topic_106.bed
cp  ../Topics_otsu/topic_149.bed S5.topic_149.bed
cp  ../Topics_otsu/topic_107.bed S6.topic_107.bed
cp  ../Topics_otsu/topic_134.bed S7.topic_134.bed
cp  ../Topics_otsu/topic_31.bed S8.topic_31.bed
cp  ../Topics_otsu/topic_167.bed S9.topic_167.bed
cp  ../Topics_otsu/topic_194.bed S10.topic_194.bed
cp  ../Topics_otsu/topic_64.bed S11.topic_64.bed
cp  ../Topics_otsu/topic_126.bed S12.topic_126.bed
cp  ../Topics_otsu/topic_182.bed S13.topic_82.bed
cp  ../Topics_otsu/topic_196.bed S14.topic_96.bed
cp  ../Topics_otsu/topic_157.bed S15.topic_157.bed
cp  ../Topics_otsu/topic_138.bed S16.topic_138.bed
#
cat S*bed > cat.bed
sort -k1,1 -k2,2n -k3,3n  cat.bed  | uniq > cat.unique.bed
#
#
chain_files=/home/sunkeyong/softeware/LiftOver/Chain_files/TS.ToHg38.over.fixed.chain
/home/sunkeyong/softeware/LiftOver/liftOver  ./S1.topic_62.bed $chain_files S1_topic_62_to_hg38.bed  S1_topic_62_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S2.topic_37.bed $chain_files S2_topic_37_to_hg38.bed  S2_topic_37_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S3.topic_85.bed $chain_files S3_topic_85_to_hg38.bed  S3_topic_85_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S4.topic_106.bed $chain_files S4_topic_106_to_hg38.bed  S4_topic_106_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S5.topic_149.bed $chain_files S5_topic_149_to_hg38.bed  S5_topic_149_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S6.topic_107.bed $chain_files S6_topic_107_to_hg38.bed  S6_topic_107_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S7.topic_134.bed $chain_files S7_topic_134_to_hg38.bed  S7_topic_134_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S8.topic_31.bed $chain_files S8_topic_31_to_hg38.bed  S8_topic_31_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S9.topic_167.bed $chain_files S9_topic_167_to_hg38.bed  S9_topic_167_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S10.topic_194.bed $chain_files S10_topic_194_to_hg38.bed  S10_topic_194_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S11.topic_64.bed $chain_files S11_topic_64_to_hg38.bed  S11_topic_64_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S12.topic_126.bed $chain_files S12_topic_126_to_hg38.bed  S12_topic_126_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S13.topic_82.bed $chain_files S13_topic_82_to_hg38.bed  S13_topic_82_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S14.topic_96.bed $chain_files S14_topic_96_to_hg38.bed  S14_topic_96_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S15.topic_157.bed $chain_files S15_topic_157_to_hg38.bed  S15_topic_157_to_hg38.unmap.bed  -minMatch=0.5
/home/sunkeyong/softeware/LiftOver/liftOver  ./S16.topic_138.bed $chain_files S16_topic_138_to_hg38.bed  S16_topic_138_to_hg38.unmap.bed  -minMatch=0.5
# 
/home/sunkeyong/softeware/LiftOver/liftOver  ./cat.unique.bed $chain_files CatUnique_to_hg38.bed  CatUnique_to_hg38.unmap.bed  -minMatch=0.5
#
#
grep -v '^[^[:space:]]*_' CatUnique_to_hg38.bed | sort -k1,1 -k2,2n > CatUnique_to_hg38.Final.bed
#
grep -v '^[^[:space:]]*_' S1_topic_62_to_hg38.bed | sort -k1,1 -k2,2n > S1_topic_62_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S2_topic_37_to_hg38.bed | sort -k1,1 -k2,2n > S2_topic_37_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S3_topic_85_to_hg38.bed | sort -k1,1 -k2,2n > S3_topic_85_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S4_topic_106_to_hg38.bed | sort -k1,1 -k2,2n > S4_topic_106_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S5_topic_149_to_hg38.bed | sort -k1,1 -k2,2n > S5_topic_149_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S6_topic_107_to_hg38.bed | sort -k1,1 -k2,2n > S6_topic_107_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S7_topic_134_to_hg38.bed | sort -k1,1 -k2,2n > S7_topic_134_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S8_topic_31_to_hg38.bed | sort -k1,1 -k2,2n > S8_topic_31_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S9_topic_167_to_hg38.bed | sort -k1,1 -k2,2n > S9_topic_167_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S10_topic_194_to_hg38.bed | sort -k1,1 -k2,2n > S10_topic_194_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S11_topic_64_to_hg38.bed | sort -k1,1 -k2,2n > S11_topic_64_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S12_topic_126_to_hg38.bed | sort -k1,1 -k2,2n > S12_topic_126_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S13_topic_82_to_hg38.bed | sort -k1,1 -k2,2n > S13_topic_82_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S14_topic_96_to_hg38.bed | sort -k1,1 -k2,2n > S14_topic_96_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S15_topic_157_to_hg38.bed | sort -k1,1 -k2,2n > S15_topic_157_to_hg38.Final.bed
grep -v '^[^[:space:]]*_' S16_topic_138_to_hg38.bed | sort -k1,1 -k2,2n > S16_topic_138_to_hg38.Final.bed
