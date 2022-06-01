# !/bin/bash
# use bowtie2 to filter rRNA by mapping to index of human rRNA

cd ../raw

array=(003_SG_Gale_HIV_Lat_Jurkat_4h_RNA_mRNA_run01_S3	005_SG_Gale_HIV_Lat_IFNB_Jurkat_4h_RNA_mRNA_merged	006_SG_Gale_HIV_Lat_IFNB_Jurkat_4h_RNA_mRNA_merged	007_SG_Gale_HIV_Lat_IFNB_Jurkat_8h_RNA_mRNA_merged	008_SG_Gale_HIV_Lat_IFNB_Jurkat_8h_RNA_mRNA_merged	009_SG_Gale_HIV_Lat_IFNB_Jurkat_8h_RNA_mRNA_merged	010_SG_Gale_HIV_Lat_IFNB_Jurkat_12h_RNA_mRNA_merged	011_SG_Gale_HIV_Lat_IFNB_Jurkat_12h_RNA_mRNA_merged	012_SG_Gale_HIV_Lat_IFNB_Jurkat_12h_RNA_mRNA_merged	013_SG_Gale_HIV_Lat_IFNy_Jurkat_4h_RNA_mRNA_run01_S13	014_SG_Gale_HIV_Lat_IFNy_Jurkat_4h_RNA_mRNA_merged	015_SG_Gale_HIV_Lat_IFNy_Jurkat_4h_RNA_mRNA_run01_S15	016_SG_Gale_HIV_Lat_IFNy_Jurkat_8h_RNA_mRNA_run01_S16	017_SG_Gale_HIV_Lat_IFNy_Jurkat_8h_RNA_mRNA_run01_S17	018_SG_Gale_HIV_Lat_IFNy_Jurkat_8h_RNA_mRNA_merged	019_SG_Gale_HIV_Lat_IFNy_Jurkat_12h_RNA_mRNA_merged	020_SG_Gale_HIV_Lat_IFNy_Jurkat_12h_RNA_mRNA_run02_S2	021_SG_Gale_HIV_Lat_IFNy_Jurkat_12h_RNA_mRNA_run02_S3	022_SG_Gale_HIV_Lat_JLat92_4h_RNA_mRNA_merged	023_SG_Gale_HIV_Lat_JLat92_4h_RNA_mRNA_merged	024_SG_Gale_HIV_Lat_JLat92_4h_RNA_mRNA_merged	025_SG_Gale_HIV_Lat_IFNB_JLat92_4h_RNA_mRNA_merged	026_SG_Gale_HIV_Lat_IFNB_JLat92_4h_RNA_mRNA_merged	027_SG_Gale_HIV_Lat_IFNB_JLat92_4h_RNA_mRNA_run02_S9	028_SG_Gale_HIV_Lat_IFNB_JLat92_8h_RNA_mRNA_merged	029_SG_Gale_HIV_Lat_IFNB_JLat92_8h_RNA_mRNA_run02_S11	030_SG_Gale_HIV_Lat_IFNB_JLat92_8h_RNA_mRNA_run02_S12	031_SG_Gale_HIV_Lat_IFNB_JLat92_12h_RNA_mRNA_run02_S13	032_SG_Gale_HIV_Lat_IFNB_JLat92_12h_RNA_mRNA_merged	033_SG_Gale_HIV_Lat_IFNB_JLat92_12h_RNA_mRNA_run02_S15	034_SG_Gale_HIV_Lat_IFNy_JLat92_4h_RNA_mRNA_merged	035_SG_Gale_HIV_Lat_IFNy_JLat92_4h_RNA_mRNA_merged	036_SG_Gale_HIV_Lat_IFNy_JLat92_4h_RNA_mRNA_merged	037_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S1	038_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S2	039_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S3	040_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S4	041_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S5	042_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S6)


echo "Array size: ${#array[*]}"
echo "Array items:"

for item in ${array[*]}
do
        printf "   %s\n" $item

        #bowtie2                                                                                                     
        cmd='/usr/local/bin/bowtie2 -p 50 -x /vol01/genome/rRNA/bowtie2_index/hmrRNA --un-conc-gz ../nohmrRNA/'$item'_nohmrRNA.fastq.gz -1 ../raw/'$item'_R1_001.fastq.gz -2 ../raw/'$item'_R2_001.fastq.gz -S '$item'.sam 1>>../logs/rRNAfilter/'$item'_rRNAfilter.log 2>&1'
        
	echo $cmd
        echo $cmd > ../logs/rRNAfilter/$item\_rRNAfilter.log
        eval $cmd

done
