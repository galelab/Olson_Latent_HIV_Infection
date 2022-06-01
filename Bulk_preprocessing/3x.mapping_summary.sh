#!/bin/bash                                                                                                                                          
# align reads to the genome

cd ../logs/

array=(001_SG_Gale_HIV_Lat_Jurkat_4h_RNA_mRNA_merged	002_SG_Gale_HIV_Lat_Jurkat_4h_RNA_mRNA_merged	003_SG_Gale_HIV_Lat_Jurkat_4h_RNA_mRNA_run01_S3	004_SG_Gale_HIV_Lat_IFNB_Jurkat_4h_RNA_mRNA_merged	005_SG_Gale_HIV_Lat_IFNB_Jurkat_4h_RNA_mRNA_merged	006_SG_Gale_HIV_Lat_IFNB_Jurkat_4h_RNA_mRNA_merged	007_SG_Gale_HIV_Lat_IFNB_Jurkat_8h_RNA_mRNA_merged	008_SG_Gale_HIV_Lat_IFNB_Jurkat_8h_RNA_mRNA_merged	009_SG_Gale_HIV_Lat_IFNB_Jurkat_8h_RNA_mRNA_merged	010_SG_Gale_HIV_Lat_IFNB_Jurkat_12h_RNA_mRNA_merged	011_SG_Gale_HIV_Lat_IFNB_Jurkat_12h_RNA_mRNA_merged	012_SG_Gale_HIV_Lat_IFNB_Jurkat_12h_RNA_mRNA_merged	013_SG_Gale_HIV_Lat_IFNy_Jurkat_4h_RNA_mRNA_run01_S13	014_SG_Gale_HIV_Lat_IFNy_Jurkat_4h_RNA_mRNA_merged	015_SG_Gale_HIV_Lat_IFNy_Jurkat_4h_RNA_mRNA_run01_S15	016_SG_Gale_HIV_Lat_IFNy_Jurkat_8h_RNA_mRNA_run01_S16	017_SG_Gale_HIV_Lat_IFNy_Jurkat_8h_RNA_mRNA_run01_S17	018_SG_Gale_HIV_Lat_IFNy_Jurkat_8h_RNA_mRNA_merged	019_SG_Gale_HIV_Lat_IFNy_Jurkat_12h_RNA_mRNA_merged	020_SG_Gale_HIV_Lat_IFNy_Jurkat_12h_RNA_mRNA_run02_S2	021_SG_Gale_HIV_Lat_IFNy_Jurkat_12h_RNA_mRNA_run02_S3	022_SG_Gale_HIV_Lat_JLat92_4h_RNA_mRNA_merged	023_SG_Gale_HIV_Lat_JLat92_4h_RNA_mRNA_merged	024_SG_Gale_HIV_Lat_JLat92_4h_RNA_mRNA_merged	025_SG_Gale_HIV_Lat_IFNB_JLat92_4h_RNA_mRNA_merged	026_SG_Gale_HIV_Lat_IFNB_JLat92_4h_RNA_mRNA_merged	027_SG_Gale_HIV_Lat_IFNB_JLat92_4h_RNA_mRNA_run02_S9	028_SG_Gale_HIV_Lat_IFNB_JLat92_8h_RNA_mRNA_merged	029_SG_Gale_HIV_Lat_IFNB_JLat92_8h_RNA_mRNA_run02_S11	030_SG_Gale_HIV_Lat_IFNB_JLat92_8h_RNA_mRNA_run02_S12	031_SG_Gale_HIV_Lat_IFNB_JLat92_12h_RNA_mRNA_run02_S13	032_SG_Gale_HIV_Lat_IFNB_JLat92_12h_RNA_mRNA_merged	033_SG_Gale_HIV_Lat_IFNB_JLat92_12h_RNA_mRNA_run02_S15	034_SG_Gale_HIV_Lat_IFNy_JLat92_4h_RNA_mRNA_merged	035_SG_Gale_HIV_Lat_IFNy_JLat92_4h_RNA_mRNA_merged	036_SG_Gale_HIV_Lat_IFNy_JLat92_4h_RNA_mRNA_merged	037_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S1	038_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S2	039_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S3	040_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S4	041_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S5	042_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S6)

echo -e 'File\t# reads\t% uniquely mapped\t# of splices\t% mismatch rate per base\t% deletion rate per base\t% insertion rate per base\t% reads mapped to multiple loci\t% reads unmapped' > mapping_to_human_summary.txt

echo "Array size: ${#array[*]}"
echo "Array items:"

for item in ${array[*]}
do

    printf "   %s\n" $item

    # write log                                                                                                                                       
    cd ../mapping/
    item="$item"Log
    OUTPUT1=$(cat $item.final.out | grep "Number of input reads" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT2=$(cat $item.final.out | grep "Uniquely mapped reads %" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT3=$(cat $item.final.out | grep "Number of splices: Total" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT4=$(cat $item.final.out | grep "Mismatch rate per base, %" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT5=$(cat $item.final.out | grep "Deletion rate per base" | sed 's/.*|//' | awk '{$1=$1}{ print }')
    OUTPUT6=$(cat $item.final.out | grep "Insertion rate per base" | sed 's/.*|//' | awk '{$1=$1}{ print }')

    OUTPUT7=$(cat $item.final.out | grep "% of reads mapped to multiple loci" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')
    OUTPUT8=$(cat $item.final.out | grep "% of reads mapped to too many loci" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')

    OUTPUT9=$(cat $item.final.out | grep "% of reads unmapped: too many mismatches" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')
    OUTPUT10=$(cat $item.final.out | grep "% of reads unmapped: too short" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')
    OUTPUT11=$(cat $item.final.out | grep "% of reads unmapped: other" | sed 's/.*|//' | sed 's/%//' | awk '{$1=$1}{ print }')

    OUTPUT12=$(echo $OUTPUT7 + $OUTPUT8 | bc)
    OUTPUT13=$(echo $OUTPUT9 + $OUTPUT10 + $OUTPUT11 | bc)

    echo -e ''$item'\t'${OUTPUT1}'\t'${OUTPUT2}'\t'${OUTPUT3}'\t'${OUTPUT4}'\t'${OUTPUT5}'\t'${OUTPUT6}'\t'${OUTPUT12}'%\t'${OUTPUT13}'%\t' >> ../logs/mapping_to_human_summary.txt

    cd ../nohmrRNA/

done
