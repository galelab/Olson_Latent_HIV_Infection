# !/bin/bash
# generate summary of rRNA filtering
              
array=(001_SG_Gale_HIV_Lat_Jurkat_4h_RNA_mRNA_merged	002_SG_Gale_HIV_Lat_Jurkat_4h_RNA_mRNA_merged	003_SG_Gale_HIV_Lat_Jurkat_4h_RNA_mRNA_run01_S3	004_SG_Gale_HIV_Lat_IFNB_Jurkat_4h_RNA_mRNA_merged	005_SG_Gale_HIV_Lat_IFNB_Jurkat_4h_RNA_mRNA_merged	006_SG_Gale_HIV_Lat_IFNB_Jurkat_4h_RNA_mRNA_merged	007_SG_Gale_HIV_Lat_IFNB_Jurkat_8h_RNA_mRNA_merged	008_SG_Gale_HIV_Lat_IFNB_Jurkat_8h_RNA_mRNA_merged	009_SG_Gale_HIV_Lat_IFNB_Jurkat_8h_RNA_mRNA_merged	010_SG_Gale_HIV_Lat_IFNB_Jurkat_12h_RNA_mRNA_merged	011_SG_Gale_HIV_Lat_IFNB_Jurkat_12h_RNA_mRNA_merged	012_SG_Gale_HIV_Lat_IFNB_Jurkat_12h_RNA_mRNA_merged	013_SG_Gale_HIV_Lat_IFNy_Jurkat_4h_RNA_mRNA_run01_S13	014_SG_Gale_HIV_Lat_IFNy_Jurkat_4h_RNA_mRNA_merged	015_SG_Gale_HIV_Lat_IFNy_Jurkat_4h_RNA_mRNA_run01_S15	016_SG_Gale_HIV_Lat_IFNy_Jurkat_8h_RNA_mRNA_run01_S16	017_SG_Gale_HIV_Lat_IFNy_Jurkat_8h_RNA_mRNA_run01_S17	018_SG_Gale_HIV_Lat_IFNy_Jurkat_8h_RNA_mRNA_merged	019_SG_Gale_HIV_Lat_IFNy_Jurkat_12h_RNA_mRNA_merged	020_SG_Gale_HIV_Lat_IFNy_Jurkat_12h_RNA_mRNA_run02_S2	021_SG_Gale_HIV_Lat_IFNy_Jurkat_12h_RNA_mRNA_run02_S3	022_SG_Gale_HIV_Lat_JLat92_4h_RNA_mRNA_merged	023_SG_Gale_HIV_Lat_JLat92_4h_RNA_mRNA_merged	024_SG_Gale_HIV_Lat_JLat92_4h_RNA_mRNA_merged	025_SG_Gale_HIV_Lat_IFNB_JLat92_4h_RNA_mRNA_merged	026_SG_Gale_HIV_Lat_IFNB_JLat92_4h_RNA_mRNA_merged	027_SG_Gale_HIV_Lat_IFNB_JLat92_4h_RNA_mRNA_run02_S9	028_SG_Gale_HIV_Lat_IFNB_JLat92_8h_RNA_mRNA_merged	029_SG_Gale_HIV_Lat_IFNB_JLat92_8h_RNA_mRNA_run02_S11	030_SG_Gale_HIV_Lat_IFNB_JLat92_8h_RNA_mRNA_run02_S12	031_SG_Gale_HIV_Lat_IFNB_JLat92_12h_RNA_mRNA_run02_S13	032_SG_Gale_HIV_Lat_IFNB_JLat92_12h_RNA_mRNA_merged	033_SG_Gale_HIV_Lat_IFNB_JLat92_12h_RNA_mRNA_run02_S15	034_SG_Gale_HIV_Lat_IFNy_JLat92_4h_RNA_mRNA_merged	035_SG_Gale_HIV_Lat_IFNy_JLat92_4h_RNA_mRNA_merged	036_SG_Gale_HIV_Lat_IFNy_JLat92_4h_RNA_mRNA_merged	037_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S1	038_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S2	039_SG_Gale_HIV_Lat_IFNy_JLat92_8h_RNA_mRNA_run03_S3	040_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S4	041_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S5	042_SG_Gale_HIV_Lat_IFNy_JLat92_12h_RNA_mRNA_run03_S6)

cd ../logs/

echo -e 'file\t# read pairs\t# concordant once\t% concordant once\t# concordant more than once\t% concordant more than once\t# of aligned discordant\t% of aligned disdorant\t# of mates condordant once\t% of mates concordant once\t# of mates concordant more than once\t% of mates concordant more than once\toverall % alignment' > rRNA_summary.txt

echo "Array size: ${#array[*]}"
echo "Array items:"

for item in ${array[*]}
do
        printf "   %s\n" $item

	item="$item"_rRNAfilter
        OUTPUT1=$(cat rRNAfilter/$item.log | grep -E 'reads' | sed 's/ reads.*//')
        OUTPUT2=$(cat rRNAfilter/$item.log | grep -E 'concordantly exactly 1 time' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT3=$(cat rRNAfilter/$item.log | grep -E 'concordantly exactly 1 time' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT4=$(cat rRNAfilter/$item.log | grep -E 'aligned concordantly >1 times' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT5=$(cat rRNAfilter/$item.log | grep -E 'aligned concordantly >1 times' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT6=$(cat rRNAfilter/$item.log | grep -E 'aligned discordantly 1 time' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT7=$(cat rRNAfilter/$item.log | grep -E 'aligned discordantly 1 time' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT8=$(cat rRNAfilter/$item.log | grep -E 'aligned exactly 1 time' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT9=$(cat rRNAfilter/$item.log | grep -E 'aligned exactly 1 time' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT10=$(cat rRNAfilter/$item.log | grep -E 'aligned >1 times' | sed 's/ (.*//' | awk '{$1=$1}{ print }')
        OUTPUT11=$(cat rRNAfilter/$item.log | grep -E 'aligned >1 times' | sed 's/.*(//' | sed 's/%).*//')
        OUTPUT12=$(cat rRNAfilter/$item.log | grep -E 'overall alignment rate' | sed 's/%.*//')

        echo -e ''$item'\t'${OUTPUT1}'\t'${OUTPUT2}'\t'${OUTPUT3}'\t'${OUTPUT4}'\t'${OUTPUT5}'\t'${OUTPUT6}'\t'${OUTPUT7}'\t'${OUTPUT8}'\t'${OUTPUT9}'\t'${OUTPUT10}'\t'${OUTPUT11}'\t'${OUTPUT12}'' >> rRNA_summary.txt
         
done
